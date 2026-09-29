function qpo = generate_gmos_qpo_torus_from_po(po, sys, cfg, odeOpts)
%GENERATE_GMOS_QPO_TORUS_FROM_PO Build and correct one GMOS QPO torus.
%
% The output qpo.XCurve is a 6xN corrected invariant curve. Three
% collectors can then be selected from its columns.
%
% cfg fields:
%   N                 number of invariant-curve samples
%   mode              'planar' or 'vertical'
%   initialAmplitude  small center-mode amplitude in nondimensional units
%   corrector         'fixT', 'fixC', or 'none'
%   maxIter           Newton iterations
%   tolInf            Newton infinity-norm tolerance
%   verbose           print Newton progress
%   plotInitialGuess  true/false
%   CtargetMode       for fixC: 'parent' or 'initialMean'

    if nargin < 4
        odeOpts = [];
    end
    if nargin < 3 || isempty(cfg)
        cfg = struct();
    end

    defaults = struct();
    defaults.N = 51;
    defaults.mode = 'vertical';
    defaults.initialAmplitude = 1e-5;
    defaults.corrector = 'fixT';
    defaults.maxIter = 10;
    defaults.tolInf = 1e-10;
    defaults.verbose = true;
    defaults.plotInitialGuess = false;
    defaults.CtargetMode = 'parent';
    cfg = apply_struct_defaults(cfg, defaults);

    fprintf('\nGenerating GMOS QPO torus from parent periodic orbit\n');
    fprintf('  center mode       = %s\n', cfg.mode);
    fprintf('  N curve samples   = %d\n', cfg.N);
    fprintf('  initial amplitude = %.6e nondim.\n', cfg.initialAmplitude);
    fprintf('  corrector         = %s\n', cfg.corrector);

    cen = gmos_center_eigpair_from_monodromy([0, po.T], po.x0, sys.mu, cfg.mode);
    init = gmos_init_curve_from_center_mode( ...
        po.x0, po.T, cen, cfg.N, cfg.initialAmplitude, cfg.plotInitialGuess);

    z0 = gmos_pack_z(init.X0, init.T0, init.rho0);

    newtonOpts = struct();
    newtonOpts.maxIter = cfg.maxIter;
    newtonOpts.tolInf = cfg.tolInf;
    newtonOpts.verbose = cfg.verbose;

    switch lower(cfg.corrector)
        case 'none'
            z = z0;
            info = struct('normInf', NaN, 'norm2', NaN, 'converged', false, ...
                'iterations', 0, 'note', 'uncorrected initial center-mode curve');

        case {'fixt','fix_t','period','fixedperiod'}
            [z, info] = gmos_newton_fixT_analytic( ...
                z0, cfg.N, sys.mu, init.X0, init.T0, init.rho0, po.T, odeOpts, newtonOpts);

        case {'fixc','fix_c','jacobi','fixedjacobi'}
            switch lower(cfg.CtargetMode)
                case 'parent'
                    Ctarget = Jacobi_Constant(po.x0(:).', sys.mu);
                case 'initialmean'
                    Ctarget = mean(Jacobi_Constant(init.X0.', sys.mu));
                otherwise
                    error('Unknown CtargetMode. Use parent or initialMean.');
            end
            [z, info] = gmos_newton_fixC_analytic( ...
                z0, cfg.N, sys.mu, init.X0, init.T0, init.rho0, Ctarget, odeOpts, newtonOpts);
            info.Ctarget = Ctarget;

        otherwise
            error('Unknown GMOS corrector "%s".', cfg.corrector);
    end

    [XCurve, TQPO, rhoQPO] = gmos_unpack_z(z, cfg.N);
    CCurve = Jacobi_Constant(XCurve.', sys.mu);

    qpo = struct();
    qpo.source = 'generated';
    qpo.z = z;
    qpo.XCurve = XCurve;
    qpo.T = TQPO;
    qpo.rho = rhoQPO;
    qpo.N = cfg.N;
    qpo.mode = cfg.mode;
    qpo.corrector = cfg.corrector;
    qpo.initialAmplitude = cfg.initialAmplitude;
    qpo.centerEig = cen;
    qpo.init = init;
    qpo.info = info;
    qpo.CCurve = CCurve;
    qpo.CMean = mean(CCurve);
    qpo.CMin = min(CCurve);
    qpo.CMax = max(CCurve);
    qpo.parent = po;
    qpo.system = sys;

    fprintf('  corrected T       = %.16e nondim.\n', qpo.T);
    fprintf('  corrected rho     = %.16e rad = %.16e cycles\n', qpo.rho, qpo.rho/(2*pi));
    fprintf('  |T - T_PO|/T_PO   = %.3e\n', abs(qpo.T - po.T)/po.T);
    fprintf('  Jacobi mean/range = %.16e / [%.16e, %.16e]\n', qpo.CMean, qpo.CMin, qpo.CMax);
end
