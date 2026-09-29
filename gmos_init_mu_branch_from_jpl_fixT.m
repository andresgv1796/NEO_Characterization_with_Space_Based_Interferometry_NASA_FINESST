function out = gmos_init_mu_branch_from_jpl_fixT( ...
    catalogFile, orbitIndex, muSeed, mode, N, KSeeds, dMuInit, ...
    nStartupSteps, ode_opts, optsNewton, optsPALC)
%GMOS_INIT_MU_BRANCH_FROM_JPL
% Minimal initializer for fixed-T continuation in mu.
%
% Purpose
%   Build only what is actually needed to start
%   gmos_continue_fixT_mu_palc_weighted and then feed the result to
%   gmos_extend_fixT_mu_branch_weighted_bidirectional.
%
% Required inputs
%   catalogFile    : JPL CSV file for pick_po_from_JPL
%   orbitIndex     : orbit index in the catalog
%   muSeed         : catalog mass ratio
%   mode           : 'planar' or 'vertical'
%   N              : number of invariant-curve nodes
%   KSeeds         : [K0 K1] for two nearby center-mode initializations
%   dMuInit        : initial mu step used to create the second mu seed
%   nStartupSteps  : number of weighted mu-PALC startup steps
%
% Optional inputs
%   ode_opts       : ODE options
%   optsNewton     : Newton options for correction
%   optsPALC       : options for weighted mu-PALC
%
% Output
%   out.branchMaster : initial mu branch, ready for bidirectional extension
%   out.Ttarget      : fixed period used by the mu branch
%   out.Ctarget      : Jacobi constant inherited from the seed PO
%   out.restart      : startup seed data
%   out.seed         : corrected same-mu seed data
%   out.branchSimple : one-step simple mu startup branch
%   out.po           : seed PO data

    if nargin < 9 || isempty(ode_opts)
        ode_opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    end
    if nargin < 10 || isempty(optsNewton)
        optsNewton = struct();
    end
    if nargin < 11 || isempty(optsPALC)
        optsPALC = struct();
    end

    if ~isfield(optsNewton, 'maxIter'); optsNewton.maxIter = 12; end
    if ~isfield(optsNewton, 'tolInf');  optsNewton.tolInf  = 1e-10; end
    if ~isfield(optsNewton, 'verbose'); optsNewton.verbose = true; end

    if ~isfield(optsPALC, 'maxIter');         optsPALC.maxIter = 12; end
    if ~isfield(optsPALC, 'tolInf');          optsPALC.tolInf = 1e-10; end
    if ~isfield(optsPALC, 'verbose');         optsPALC.verbose = true; end
    if ~isfield(optsPALC, 'computeEGMOS');    optsPALC.computeEGMOS = true; end
    if ~isfield(optsPALC, 'fdMu');            optsPALC.fdMu = 1e-8; end
    if ~isfield(optsPALC, 'maxBacktrack');    optsPALC.maxBacktrack = 8; end
    if ~isfield(optsPALC, 'backtrackFactor'); optsPALC.backtrackFactor = 0.5; end
    if ~isfield(optsPALC, 'growFactor');      optsPALC.growFactor = 1.10; end
    if ~isfield(optsPALC, 'enforceMuDirection')
        optsPALC.enforceMuDirection = true;
    end

    validateattributes(KSeeds, {'numeric'}, {'vector','numel',2,'real','finite'});
    validateattributes(dMuInit, {'numeric'}, {'scalar','real','finite','nonzero'});
    validateattributes(nStartupSteps, {'numeric'}, {'scalar','integer','>=',0});

    mode = char(string(mode));
    if ~ismember(lower(mode), {'planar','vertical'})
        error('mode must be ''planar'' or ''vertical''.');
    end

    % ------------------------------------------------------------
    % Seed periodic orbit from JPL
    % ------------------------------------------------------------
    po = pick_po_from_JPL(catalogFile, orbitIndex);
    xPO0 = po.x0(:);
    TPO  = po.T;

    solPO = ode89(@(t,y) CR3BP_mu(t,y,muSeed), [0, TPO], xPO0, ode_opts);
    closureError = norm(solPO.y(:,end) - xPO0);

    Ctarget = Jacobi_Constant(xPO0.', muSeed);   % keep only as metadata/diagnostic
    Ttarget = TPO;                               % this is the actual fixed-T target

    % ------------------------------------------------------------
    % Build two nearby GMOS seeds at the same mu and correct them at fixed T
    % ------------------------------------------------------------
    demo0 = gmos_demo_invariance_only(muSeed, xPO0, TPO, N, KSeeds(1), false, mode);
    demo1 = gmos_demo_invariance_only(muSeed, xPO0, TPO, N, KSeeds(2), false, mode);

    init0 = demo0.init;
    init1 = demo1.init;

    Xref   = init0.X0;
    Tref   = init0.T0;
    rhoref = init0.rho0;

    z0 = gmos_pack_z(init0.X0, init0.T0, init0.rho0);
    z1 = gmos_pack_z(init1.X0, init1.T0, init1.rho0);

    [z0_sol, info0] = gmos_newton_fixT_analytic( ...
    z0, N, muSeed, Xref, Tref, rhoref, Ttarget, ode_opts, optsNewton);

    [Xref1, Tref1, rhoref1] = gmos_unpack_z(z0_sol, N);

    [z1_sol, info1] = gmos_newton_fixT_analytic( ...
    z1, N, muSeed, Xref1, Tref1, rhoref1, Ttarget, ode_opts, optsNewton);

    % ------------------------------------------------------------
    % Choose one corrected torus as the base torus for fixed-T mu continuation
    % ------------------------------------------------------------
    zBase = z1_sol;
    [~, Ttarget, ~] = gmos_unpack_z(zBase, N);

    % ------------------------------------------------------------
    % Make the second mu seed with the simple fixed-T solver
    % ------------------------------------------------------------
    muValsSimple = [muSeed; muSeed + dMuInit];

    branchSimple = gmos_continue_fixT_mu_simple( ...
        zBase, muValsSimple, N, Ttarget, ode_opts, optsNewton);

    if numel(branchSimple.z) < 2
        error(['Simple mu startup produced fewer than two members. ', ...
               'Reduce |dMuInit| and try again.']);
    end

    y0 = [branchSimple.z{1}(:); branchSimple.mu(1)];
    y1 = [branchSimple.z{2}(:); branchSimple.mu(2)];

    % ------------------------------------------------------------
    % Build the weighted PALC step size from the startup pair
    % ------------------------------------------------------------
    scale = gmos_build_palc_scaling_mu(y0, y1, N, struct());

    dySeedScaled = scale.w .* (y1 - y0);
    tHat0 = dySeedScaled / norm(dySeedScaled);

    deltaMuSeed = y1(end) - y0(end);
    dsStartup = (scale.w(end) * deltaMuSeed) / tHat0(end);

    % ------------------------------------------------------------
    % Initial weighted mu-PALC branch
    % ------------------------------------------------------------
    branchMaster = gmos_continue_fixT_mu_palc_weighted( ...
        y0, y1, nStartupSteps, dsStartup, N, Ttarget, ode_opts, optsPALC);

    % ------------------------------------------------------------
    % Package output
    % ------------------------------------------------------------
    out = struct();

    out.po = struct();
    out.po.data = po;
    out.po.xPO0 = xPO0;
    out.po.TPO = TPO;
    out.po.closureError = closureError;

    out.Ctarget = Ctarget;
    out.Ttarget = Ttarget;

    out.seed = struct();
    out.seed.z0_sol = z0_sol;
    out.seed.z1_sol = z1_sol;
    out.seed.info0 = info0;
    out.seed.info1 = info1;

    out.branchSimple = branchSimple;
    out.branchMaster = branchMaster;

    out.restart = struct();
    out.restart.y0 = y0;
    out.restart.y1 = y1;
    out.restart.scale = scale;
    out.restart.deltaMuSeed = deltaMuSeed;
    out.restart.dsStartup = dsStartup;
end