function fam = gmos_continue_fixC_family_adaptive( ...
    z0, z1, nSteps, ds, N, mu, Ctarget, ode_opts, opts)
%GMOS_CONTINUE_FIXC_FAMILY_ADAPTIVE Fixed-C pseudo arclength continuation of GMOS tori
% with adaptive PALC step size and backtracking.
%
% At each successful continuation step, the previously converged family
% member is used as the reference torus for the phase conditions.
%
% Inputs
%   z0, z1   : two converged solutions to initialize the secant direction
%   nSteps   : number of continuation steps after z1
%   ds       : initial pseudo arclength step size
%   N, mu    : discretization size and CR3BP parameter
%   Ctarget  : fixed Jacobi constant target (average along the curve)
%   ode_opts : ODE options passed to the flow / STM propagation
%   opts     : Newton / continuation options structure
%
% opts fields with defaults
%   tolInf          = 1e-10
%   verbose         = true
%   verboseStep     = true
%   maxBacktrack    = 6
%   backtrackFactor = 0.5
%   growFactor      = 1.15
%   dsMax           = Inf
%   dsMin           = 0
%
% Output
%   fam      : struct with fields
%                z        : converged family members
%                rho      : rho history
%                T        : T history
%                Cavg     : average Jacobi constant history
%                normInf  : final infinity norm at each accepted step
%                info     : Newton info for each accepted step
%                success  : logical convergence flags
%                dsUsed   : accepted PALC step used for each member

    if nargin < 8 || isempty(ode_opts)
        ode_opts = [];
    end
    if nargin < 9 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts, 'tolInf');          opts.tolInf = 1e-10; end
    if ~isfield(opts, 'verbose');         opts.verbose = true; end
    if ~isfield(opts, 'verboseStep');     opts.verboseStep = true; end
    if ~isfield(opts, 'maxBacktrack');    opts.maxBacktrack = 6; end
    if ~isfield(opts, 'backtrackFactor'); opts.backtrackFactor = 0.5; end
    if ~isfield(opts, 'growFactor');      opts.growFactor = 1.15; end
    if ~isfield(opts, 'dsMax');           opts.dsMax = Inf; end
    if ~isfield(opts, 'dsMin');           opts.dsMin = 0; end

    if ~isfinite(ds) || ds == 0
        error('Initial ds must be a finite nonzero scalar.');
    end

    nStore = nSteps + 2;

    fam = struct();
    fam.z = cell(nStore, 1);
    fam.rho = nan(nStore, 1);
    fam.T = nan(nStore, 1);
    fam.Cavg = nan(nStore, 1);
    fam.normInf = nan(nStore, 1);
    fam.info = cell(nStore, 1);
    fam.success = false(nStore, 1);
    fam.dsUsed = nan(nStore, 1);

    fam.z{1} = z0;
    fam.z{2} = z1;

    [X0_seed, T0, rho0] = gmos_unpack_z(z0, N);
    [X1_seed, T1, rho1] = gmos_unpack_z(z1, N);

    fam.T(1) = T0;
    fam.T(2) = T1;
    fam.rho(1) = rho0;
    fam.rho(2) = rho1;
    fam.Cavg(1) = mean(Jacobi_Constant(X0_seed.', mu));
    fam.Cavg(2) = mean(Jacobi_Constant(X1_seed.', mu));

    fam.success(1) = true;
    fam.success(2) = true;
    fam.dsUsed(1) = NaN;
    fam.dsUsed(2) = norm(z1 - z0);

    lastAccepted = 2;
    dsCurrent = ds;

    for k = 3:nStore
        z_prevprev = fam.z{k-2};
        z_prev = fam.z{k-1};

        [Xref_k, Tref_k, rhoref_k] = gmos_unpack_z(z_prev, N);

        accept = false;
        dsTry = dsCurrent;
        z_try = [];
        info = struct();

        for retry = 1:opts.maxBacktrack
            if abs(dsTry) <= opts.dsMin
                break
            end

            z_pred = gmos_predictor_pseudo_arclength(z_prev, z_prevprev, dsTry);

            [z_try_local, info_local] = gmos_newton_fixC_palc_analytic( ...
                z_pred, N, mu, Xref_k, Tref_k, rhoref_k, Ctarget, ...
                z_prev, z_prevprev, dsTry, ode_opts, opts);

            if isfield(info_local, 'normInf') && ~isempty(info_local.normInf)
                finalNormInf = info_local.normInf(end);
            else
                finalNormInf = NaN;
            end

            if isfinite(finalNormInf) && finalNormInf < opts.tolInf
                accept = true;
                z_try = z_try_local;
                info = info_local;
                break
            end

            info = info_local;
            dsTry = dsTry * opts.backtrackFactor;
        end

        if ~accept
            if isfield(info, 'normInf') && ~isempty(info.normInf)
                finalNormInf = info.normInf(end);
            else
                finalNormInf = NaN;
            end

            warning(['GMOS fixed-C continuation failed at step %d: ', ...
                     'final ||F||inf = %.3e after backtracking. ', ...
                     'Stopping family here.'], ...
                    k-2, finalNormInf);
            break
        end

        fam.z{k} = z_try;
        fam.info{k} = info;
        fam.normInf(k) = info.normInf(end);
        fam.success(k) = true;
        fam.dsUsed(k) = dsTry;

        [X_new, T_new, rho_new] = gmos_unpack_z(z_try, N);
        fam.T(k) = T_new;
        fam.rho(k) = rho_new;
        fam.Cavg(k) = mean(Jacobi_Constant(X_new.', mu));

        if opts.verboseStep
            fprintf(['Continuation step %d/%d: T = %.15f, rho = %.15f, ', ...
                     'Cavg = %.15f, ds = %.3e, final ||F||inf = %.3e\n'], ...
                k-2, nSteps, T_new, rho_new, fam.Cavg(k), dsTry, info.normInf(end));
        end

        lastAccepted = k;
        dsCurrent = sign(dsTry) * min(abs(opts.growFactor * dsTry), abs(opts.dsMax));
    end

    fam.z = fam.z(1:lastAccepted);
    fam.rho = fam.rho(1:lastAccepted);
    fam.T = fam.T(1:lastAccepted);
    fam.Cavg = fam.Cavg(1:lastAccepted);
    fam.normInf = fam.normInf(1:lastAccepted);
    fam.info = fam.info(1:lastAccepted);
    fam.success = fam.success(1:lastAccepted);
    fam.dsUsed = fam.dsUsed(1:lastAccepted);
end