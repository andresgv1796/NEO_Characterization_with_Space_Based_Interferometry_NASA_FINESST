function fam = gmos_continue_fixT_family( ...
    z0, z1, nSteps, ds, N, mu, Ttarget, ode_opts, opts)
%GMOS_CONTINUE_FIXT_FAMILY Fixed-T pseudo arclength continuation of GMOS tori.
%
% Paper-faithful continuation wrapper:
%   At each successful continuation step, the previously converged family
%   member is used as the reference torus for the phase conditions.
%
% Inputs
%   z0, z1   : two converged solutions to initialize the secant direction
%   nSteps   : number of continuation steps after z1
%   ds       : pseudo arclength step size
%   N, mu    : discretization size and CR3BP parameter
%   Ttarget  : fixed stroboscopic time
%   ode_opts : ODE options passed to the flow / STM propagation
%   opts     : Newton options structure
%
% Output
%   fam      : struct with fields
%                z        : converged family members
%                rho      : rho history
%                T        : T history
%                normInf  : final infinity norm at each accepted step
%                info     : Newton info for each accepted step
%                success  : logical convergence flags

    if nargin < 8 || isempty(ode_opts)
        ode_opts = [];
    end
    if nargin < 9 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'tolInf')
        opts.tolInf = 1e-10;
    end
    if ~isfield(opts, 'verbose')
        opts.verbose = true;
    end

    nStore = nSteps + 2;

    fam = struct();
    fam.z = cell(nStore, 1);
    fam.rho = nan(nStore, 1);
    fam.T = nan(nStore, 1);
    fam.normInf = nan(nStore, 1);
    fam.info = cell(nStore, 1);
    fam.success = false(nStore, 1);

    % Store the two seed solutions
    fam.z{1} = z0;
    fam.z{2} = z1;

    [~, T0, rho0] = gmos_unpack_z(z0, N);
    [~, T1, rho1] = gmos_unpack_z(z1, N);

    fam.T(1) = T0;
    fam.T(2) = T1;
    fam.rho(1) = rho0;
    fam.rho(2) = rho1;

    % Assume supplied seeds are already converged
    fam.success(1) = true;
    fam.success(2) = true;

    lastAccepted = 2;

    for k = 3:nStore
        z_prevprev = fam.z{k-2};
        z_prev = fam.z{k-1};

        % Paper-faithful reference: previous converged family member
        [Xref_k, Tref_k, rhoref_k] = gmos_unpack_z(z_prev, N);

        % Predictor
        z_pred = gmos_predictor_pseudo_arclength(z_prev, z_prevprev, ds);

        % Corrector
        [z_try, info] = gmos_newton_fixT_palc_analytic( ...
            z_pred, N, mu, Xref_k, Tref_k, rhoref_k, Ttarget, ...
            z_prev, z_prevprev, ds, ode_opts, opts);

        if isfield(info, 'normInf') && ~isempty(info.normInf)
            finalNormInf = info.normInf(end);
        else
            finalNormInf = NaN;
        end

        if ~(isfinite(finalNormInf) && finalNormInf < opts.tolInf)
            warning(['GMOS continuation step %d failed: final ||F||inf = %.3e. ', ...
                     'Stopping family here. Reduce ds and restart from the last ', ...
                     'two accepted members.'], k-2, finalNormInf);
            break
        end

        % Accept the new member only after convergence
        fam.z{k} = z_try;
        fam.info{k} = info;
        fam.normInf(k) = finalNormInf;
        fam.success(k) = true;

        [~, T_new, rho_new] = gmos_unpack_z(z_try, N);
        fam.T(k) = T_new;
        fam.rho(k) = rho_new;

        fprintf(['Continuation step %d/%d: T = %.15f, rho = %.15f, ', ...
                 'final ||F||inf = %.3e\n'], ...
            k-2, nSteps, T_new, rho_new, finalNormInf);

        lastAccepted = k;
    end

    % Trim unused tail
    fam.z = fam.z(1:lastAccepted);
    fam.rho = fam.rho(1:lastAccepted);
    fam.T = fam.T(1:lastAccepted);
    fam.normInf = fam.normInf(1:lastAccepted);
    fam.info = fam.info(1:lastAccepted);
    fam.success = fam.success(1:lastAccepted);
end