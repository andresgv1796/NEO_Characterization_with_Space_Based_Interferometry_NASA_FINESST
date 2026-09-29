function branch = gmos_continue_fixT_mu_simple( ...
    z_start, muVals, N, Ttarget, ode_opts, opts)
%GMOS_CONTINUE_FIXT_MU_SIMPLE Simple parameter stepping in mu for one fixed-T torus branch.
%
% The torus is corrected independently at each mu_k using the previously
% converged torus as both:
%   1) initial guess, and
%   2) reference torus for the phase conditions.
%
% IMPORTANT
%   z_start is assumed to be a converged fixed-T torus at muVals(1).
%   The first branch point is stored directly without re-solving.
%
% Inputs
%   z_start  : converged torus at muVals(1), packed as [X0(:); T; rho]
%   muVals   : vector of mu values to march through, in desired order
%   N        : number of invariant-curve nodes
%   Ttarget  : fixed stroboscopic time used in gmos_newton_fixT_analytic
%   ode_opts : ODE options for flow / STM propagation
%   opts     : options struct, same style as gmos_newton_fixT_analytic
%              Additional optional fields:
%                computeEGMOS (default true)
%                verboseStep  (default true)
%
% Output
%   branch   : struct with fields
%                z        : converged branch members
%                mu       : mu values actually accepted
%                rho      : rho history
%                T        : T history
%                normInf  : final Newton infinity norm at each accepted step
%                info     : Newton info for each accepted step
%                success  : logical convergence flags
%                EGMOS    : midpoint accuracy metric (if requested)
%
% Notes
%   - This is not PALC in extended (z,mu) space.
%   - It is meant for smooth branches and moderate mu steps.
%   - If this fails, reduce the mu step or switch to an augmented PALC solver.

    if nargin < 5 || isempty(ode_opts)
        ode_opts = [];
    end
    if nargin < 6 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'tolInf');       opts.tolInf = 1e-10; end
    if ~isfield(opts, 'verbose');      opts.verbose = true; end
    if ~isfield(opts, 'computeEGMOS'); opts.computeEGMOS = true; end
    if ~isfield(opts, 'verboseStep');  opts.verboseStep = true; end

    muVals = muVals(:);
    nMu = numel(muVals);
    if nMu < 1
        error('muVals must contain at least one value.');
    end

    branch = struct();
    branch.z = cell(nMu, 1);
    branch.mu = nan(nMu, 1);
    branch.rho = nan(nMu, 1);
    branch.T = nan(nMu, 1);
    branch.normInf = nan(nMu, 1);
    branch.info = cell(nMu, 1);
    branch.success = false(nMu, 1);
    if opts.computeEGMOS
        branch.EGMOS = nan(nMu, 1);
    else
        branch.EGMOS = [];
    end

    % Store the initial converged torus.
    branch.z{1} = z_start;
    branch.mu(1) = muVals(1);
    [~, T1, rho1] = gmos_unpack_z(z_start, N);
    branch.T(1) = T1;
    branch.rho(1) = rho1;
    branch.normInf(1) = 0.0;
    branch.info{1} = struct('normInf', 0.0, 'norm2', 0.0, ...
                            'message', 'Initial converged torus supplied by user.');
    branch.success(1) = true;
    if opts.computeEGMOS
        [branch.EGMOS(1), ~] = gmos_accuracy_test_paper(z_start, N, muVals(1), ode_opts);
    end

    z_guess = z_start;
    lastAccepted = 1;

    for k = 2:nMu
        mu_k = muVals(k);

        % Use previous converged torus as phase reference.
        [Xref_k, Tref_k, rhoref_k] = gmos_unpack_z(branch.z{k-1}, N);

        [z_try, info] = gmos_newton_fixT_analytic( ...
            z_guess, N, mu_k, Xref_k, Tref_k, rhoref_k, Ttarget, ode_opts, opts);

        if isfield(info, 'normInf') && ~isempty(info.normInf)
            finalNormInf = info.normInf(end);
        else
            finalNormInf = NaN;
        end

        if ~(isfinite(finalNormInf) && finalNormInf < opts.tolInf)
            warning(['GMOS mu-continuation failed at k = %d (mu = %.16g): ', ...
                     'final ||F||inf = %.3e. Reduce the mu step and restart ', ...
                     'from the last accepted torus.'], ...
                     k, mu_k, finalNormInf);
            break
        end

        branch.z{k} = z_try;
        branch.mu(k) = mu_k;
        branch.info{k} = info;
        branch.normInf(k) = finalNormInf;
        branch.success(k) = true;

        [~, T_k, rho_k] = gmos_unpack_z(z_try, N);
        branch.T(k) = T_k;
        branch.rho(k) = rho_k;

        if opts.computeEGMOS
            [branch.EGMOS(k), ~] = gmos_accuracy_test_paper(z_try, N, mu_k, ode_opts);
        end

        if opts.verboseStep
            if opts.computeEGMOS
                fprintf(['mu-step %d/%d: mu = %.16g, T = %.15f, rho = %.15f, ', ...
                         'EGMOS = %.3e, final ||F||inf = %.3e\n'], ...
                        k-1, nMu-1, mu_k, T_k, rho_k, branch.EGMOS(k), finalNormInf);
            else
                fprintf(['mu-step %d/%d: mu = %.16g, T = %.15f, rho = %.15f, ', ...
                         'final ||F||inf = %.3e\n'], ...
                        k-1, nMu-1, mu_k, T_k, rho_k, finalNormInf);
            end
        end

        z_guess = z_try;
        lastAccepted = k;
    end

    % Trim unused tail.
    branch.z = branch.z(1:lastAccepted);
    branch.mu = branch.mu(1:lastAccepted);
    branch.rho = branch.rho(1:lastAccepted);
    branch.T = branch.T(1:lastAccepted);
    branch.normInf = branch.normInf(1:lastAccepted);
    branch.info = branch.info(1:lastAccepted);
    branch.success = branch.success(1:lastAccepted);
    if opts.computeEGMOS
        branch.EGMOS = branch.EGMOS(1:lastAccepted);
    end
end
