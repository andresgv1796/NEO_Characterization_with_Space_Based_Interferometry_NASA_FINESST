function branch = gmos_continue_fixT_mu_palc( ...
    y0_sol, y1_sol, nSteps, ds, N, Ttarget, ode_opts, opts)
%GMOS_CONTINUE_FIXT_MU_PALC Pseudo arclength continuation in the extended space [z; mu].
%
% Inputs
%   y0_sol   : first converged extended solution [X0(:); T; rho; mu]
%   y1_sol   : second converged extended solution [X0(:); T; rho; mu]
%   nSteps   : number of continuation steps after y1_sol
%   ds       : pseudo arclength step size in the extended space
%   N        : number of curve nodes
%   Ttarget  : fixed stroboscopic time
%   ode_opts : ODE options
%   opts     : Newton / FD options passed to gmos_newton_fixT_mu_palc_analytic
%
% Output
%   branch   : struct with fields y, z, mu, T, rho, info, success, EGMOS

    if nargin < 7
        ode_opts = [];
    end
    if nargin < 8 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'tolInf');       opts.tolInf = 1e-10; end
    if ~isfield(opts, 'verbose');      opts.verbose = true; end
    if ~isfield(opts, 'computeEGMOS'); opts.computeEGMOS = true; end

    nKeep = nSteps + 2;

    branch = struct();
    branch.y = cell(nKeep, 1);
    branch.z = cell(nKeep, 1);
    branch.mu = nan(nKeep, 1);
    branch.T = nan(nKeep, 1);
    branch.rho = nan(nKeep, 1);
    branch.info = cell(nKeep, 1);
    branch.success = false(nKeep, 1);
    if opts.computeEGMOS
        branch.EGMOS = nan(nKeep, 1);
    else
        branch.EGMOS = [];
    end

    % Store the two initial corrected members.
    branch.y{1} = y0_sol;
    [X0_0, T0, rho0, mu0] = gmos_unpack_y_mu(y0_sol, N);
    branch.z{1} = gmos_pack_z(X0_0, T0, rho0);
    branch.mu(1) = mu0;
    branch.T(1) = T0;
    branch.rho(1) = rho0;
    branch.info{1} = struct('message', 'User supplied first converged extended solution.');
    branch.success(1) = true;
    if opts.computeEGMOS
        [branch.EGMOS(1), ~] = gmos_accuracy_test_paper(branch.z{1}, N, mu0, ode_opts);
    end

    branch.y{2} = y1_sol;
    [X0_1, T1, rho1, mu1] = gmos_unpack_y_mu(y1_sol, N);
    branch.z{2} = gmos_pack_z(X0_1, T1, rho1);
    branch.mu(2) = mu1;
    branch.T(2) = T1;
    branch.rho(2) = rho1;
    branch.info{2} = struct('message', 'User supplied second converged extended solution.');
    branch.success(2) = true;
    if opts.computeEGMOS
        [branch.EGMOS(2), ~] = gmos_accuracy_test_paper(branch.z{2}, N, mu1, ode_opts);
    end

    y_prevprev = y0_sol;
    y_prev = y1_sol;
    lastAccepted = 2;

    for k = 3:nKeep
        y_pred = gmos_predictor_pseudo_arclength(y_prev, y_prevprev, ds);

        [Xref, Tref, rhoref, muRef] = gmos_unpack_y_mu(y_prev, N);

        [y_corr, info] = gmos_newton_fixT_mu_palc_analytic( ...
            y_pred, N, Xref, Tref, rhoref, muRef, Ttarget, ...
            y_prev, y_prevprev, ds, ode_opts, opts);

        finalNormInf = info.normInf(end);
        if ~(isfinite(finalNormInf) && finalNormInf < opts.tolInf)
            warning(['Extended PALC failed at continuation step %d: ', ...
                     'final ||F||inf = %.3e.'], k - 2, finalNormInf);
            break
        end

        branch.y{k} = y_corr;
        [X0k, Tk, rhok, muk] = gmos_unpack_y_mu(y_corr, N);
        branch.z{k} = gmos_pack_z(X0k, Tk, rhok);
        branch.mu(k) = muk;
        branch.T(k) = Tk;
        branch.rho(k) = rhok;
        branch.info{k} = info;
        branch.success(k) = true;
        if opts.computeEGMOS
            [branch.EGMOS(k), ~] = gmos_accuracy_test_paper(branch.z{k}, N, muk, ode_opts);
        end

        if opts.verbose
            if opts.computeEGMOS
                fprintf(['PALC step %d/%d: mu = %.16g, T = %.15f, rho = %.15f, ', ...
                         'EGMOS = %.3e, final ||F||inf = %.3e\n'], ...
                        k - 2, nSteps, muk, Tk, rhok, branch.EGMOS(k), finalNormInf);
            else
                fprintf(['PALC step %d/%d: mu = %.16g, T = %.15f, rho = %.15f, ', ...
                         'final ||F||inf = %.3e\n'], ...
                        k - 2, nSteps, muk, Tk, rhok, finalNormInf);
            end
        end

        y_prevprev = y_prev;
        y_prev = y_corr;
        lastAccepted = k;
    end

    branch.y = branch.y(1:lastAccepted);
    branch.z = branch.z(1:lastAccepted);
    branch.mu = branch.mu(1:lastAccepted);
    branch.T = branch.T(1:lastAccepted);
    branch.rho = branch.rho(1:lastAccepted);
    branch.info = branch.info(1:lastAccepted);
    branch.success = branch.success(1:lastAccepted);
    if opts.computeEGMOS
        branch.EGMOS = branch.EGMOS(1:lastAccepted);
    end
end
