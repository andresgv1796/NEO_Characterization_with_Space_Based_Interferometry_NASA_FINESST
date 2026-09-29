function branch = gmos_continue_fixT_mu_palc_weighted( ...
    y0_sol, y1_sol, nSteps, ds, N, Ttarget, ode_opts, opts)
%GMOS_CONTINUE_FIXT_MU_PALC_WEIGHTED Weighted PALC continuation in the extended space [z; mu].
%
% This is a robust replacement for the unweighted extended PALC routine.
% The pseudo arclength metric is applied in the scaled variable ys = D*y,
% with diagonal D chosen from the initial seed secant unless the user
% overrides the block scales.
%
% Inputs
%   y0_sol, y1_sol : first two converged extended solutions [X0(:); T; rho; mu]
%                    on the same branch, ordered in the desired traversal direction
%   nSteps         : number of continuation steps after y1_sol
%   ds             : PALC step size in the scaled space. A good first choice is
%                      ds0 = norm(scale.w .* (y1_sol - y0_sol));
%                    and then ds = alpha * ds0 with alpha in [0.25, 1].
%   N              : number of invariant-curve nodes
%   Ttarget        : fixed stroboscopic time
%   ode_opts       : ODE options
%   opts           : optional struct. In addition to Newton options, supports
%                      computeEGMOS         (default true)
%                      verboseStep          (default true)
%                      maxBacktrack        (default 6)
%                      backtrackFactor     (default 0.5)
%                      growFactor          (default 1.15)
%                      dsMax               (default Inf)
%                      enforceMuDirection  (default false)
%                      muDirectionTol      (default 0)
%                      xScale, TScale, rhoScale, muScale
%
% Output
%   branch         : struct with fields
%                      y, z, mu, T, rho, info, success, EGMOS,
%                      dsUsed, scale

    if nargin < 7
        ode_opts = [];
    end
    if nargin < 8 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'tolInf');              opts.tolInf = 1e-10; end
    if ~isfield(opts, 'verbose');             opts.verbose = true; end
    if ~isfield(opts, 'verboseStep');         opts.verboseStep = true; end
    if ~isfield(opts, 'computeEGMOS');        opts.computeEGMOS = true; end
    if ~isfield(opts, 'maxBacktrack');        opts.maxBacktrack = 6; end
    if ~isfield(opts, 'backtrackFactor');     opts.backtrackFactor = 0.5; end
    if ~isfield(opts, 'growFactor');          opts.growFactor = 1.15; end
    if ~isfield(opts, 'dsMax');               opts.dsMax = Inf; end
    if ~isfield(opts, 'enforceMuDirection');  opts.enforceMuDirection = false; end
    if ~isfield(opts, 'muDirectionTol');      opts.muDirectionTol = 0.0; end

    y0_sol = y0_sol(:);
    y1_sol = y1_sol(:);

    if numel(y0_sol) ~= 6*N + 3 || numel(y1_sol) ~= 6*N + 3
        error('y0_sol and y1_sol must both have length 6*N + 3.');
    end

    scale = gmos_build_palc_scaling_mu(y0_sol, y1_sol, N, opts);
    w = scale.w;

    dsCurrent = ds;
    if ~isfinite(dsCurrent) || dsCurrent == 0
        error('ds must be a finite nonzero scalar.');
    end

    [~, ~, ~, mu0] = gmos_unpack_y_mu(y0_sol, N);
    [~, ~, ~, mu1] = gmos_unpack_y_mu(y1_sol, N);
    muSeedSign = sign(mu1 - mu0);

    nKeep = nSteps + 2;

    branch = struct();
    branch.y = cell(nKeep, 1);
    branch.z = cell(nKeep, 1);
    branch.mu = nan(nKeep, 1);
    branch.T = nan(nKeep, 1);
    branch.rho = nan(nKeep, 1);
    branch.info = cell(nKeep, 1);
    branch.success = false(nKeep, 1);
    branch.dsUsed = nan(nKeep, 1);
    branch.scale = scale;
    if opts.computeEGMOS
        branch.EGMOS = nan(nKeep, 1);
    else
        branch.EGMOS = [];
    end

    branch.y{1} = y0_sol;
    [X0_0, T0, rho0, mu0] = gmos_unpack_y_mu(y0_sol, N);
    branch.z{1} = gmos_pack_z(X0_0, T0, rho0);
    branch.mu(1) = mu0;
    branch.T(1) = T0;
    branch.rho(1) = rho0;
    branch.info{1} = struct('message', 'User supplied first converged extended solution.');
    branch.success(1) = true;
    branch.dsUsed(1) = NaN;
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
    branch.dsUsed(2) = norm(w .* (y1_sol - y0_sol));
    if opts.computeEGMOS
        [branch.EGMOS(2), ~] = gmos_accuracy_test_paper(branch.z{2}, N, mu1, ode_opts);
    end

    y_prevprev = y0_sol;
    y_prev = y1_sol;
    lastAccepted = 2;

    for k = 3:nKeep
        accept = false;
        dsTry = dsCurrent;
        lastInfo = struct();
        y_corr = [];

        for retry = 1:opts.maxBacktrack
            [y_pred, ~] = gmos_predictor_pseudo_arclength_weighted(y_prev, y_prevprev, dsTry, w);
            [Xref, Tref, rhoref, muRef] = gmos_unpack_y_mu(y_prev, N);

            newtonOpts = opts;
            if isfield(newtonOpts, 'verboseStep')
                newtonOpts = rmfield(newtonOpts, 'verboseStep');
            end
            if isfield(newtonOpts, 'computeEGMOS')
                newtonOpts = rmfield(newtonOpts, 'computeEGMOS');
            end
            if isfield(newtonOpts, 'maxBacktrack')
                newtonOpts = rmfield(newtonOpts, 'maxBacktrack');
            end
            if isfield(newtonOpts, 'backtrackFactor')
                newtonOpts = rmfield(newtonOpts, 'backtrackFactor');
            end
            if isfield(newtonOpts, 'growFactor')
                newtonOpts = rmfield(newtonOpts, 'growFactor');
            end
            if isfield(newtonOpts, 'dsMax')
                newtonOpts = rmfield(newtonOpts, 'dsMax');
            end
            if isfield(newtonOpts, 'enforceMuDirection')
                newtonOpts = rmfield(newtonOpts, 'enforceMuDirection');
            end
            if isfield(newtonOpts, 'muDirectionTol')
                newtonOpts = rmfield(newtonOpts, 'muDirectionTol');
            end
            if isfield(newtonOpts, 'xScale')
                newtonOpts = rmfield(newtonOpts, 'xScale');
            end
            if isfield(newtonOpts, 'TScale')
                newtonOpts = rmfield(newtonOpts, 'TScale');
            end
            if isfield(newtonOpts, 'rhoScale')
                newtonOpts = rmfield(newtonOpts, 'rhoScale');
            end
            if isfield(newtonOpts, 'muScale')
                newtonOpts = rmfield(newtonOpts, 'muScale');
            end

            [y_corr_try, infoTry] = gmos_newton_fixT_mu_palc_weighted( ...
                y_pred, N, Xref, Tref, rhoref, muRef, Ttarget, ...
                y_prev, y_prevprev, dsTry, w, ode_opts, newtonOpts);

            if isfield(infoTry, 'normInf') && ~isempty(infoTry.normInf)
                finalNormInf = infoTry.normInf(end);
            else
                finalNormInf = NaN;
            end

            muPrev = y_prev(end);
            muCorr = y_corr_try(end);
            muStep = muCorr - muPrev;
            dirOK = true;
            if opts.enforceMuDirection && muSeedSign ~= 0
                dirOK = (muStep * muSeedSign) > opts.muDirectionTol;
            end

            if isfinite(finalNormInf) && finalNormInf < opts.tolInf && dirOK
                accept = true;
                y_corr = y_corr_try;
                lastInfo = infoTry;
                break;
            end

            lastInfo = infoTry;
            dsTry = dsTry * opts.backtrackFactor;
        end

        if ~accept
            if isfield(lastInfo, 'normInf') && ~isempty(lastInfo.normInf)
                finalNormInf = lastInfo.normInf(end);
            else
                finalNormInf = NaN;
            end
            warning(['Weighted PALC failed at continuation step %d: ', ...
                     'final ||F||inf = %.3e after backtracking.'], ...
                    k - 2, finalNormInf);
            break
        end

        branch.y{k} = y_corr;
        [X0k, Tk, rhok, muk] = gmos_unpack_y_mu(y_corr, N);
        branch.z{k} = gmos_pack_z(X0k, Tk, rhok);
        branch.mu(k) = muk;
        branch.T(k) = Tk;
        branch.rho(k) = rhok;
        branch.info{k} = lastInfo;
        branch.success(k) = true;
        branch.dsUsed(k) = dsTry;
        if opts.computeEGMOS
            [branch.EGMOS(k), ~] = gmos_accuracy_test_paper(branch.z{k}, N, muk, ode_opts);
        end

        if opts.verboseStep
            if opts.computeEGMOS
                fprintf(['PALC step %d/%d: mu = %.16g, T = %.15f, rho = %.15f, ', ...
                         'EGMOS = %.3e, ds = %.3e, final ||F||inf = %.3e\n'], ...
                        k - 2, nSteps, muk, Tk, rhok, branch.EGMOS(k), dsTry, lastInfo.normInf(end));
            else
                fprintf(['PALC step %d/%d: mu = %.16g, T = %.15f, rho = %.15f, ', ...
                         'ds = %.3e, final ||F||inf = %.3e\n'], ...
                        k - 2, nSteps, muk, Tk, rhok, dsTry, lastInfo.normInf(end));
            end
        end

        y_prevprev = y_prev;
        y_prev = y_corr;
        lastAccepted = k;
        dsCurrent = min(opts.growFactor * dsTry, opts.dsMax);
    end

    branch.y = branch.y(1:lastAccepted);
    branch.z = branch.z(1:lastAccepted);
    branch.mu = branch.mu(1:lastAccepted);
    branch.T = branch.T(1:lastAccepted);
    branch.rho = branch.rho(1:lastAccepted);
    branch.info = branch.info(1:lastAccepted);
    branch.success = branch.success(1:lastAccepted);
    branch.dsUsed = branch.dsUsed(1:lastAccepted);
    if opts.computeEGMOS
        branch.EGMOS = branch.EGMOS(1:lastAccepted);
    end
end
