function [y, info] = gmos_newton_fixT_mu_palc_weighted( ...
    y0, N, Xref, Tref, rhoref, muRef, Ttarget, y_prev, y_prevprev, ds, w, ode_opts, opts)
%GMOS_NEWTON_FIXT_MU_PALC_WEIGHTED Weighted Newton corrector for fixed-T continuation in [z; mu].
%
% Unknown
%   y = [X0(:); T; rho; mu]
%
% Residual system
%   Finv(X0,T,rho;mu) = 0
%   p0(X0; Xref,Tref,rhoref,muRef) = 0
%   p1(X0; Xref,Tref,rhoref,muRef) = 0
%   s0(T) = T - Ttarget = 0
%   s1(y) = 0, weighted pseudo arclength constraint in ys = D*y
%
% Notes
%   1) The invariance Jacobian in [X0(:); T; rho] is analytic through the STM.
%   2) The extra mu column is a centered finite difference.
%   3) The phase-condition tangents are frozen on the supplied reference torus.
%   4) The PALC constraint is applied in the weighted scaled space to avoid the
%      6N state components overwhelming the single mu component.

    if nargin < 12
        ode_opts = [];
    end
    if nargin < 13 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'maxIter'); opts.maxIter = 10; end
    if ~isfield(opts, 'tolInf');  opts.tolInf = 1e-10; end
    if ~isfield(opts, 'verbose'); opts.verbose = true; end
    if ~isfield(opts, 'fdMu');    opts.fdMu = []; end

    y = y0(:);
    w = w(:);

    hist_inf = zeros(opts.maxIter, 1);
    hist_2 = zeros(opts.maxIter, 1);
    hist_mu = zeros(opts.maxIter, 1);

    for k = 1:opts.maxIter
        [X0, T, rho, mu] = gmos_unpack_y_mu(y, N);

        [Finv_vec, Jinv_ext, ~] = gmos_invariance_jacobian_analytic_muFD( ...
            X0, T, rho, mu, ode_opts, opts.fdMu);

        [p0, p1, tang0, tang1] = gmos_phase_conditions(X0, Xref, muRef, Tref, rhoref);
        s0 = T - Ttarget;
        [s1, Jpalc, ~] = gmos_constraint_pseudo_arclength_weighted( ...
            y, y_prev, y_prevprev, ds, w);

        F = [Finv_vec; p0; p1; s0; s1];

        nX = 6 * N;
        nY = nX + 3;
        J = spalloc(nX + 4, nY, nnz(Jinv_ext) + 2*nX + nY + 1);

        J(1:nX, :) = Jinv_ext;
        J(nX + 1, 1:nX) = tang0(:).' / N;
        J(nX + 2, 1:nX) = tang1(:).' / N;
        J(nX + 3, nX + 1) = 1.0;
        J(nX + 4, :) = Jpalc;

        nInf = norm(F, inf);
        n2 = norm(F, 2);
        hist_inf(k) = nInf;
        hist_2(k) = n2;
        hist_mu(k) = mu;

        if opts.verbose
            fprintf(['Newton iter %2d: ||F||inf = %.3e, ||F||2 = %.3e, ', ...
                     'mu = %.16g\n'], k, nInf, n2, mu);
        end

        if nInf < opts.tolInf
            hist_inf = hist_inf(1:k);
            hist_2 = hist_2(1:k);
            hist_mu = hist_mu(1:k);
            break;
        end

        dy = lsqminnorm(J, F);
        y = y - dy;
    end

    info = struct();
    info.normInf = hist_inf;
    info.norm2 = hist_2;
    info.muHist = hist_mu;
end
