function [y, info] = gmos_newton_fixT_mu_palc_analytic( ...
    y0, N, Xref, Tref, rhoref, muRef, Ttarget, y_prev, y_prevprev, ds, ode_opts, opts)
%GMOS_NEWTON_FIXT_MU_PALC_ANALYTIC Newton corrector for fixed T with mu as an unknown.
%
% Unknown:
%   y = [X0(:); T; rho; mu]
%
% Solves:
%   Finv(X0,T,rho;mu) = 0
%   p0(X0; Xref,Tref,rhoref,muRef) = 0
%   p1(X0; Xref,Tref,rhoref,muRef) = 0
%   s0(T) = T - Ttarget = 0
%   s1(y) = t_hat' * (y - y_prev) - ds = 0
%
% Notes
%   - The phase-condition tangents are computed from the reference torus
%     using its own parameter muRef, so they remain frozen during the
%     corrector.
%   - The mu Jacobian column for the invariance block is obtained with a
%     centered finite difference.

    if nargin < 11
        ode_opts = [];
    end
    if nargin < 12 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'maxIter'); opts.maxIter = 10; end
    if ~isfield(opts, 'tolInf');  opts.tolInf = 1e-10; end
    if ~isfield(opts, 'verbose'); opts.verbose = true; end
    if ~isfield(opts, 'fdMu');    opts.fdMu = []; end

    y = y0;

    hist_inf = zeros(opts.maxIter, 1);
    hist_2 = zeros(opts.maxIter, 1);

    for k = 1:opts.maxIter
        [X0, T, rho, mu] = gmos_unpack_y_mu(y, N);

        % Invariance residual and extended Jacobian.
        [Finv_vec, Jinv_ext, ~] = gmos_invariance_jacobian_analytic_muFD( ...
            X0, T, rho, mu, ode_opts, opts.fdMu);

        % Phase conditions frozen from the reference torus.
        [p0, p1, tang0, tang1] = gmos_phase_conditions(X0, Xref, muRef, Tref, rhoref);

        % Fixed T.
        s0 = T - Ttarget;

        % Pseudo arclength in the extended space.
        [s1, t_hat] = gmos_constraint_pseudo_arclength(y, y_prev, y_prevprev, ds);

        % Residual stack.
        F = [Finv_vec; p0; p1; s0; s1];

        nX = 6 * N;
        nY = nX + 3;

        J = spalloc(nX + 4, nY, nnz(Jinv_ext) + 2*nX + nY + 1);

        % Invariance rows.
        J(1:nX, :) = Jinv_ext;

        % p0 row.
        J(nX + 1, 1:nX) = (tang0(:).' / N);

        % p1 row.
        J(nX + 2, 1:nX) = (tang1(:).' / N);

        % s0 row: derivative wrt T only.
        J(nX + 3, nX + 1) = 1.0;

        % s1 row in the full extended space.
        J(nX + 4, :) = t_hat.';

        nInf = norm(F, inf);
        n2 = norm(F, 2);
        hist_inf(k) = nInf;
        hist_2(k) = n2;

        if opts.verbose
            fprintf(['Newton iter %2d: ||F||inf = %.3e, ||F||2 = %.3e, ', ...
                     'mu = %.16g\n'], k, nInf, n2, mu);
        end

        if nInf < opts.tolInf
            hist_inf = hist_inf(1:k);
            hist_2 = hist_2(1:k);
            break;
        end

        dy = lsqminnorm(J, F);
        y = y - dy;
    end

    info = struct();
    info.normInf = hist_inf;
    info.norm2 = hist_2;
end
