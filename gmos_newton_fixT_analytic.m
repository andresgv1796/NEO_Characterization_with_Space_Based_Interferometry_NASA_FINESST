function [z, info] = gmos_newton_fixT_analytic(z0, N, mu, Xref, Tref, rhoref, Ttarget, ode_opts, opts)
%GMOS_NEWTON_FIXT_ANALYTIC Paper-style Newton corrector with analytic Jacobian (fixed T).
%
% Solves F(z)=0 in least-squares sense using:
%   dz = argmin ||J*dz - F||_2  (lsqminnorm)
%   z <- z - dz
%
% Unknowns are z = [X0(:); T; rho], but T is constrained by s0(T)=T-Ttarget.
%
% Inputs
%   z0      : initial guess [X0(:); T; rho]
%   N       : number of curve points
%   mu      : CR3BP mass parameter
%   Xref    : 6xN reference curve for phase conditions
%   Tref    : scalar reference stroboscopic time
%   rhoref  : scalar reference rotation number
%   Ttarget : scalar target stroboscopic time (fixed)
%   ode_opts: (optional) odeset options
%   opts    : (optional) struct
%               maxIter (default 10)
%               tolInf  (default 1e-10)
%               verbose (default true)
%
% Outputs
%   z       : solution (or last iterate)
%   info    : iteration history

    if nargin < 8
        ode_opts = [];
    end
    if nargin < 9 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'maxIter'); opts.maxIter = 10; end
    if ~isfield(opts, 'tolInf');  opts.tolInf = 1e-10; end
    if ~isfield(opts, 'verbose'); opts.verbose = true; end

    z = z0;

    hist_inf = zeros(opts.maxIter, 1);
    hist_2 = zeros(opts.maxIter, 1);

    for k = 1:opts.maxIter
        % Unpack
        [X0, T, rho] = gmos_unpack_z(z, N);

        % Invariance residual and Jacobian
        [Finv_vec, Jinv, ~] = gmos_invariance_jacobian_analytic_fast(X0, T, rho, mu, ode_opts);

        % Phase conditions (paper)
        [p0, p1, tang0, tang1] = gmos_phase_conditions(X0, Xref, mu, Tref, rhoref);

        % Fixed T constraint
        s0 = T - Ttarget;

        % Full residual
        F = [Finv_vec; p0; p1; s0];

        % Full Jacobian
        nX = 6*N;
        nZ = nX + 2;
        J = spalloc(nX + 3, nZ, nnz(Jinv) + 2*nX + 1);

        % Invariance rows
        J(1:nX, :) = Jinv;

        % p0 row: derivative wrt X0(:) only
        J(nX + 1, 1:nX) = (tang0(:).' / N);

        % p1 row: derivative wrt X0(:) only
        J(nX + 2, 1:nX) = (tang1(:).' / N);

        % s0 row: derivative wrt T only
        J(nX + 3, nX + 1) = 1.0;

        % Norms
        nInf = norm(F, inf);
        n2 = norm(F, 2);
        hist_inf(k) = nInf;
        hist_2(k) = n2;

        if opts.verbose
            fprintf('Newton iter %2d: ||F||inf = %.3e, ||F||2 = %.3e\n', k, nInf, n2);
        end

        if nInf < opts.tolInf
            hist_inf = hist_inf(1:k);
            hist_2 = hist_2(1:k);
            break;
        end

        % Paper-style pseudoinverse step (least squares)
        dz = lsqminnorm(J, F);

        % Update
        z = z - dz;
    end

    info = struct();
    info.normInf = hist_inf;
    info.norm2 = hist_2;
end