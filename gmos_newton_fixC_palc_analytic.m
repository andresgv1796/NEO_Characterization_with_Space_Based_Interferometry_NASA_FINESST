function [z, info] = gmos_newton_fixC_palc_analytic( ...
    z0, N, mu, Xref, Tref, rhoref, Ctarget, z_prev, z_prevprev, ds, ode_opts, opts)
%GMOS_NEWTON_FIXC_PALC_ANALYTIC Newton corrector with fixed C and pseudo arclength.
%
% Solves:
%   Finv(X0,T,rho) = 0
%   p0(X0) = 0
%   p1(X0) = 0
%   s0(X0) = mean_j C(X0(:,j)) - Ctarget = 0
%   s1(z) = t'*(z - z_prev) - ds = 0
%
% Inputs
%   z0        : initial guess
%   N         : number of curve points
%   mu        : CR3BP mass parameter
%   Xref      : 6xN reference curve for phase conditions
%   Tref      : reference time for phase conditions
%   rhoref    : reference rho for phase conditions
%   Ctarget   : fixed average Jacobi constant
%   z_prev    : previous converged solution (for s1)
%   z_prevprev: solution before z_prev
%   ds        : continuation step
%   ode_opts  : (optional) odeset options
%   opts      : (optional) struct with fields maxIter, tolInf, verbose
%
% Outputs
%   z         : corrected solution
%   info      : iteration history and convergence info

    if nargin < 11
        ode_opts = [];
    end
    if nargin < 12 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'maxIter'); opts.maxIter = 10; end
    if ~isfield(opts, 'tolInf');  opts.tolInf = 1e-10; end
    if ~isfield(opts, 'verbose'); opts.verbose = true; end

    z = z0;

    hist_inf = zeros(opts.maxIter, 1);
    hist_2 = zeros(opts.maxIter, 1);
    converged = false;

    for k = 1:opts.maxIter
        [X0, T, rho] = gmos_unpack_z(z, N);

        [Finv_vec, Jinv, ~] = gmos_invariance_jacobian_analytic_fast(X0, T, rho, mu, ode_opts);

        [p0, p1, tang0, tang1] = gmos_phase_conditions(X0, Xref, mu, Tref, rhoref);

        s0 = gmos_constraint_fix_C(X0, mu, Ctarget);

        [s1, t_hat] = gmos_constraint_pseudo_arclength(z, z_prev, z_prevprev, ds);

        F = [Finv_vec; p0; p1; s0; s1];

        nX = 6*N;
        nZ = nX + 2;
        J = spalloc(nX + 4, nZ, nnz(Jinv) + 3*nX + nZ);

        J(1:nX, :) = Jinv;
        J(nX + 1, 1:nX) = (tang0(:).' / N);
        J(nX + 2, 1:nX) = (tang1(:).' / N);

        dC = gmos_jacobi_gradient(X0, mu);
        J(nX + 3, 1:nX) = (dC(:).' / N);

        J(nX + 4, :) = t_hat.';

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
            converged = true;
            break;
        end

        dz = lsqminnorm(J, F);
        z = z - dz;
    end

    if ~converged
        hist_inf = hist_inf(1:find(hist_inf ~= 0, 1, 'last'));
        hist_2 = hist_2(1:find(hist_2 ~= 0, 1, 'last'));
    end

    info = struct();
    info.normInf = hist_inf;
    info.norm2 = hist_2;
    info.converged = converged;
    info.iterations = numel(hist_inf);
end
