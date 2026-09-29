function [Finv_vec, Jext, cache] = gmos_invariance_jacobian_analytic_muFD(X0, T, rho, mu, ode_opts, fdMu)
%GMOS_INVARIANCE_JACOBIAN_ANALYTIC_MUFD Analytic invariance Jacobian in z plus FD column in mu.
%
% Finv = R_{-rho}(Phi^T(X0;mu)) - X0
%
% Jacobian is taken with respect to the extended unknown
%   y = [X0(:); T; rho; mu].
%
% The [X0(:); T; rho] block is analytic via the STM. The mu column is
% approximated with a finite difference that preserves nonnegative mu.

    if nargin < 5
        ode_opts = [];
    end
    if nargin < 6 || isempty(fdMu)
        fdMu = min(1e-8, 0.1 * max(abs(mu), 1e-12));
    end
    if fdMu <= 0 || ~isfinite(fdMu)
        error('fdMu must be a positive finite scalar.');
    end

    [Finv_vec, Jinv, cache0] = gmos_invariance_jacobian_analytic_fast(X0, T, rho, mu, ode_opts);

    if mu - fdMu > 0
        % Centered difference when both sides stay physical.
        [Fplus, ~]  = gmos_invariance_jacobian_analytic_fast(X0, T, rho, mu + fdMu, ode_opts);
        [Fminus, ~] = gmos_invariance_jacobian_analytic_fast(X0, T, rho, mu - fdMu, ode_opts);
        J_mu = (Fplus - Fminus) / (2 * fdMu);
        fdScheme = 'centered';
    else
        % Forward difference fallback for very small mu.
        [Fplus, ~] = gmos_invariance_jacobian_analytic_fast(X0, T, rho, mu + fdMu, ode_opts);
        J_mu = (Fplus - Finv_vec) / fdMu;
        fdScheme = 'forward';
    end

    Jext = [Jinv, J_mu];

    cache = cache0;
    cache.fdMu = fdMu;
    cache.J_mu = J_mu;
    cache.fdScheme = fdScheme;
end
