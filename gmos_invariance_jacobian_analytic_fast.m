function [Finv_vec, Jinv, cache] = gmos_invariance_jacobian_analytic_fast(X0, T, rho, mu, ode_opts)
%GMOS_INVARIANCE_JACOBIAN_ANALYTIC_FAST Analytic Finv(:) and Jacobian for invariance block.
%
% Finv = R_{-rho}(Phi^T(X0)) - X0
% Jinv is Jacobian of Finv(:) with respect to z = [X0(:); T; rho]
%
% Inputs
%   X0      : 6xN curve samples
%   T       : scalar propagation time
%   rho     : scalar rotation number
%   mu      : CR3BP mass parameter
%   ode_opts: (optional) odeset options
%
% Outputs
%   Finv_vec : (6N)x1 vector Finv(:)
%   Jinv     : (6N)x(6N+2) Jacobian for [X0(:); T; rho]
%   cache    : struct with intermediate values

    if nargin < 5
        ode_opts = [];
    end

    [n, N] = size(X0);
    if n ~= 6
        error('X0 must be 6xN.');
    end

    % Propagate each point and its STM
    [X1, PhiBlocks, fEnd] = gmos_flowmap_curve_CR3BP_withSTM(X0, T, mu, ode_opts);

    % Shift operator as an explicit NxN matrix: S such that Xshift = X * S
    % We can build it by shifting the identity matrix since gmos_shift_curve_fft
    % acts row-wise along the column index.
    S = gmos_shift_curve_fft(eye(N), rho);     % NxN
    ST_kron = kron(S.', speye(6));             % (6N)x(6N), sparse

    % Invariance residual
    X1shift = X1 * S;                          % 6xN
    Finv = X1shift - X0;                       % 6xN
    Finv_vec = Finv(:);

    % Build block diagonal STM matrix: PhiBig maps vec(dX0) -> vec(dX1) (no T term)
    PhiBig = spalloc(6*N, 6*N, 36*N);
    for j = 1:N
        rows = (6*(j-1)+1):(6*j);
        cols = rows;
        PhiBig(rows, cols) = PhiBlocks(:, :, j);
    end

    % Derivative wrt X0(:): ST_kron * PhiBig - I
    J_X0 = ST_kron * PhiBig - speye(6*N);

    % Derivative wrt T: shift of fEnd
    % vec(R_{-rho}(fEnd)) = vec(fEnd * S) = (S.' ⊗ I) vec(fEnd)
    RT_fEnd = fEnd * S;                        % 6xN
    J_T = RT_fEnd(:);                          % (6N)x1

    % Derivative wrt rho: d/d rho (X1 * S(rho)) = X1 * dS/drho
    % We can compute dS/drho by applying the rho-derivative shift to I.
    dSdrho = gmos_shift_curve_fft_drho(eye(N), rho);  % NxN
    dRho_X1 = X1 * dSdrho;                     % 6xN
    J_rho = dRho_X1(:);                        % (6N)x1

    % Assemble full invariance Jacobian
    Jinv = [J_X0, J_T, J_rho];

    cache = struct();
    cache.X1 = X1;
    cache.X1shift = X1shift;
    cache.PhiBlocks = PhiBlocks;
    cache.fEnd = fEnd;
    cache.S = S;
    cache.dSdrho = dSdrho;
end