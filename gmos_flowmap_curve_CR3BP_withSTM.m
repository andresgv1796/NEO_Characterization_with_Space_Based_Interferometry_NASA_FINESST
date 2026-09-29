function [X1, PhiBlocks, fEnd] = gmos_flowmap_curve_CR3BP_withSTM(X0, T, mu, ode_opts)
%GMOS_FLOWMAP_CURVE_CR3BP_WITHSTM Propagate curve points and STMs.
%
% Inputs
%   X0      : 6xN
%   T       : scalar
%   mu      : CR3BP mass parameter
%   ode_opts: odeset options (optional)
%
% Outputs
%   X1        : 6xN
%   PhiBlocks : 6x6xN, PhiBlocks(:,:,j) = dPhi^T/dx evaluated at X0(:,j)
%   fEnd      : 6xN, fEnd(:,j) = f(X1(:,j))

    if nargin < 4 || isempty(ode_opts)
        ode_opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    end

    [n, N] = size(X0);
    if n ~= 6
        error('X0 must be 6xN.');
    end

    X1 = zeros(6, N);
    fEnd = zeros(6, N);
    PhiBlocks = zeros(6, 6, N);

    Phi0 = eye(6);
    Phi0_vec = reshape(Phi0, 36, 1);

    for j = 1:N
        % Initial combined state: [vec(Phi); x]
        Y0 = [Phi0_vec; X0(:, j)];

        % STM_vec must propagate [vec(Phi); x] using A_CR3BP and CR3BP dynamics
        sol = ode89(@(t,Y) STM_vec(t, Y, mu), [0, T], Y0, ode_opts);
        YT = sol.y(:, end);

        PhiT_vec = YT(1:36);
        xT = YT(37:42);

        PhiBlocks(:, :, j) = reshape(PhiT_vec, 6, 6);
        X1(:, j) = xT;
        fEnd(:, j) = CR3BP_mu(0, xT, mu);
    end
end