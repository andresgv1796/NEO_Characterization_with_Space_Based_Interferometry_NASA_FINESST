function qpos = gmos_generate_qpos_iterates(z, N, mu, nIter, ode_opts)
%GMOS_GENERATE_QPOS_ITERATES Generate qpos points by iterating the stroboscopic map.
%
% Inputs
%   z       : one solution vector [X0(:); T; rho]
%   N       : curve discretization
%   mu      : CR3BP mass parameter
%   nIter   : number of stroboscopic iterates
%   ode_opts: odeset options (optional)
%
% Output
%   qpos    : 3 x (N*nIter) array of positions from all iterates

    if nargin < 5
        ode_opts = [];
    end

    [X0, T, rho] = gmos_unpack_z(z, N);

    X = X0;
    qpos = zeros(3, N*nIter);

    for k = 1:nIter
        % Store current positions
        idx = (k-1)*N + (1:N);
        qpos(:, idx) = X(1:3, :);

        % Stroboscopic map and reparameterization to keep curve aligned
        X1 = gmos_flowmap_curve_CR3BP(X, T, mu, ode_opts);
        X = gmos_shift_curve_fft(X1, rho);
    end
end