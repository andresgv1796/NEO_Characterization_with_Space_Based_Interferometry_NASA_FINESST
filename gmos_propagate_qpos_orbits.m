function traj = gmos_propagate_qpos_orbits(z, N, mu, jList, nPeriods, nSavePerPeriod, ode_opts)
%GMOS_PROPAGATE_QPOS_ORBITS Propagate selected invariant-curve samples to 3D trajectories.
%
% Inputs
%   z              : GMOS solution vector [X0(:); T; rho]
%   N              : number of curve points in X0
%   mu             : CR3BP mass parameter
%   jList          : indices of curve points to propagate (e.g., 1:10:N)
%   nPeriods       : number of stroboscopic periods to propagate (total time = nPeriods*T)
%   nSavePerPeriod : number of saved points per period (for uniform output sampling)
%   ode_opts       : odeset options (optional)
%
% Output (struct)
%   traj.t         : 1xK time vector (uniform grid)
%   traj.Y         : 6xKxM array of states, M = numel(jList)
%   traj.jList     : j indices used
%   traj.T         : stroboscopic time used
%   traj.rho       : rotation number stored

    if nargin < 7 || isempty(ode_opts)
        ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    end

    [X0, T, rho] = gmos_unpack_z(z, N);

    jList = jList(:).';
    M = numel(jList);

    tEnd = nPeriods * T;
    K = nPeriods * nSavePerPeriod + 1;
    tGrid = linspace(0, tEnd, K);

    Y = zeros(6, K, M);

    for m = 1:M
        j = jList(m);
        x0 = X0(:, j);

        % Or use ode113 if robust long integrations are needed
        sol = ode89(@(t,y) CR3BP_mu(t, y, mu), [0, tEnd], x0, ode_opts);

        % Sample onto uniform grid for easy plotting/comparison
        Y(:, :, m) = deval(sol, tGrid);
    end

    traj = struct();
    traj.t = tGrid;
    traj.Y = Y;
    traj.jList = jList;
    traj.T = T;
    traj.rho = rho;
end