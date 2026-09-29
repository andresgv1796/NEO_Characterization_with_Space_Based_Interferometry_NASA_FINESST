function M_eigs = monodromy_CR3BP(tspan, IC, mu)
% monodromy_CR3BP  Compute eigenvectors and eigenvalues of the monodromy matrix
% for a periodic CR3BP orbit, returning a 6×12 matrix [V, D].
%
% Usage:
%   M_eigs = monodromy_CR3BP(tspan, IC, mu)
%
% Inputs:
%   tspan – 1×2 vector [t0, T] covering exactly one orbit period
%   IC    – 6×1 initial state [x; y; z; vx; vy; vz]
%   mu    – scalar mass parameter
%
% Output:
%   M_eigs – 6×12 matrix whose columns 1–6 are eigenvectors V and
%            whose columns 7–12 have the eigenvalues D on the diagonal

    % Validate inputs
    assert(numel(IC)==6,    'IC must be a 6×1 vector');
    assert(numel(tspan)==2, 'tspan must be a 1×2 vector');

    % Persistently store ODE options to avoid repeated allocation
    persistent ode_opts
    if isempty(ode_opts)
        ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    end

    % Initialize STM + state vector
    I     = eye(6);
    y0    = zeros(42,1);
    y0(1:36)   = I(:);    % flattened STM
    y0(37:42)  = IC(:);   % state vector

    % Integrate variational equations; only endpoint is needed
    sol = ode89(@(t,Y) STM_vec(t, Y, mu), tspan, y0, ode_opts);

    % Extract final STM and reshape to 6×6
    STM_T = reshape(sol.y(1:36,end), [6, 6]);

    % Compute eigen-decomposition
    [V, D] = eig(STM_T);

    % Combine into one output
    M_eigs = [V, D];
end
