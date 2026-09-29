function X1 = gmos_flowmap_curve_CR3BP(X0, T, mu, ode_opts)
% X1(:,j) = Phi^T(X0(:,j))


    if nargin < 4 || isempty(ode_opts)
        ode_opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
    end

    [n,N] = size(X0);
    if n ~= 6
        error('X0 must be 6xN.');
    end

    X1 = zeros(6,N);
    for j = 1:N
        sol = ode89(@(t,Y) CR3BP_mu(t,Y,mu), [0,T], X0(:,j), ode_opts);
        X1(:,j) = sol.y(:,end);
    end
end