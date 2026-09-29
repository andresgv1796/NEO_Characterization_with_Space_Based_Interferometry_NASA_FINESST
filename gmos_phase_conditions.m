function [p0, p1, tang0, tang1] = gmos_phase_conditions(X0, Xref, mu, Tref, rhoref)
%GMOS_PHASE_CONDITIONS Paper-faithful GMOS phase conditions.
%
% Inputs
%   X0     : 6xN current curve samples v(theta1_j)
%   Xref   : 6xN reference curve samples vtilde(theta1_j)
%   mu     : CR3BP mass parameter
%   Tref   : scalar stroboscopic time associated with Xref (Ttilde)
%   rhoref : scalar rotation number associated with Xref (rhotilde)
%
% Outputs
%   p0     : (1/N) * sum_j (X0(:,j)-Xref(:,j))' * dXref/dtheta0(:,j) = 0
%   p1     : (1/N) * sum_j (X0(:,j))'           * dXref/dtheta1(:,j) = 0
%   tang0  : 6xN, dXref/dtheta0 evaluated from paper identity
%   tang1  : 6xN, dXref/dtheta1 (spectral derivative)

    [n, N] = size(X0);

    if n ~= 6 || any(size(Xref) ~= [6, N])
        error('X0 and Xref must both be 6xN.');
    end
    if ~(isscalar(Tref) && isscalar(rhoref))
        error('Tref and rhoref must be scalars.');
    end

    % Paper frequencies associated with the reference solution
    omega0 = 2*pi / Tref;
    omega1 = rhoref / Tref;

    % Tangent for theta1 shift: dXref/dtheta1
    tang1 = gmos_dtheta1_spectral(Xref);

    % Evaluate flow field along reference curve
    fref = zeros(6, N);
    for j = 1:N
        fref(:, j) = CR3BP_mu(0, Xref(:, j), mu);
    end

    % Paper identity: dXref/dtheta0 = (1/omega0) * ( fref - omega1 * dXref/dtheta1 )
    tang0 = (1/omega0) * (fref - omega1 * tang1);

    % Discrete inner products (paper uses 1/N scaling)
    dX = X0 - Xref;

    p0 = 0.0;
    p1 = 0.0;
    for j = 1:N
        p0 = p0 + dX(:, j).' * tang0(:, j);
        p1 = p1 + X0(:, j).' * tang1(:, j);
    end

    p0 = p0 / N;
    p1 = p1 / N;
end