function dC = gmos_jacobi_gradient(X0, mu)
%GMOS_JACOBI_GRADIENT Gradient of the CR3BP Jacobi constant at curve nodes.
%
% Inputs
%   X0 : 6xN curve samples [x; y; z; vx; vy; vz]
%   mu : CR3BP mass parameter
%
% Output
%   dC : 6xN, where dC(:,j) = grad C evaluated at X0(:,j)
%
% The Jacobi constant used here is consistent with Jacobi_Constant.m:
%
%   C = 2*U - (vx^2 + vy^2 + vz^2),
%   U = (1-mu)/r13 + mu/r23 + 0.5*(x^2 + y^2),
%
% with normalized mean motion n = 1.

    if size(X0, 1) ~= 6
        error('X0 must be 6xN.');
    end
    if ~isscalar(mu)
        error('mu must be a scalar.');
    end

    x  = X0(1, :);
    y  = X0(2, :);
    z  = X0(3, :);
    vx = X0(4, :);
    vy = X0(5, :);
    vz = X0(6, :);

    r13 = sqrt((x + mu).^2 + y.^2 + z.^2);
    r23 = sqrt((x - 1 + mu).^2 + y.^2 + z.^2);

    dCdx = 2*x - 2*(1 - mu).*(x + mu)./r13.^3 - 2*mu.*(x - 1 + mu)./r23.^3;
    dCdy = 2*y - 2*(1 - mu).*y./r13.^3      - 2*mu.*y./r23.^3;
    dCdz =      - 2*(1 - mu).*z./r13.^3      - 2*mu.*z./r23.^3;

    dCdvx = -2*vx;
    dCdvy = -2*vy;
    dCdvz = -2*vz;

    dC = [dCdx; dCdy; dCdz; dCdvx; dCdvy; dCdvz];
end
