function z = gmos_pack_z(X0, T, rho)
%GMOS_PACK_Z Pack GMOS unknowns into a single vector.
%
% Inputs
%   X0   : 6xN curve samples
%   T    : scalar stroboscopic time
%   rho  : scalar rotation number
%
% Output
%   z    : (6*N + 2)x1 vector [X0(:); T; rho]

    z = [X0(:); T; rho];
end