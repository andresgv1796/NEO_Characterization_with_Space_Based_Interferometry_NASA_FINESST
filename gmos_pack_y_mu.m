function y = gmos_pack_y_mu(X0, T, rho, mu)
%GMOS_PACK_Y_MU Pack GMOS unknowns plus continuation parameter mu.
%
% Inputs
%   X0   : 6xN invariant-curve samples
%   T    : scalar stroboscopic time
%   rho  : scalar rotation number
%   mu   : scalar CR3BP mass parameter
%
% Output
%   y    : [X0(:); T; rho; mu]

    y = [X0(:); T; rho; mu];
end
