function [X0, T, rho, mu] = gmos_unpack_y_mu(y, N)
%GMOS_UNPACK_Y_MU Unpack extended GMOS unknown vector [X0(:); T; rho; mu].
%
% Inputs
%   y   : vector of length 6*N + 3
%   N   : number of curve nodes
%
% Outputs
%   X0  : 6xN invariant-curve samples
%   T   : scalar stroboscopic time
%   rho : scalar rotation number
%   mu  : scalar CR3BP mass parameter

    if numel(y) ~= 6*N + 3
        error('y must have length 6*N + 3.');
    end

    X0 = reshape(y(1:6*N), [6, N]);
    T = y(6*N + 1);
    rho = y(6*N + 2);
    mu = y(6*N + 3);
end
