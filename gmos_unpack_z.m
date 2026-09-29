function [X0, T, rho] = gmos_unpack_z(z, N)
%GMOS_UNPACK_Z Unpack GMOS unknown vector into X0, T, rho.
%
% Inputs
%   z   : (6*N + 2)x1 vector [X0(:); T; rho]
%   N   : number of curve points
%
% Outputs
%   X0  : 6xN curve samples
%   T   : scalar stroboscopic time
%   rho : scalar rotation number

    n = 6;

    if numel(z) ~= n*N + 2
        error('z must have length 6*N + 2.');
    end

    X0 = reshape(z(1:n*N), [n, N]);
    T = z(n*N + 1);
    rho = z(n*N + 2);
end