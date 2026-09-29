 
function dXdth = gmos_dtheta1_spectral(X)
%GMOS_DTHETA1_SPECTRAL Spectral derivative dX/dtheta1 along theta1 samples.
% Input
%   X       : 6xN samples X(:,j)=X(theta1_j)
%
% Output
%   dXdth   : 6xN samples of derivative with respect to theta1
    [~, N] = size(X);

    Xhat = fftshift(fft(X, [], 2), 2);

    if mod(N,2) == 0
        k = (-N/2):(N/2-1);
    else
        k = (-(N-1)/2):((N-1)/2);
    end

    Xhat_d = (1i * k) .* Xhat;

    dXdth_c = ifft(ifftshift(Xhat_d, 2), [], 2);

    % Keep real if the imaginary part is negligible
    if norm(imag(dXdth_c(:)), inf) <= 1e-12 * max(1, norm(real(dXdth_c(:)), inf))
        dXdth = real(dXdth_c);
    else
        dXdth = dXdth_c;
    end
end