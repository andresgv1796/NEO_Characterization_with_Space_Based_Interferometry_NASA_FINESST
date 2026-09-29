function Xshift = gmos_shift_curve_fft(X, rho)
%GMOS_SHIFT_CURVE_FFT Spectral shift: returns X(theta1 - rho) for 2pi-periodic samples.
%
% Inputs
%   X   : 6xN (or mxN) samples along theta1
%   rho : scalar shift (radians)
%
% Output
%   Xshift : same size as X

    [~, N] = size(X);

    Xhat = fftshift(fft(X, [], 2), 2);

    if mod(N,2) == 0
        k = (-N/2):(N/2-1);
    else
        k = (-(N-1)/2):((N-1)/2);
    end

    Xhat_shift = Xhat .* exp(-1i * k * rho);

    Xshift_c = ifft(ifftshift(Xhat_shift, 2), [], 2);

    if norm(imag(Xshift_c(:)), inf) <= 1e-12 * max(1, norm(real(Xshift_c(:)), inf))
        Xshift = real(Xshift_c);
    else
        Xshift = Xshift_c;
    end
end