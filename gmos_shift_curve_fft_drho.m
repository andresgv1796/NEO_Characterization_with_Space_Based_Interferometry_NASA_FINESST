function Xdrho = gmos_shift_curve_fft_drho(X, rho)
%GMOS_SHIFT_CURVE_FFT_DRHO Derivative wrt rho of R_{-rho}(X).
%
% If R_{-rho} multiplies mode k by exp(-i*k*rho),
% then d/d rho multiplies mode k by (-i*k)*exp(-i*k*rho).

    [~, N] = size(X);

    Xhat = fft(X, [], 2);
    Xhat = fftshift(Xhat, 2);

    if mod(N,2) == 0
        k = (-N/2):(N/2-1);
    else
        k = (-(N-1)/2):((N-1)/2);
    end

    phase = exp(-1i * k * rho);
    mult = (-1i * k) .* phase;

    Xhat_drho = Xhat .* mult;

    Xhat_drho = ifftshift(Xhat_drho, 2);
    Xdrho = real(ifft(Xhat_drho, [], 2));
end