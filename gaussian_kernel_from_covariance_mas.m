function K = gaussian_kernel_from_covariance_mas( ...
    Sigma_mas2,dl_mas,dm_mas,truncateSigma)
%GAUSSIAN_KERNEL_FROM_COVARIANCE_MAS
% Build a unit-sum 2-D Gaussian smoothing kernel for an arbitrary positive
% definite covariance matrix expressed in mas^2 in the (l,m) basis.

    Sigma_mas2 = 0.5*(Sigma_mas2+Sigma_mas2.');

    [Q,D] = eig(Sigma_mas2);
    d = real(diag(D));

    tol = 1e-14*max(1,max(abs(d)));

    if min(d) < -tol
        error('Gaussian-kernel covariance is not positive semidefinite.');
    end

    d = max(d,tol);

    SigmaUse = Q*diag(d)*Q.';
    invSigma = inv(SigmaUse);

    sigmaMax = sqrt(max(d));

    nL = max(1,ceil(truncateSigma*sigmaMax/abs(dl_mas)));
    nM = max(1,ceil(truncateSigma*sigmaMax/abs(dm_mas)));

    l = (-nL:nL)*abs(dl_mas);
    m = (-nM:nM)*abs(dm_mas);

    [L,M] = meshgrid(l,m);

    q = ...
        invSigma(1,1)*L.^2 + ...
        2*invSigma(1,2)*L.*M + ...
        invSigma(2,2)*M.^2;

    K = exp(-0.5*q);
    K = K/sum(K(:));

end