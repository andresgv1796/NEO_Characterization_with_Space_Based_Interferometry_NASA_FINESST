function Sigma = clean_beam_covariance_mas2( ...
    fwhmMajor_mas,fwhmMinor_mas,paDeg)
%CLEAN_BEAM_COVARIANCE_MAS2
% Convert fitted Gaussian CLEAN-beam major/minor FWHM and position angle
% into a 2x2 covariance matrix in the (l,m) image basis.
%
% PA is measured counterclockwise from +l toward +m.

    fwhmToSigma = 1/(2*sqrt(2*log(2)));

    sigmaMajor = fwhmMajor_mas*fwhmToSigma;
    sigmaMinor = fwhmMinor_mas*fwhmToSigma;

    theta = deg2rad(paDeg);

    R = [ ...
        cos(theta),-sin(theta); ...
        sin(theta), cos(theta)];

    Sigma = R*diag([sigmaMajor^2,sigmaMinor^2])*R.';

end