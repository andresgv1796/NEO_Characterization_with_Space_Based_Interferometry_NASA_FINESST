function targets = fibonacci_sphere_targets(nTargets)
%FIBONACCI_SPHERE_TARGETS Deterministic quasi-uniform unit-sphere sampling.

    narginchk(1,1);
    validateattributes(nTargets,{'numeric'}, ...
        {'real','finite','scalar','integer','>=',4},mfilename,'nTargets',1);

    k = (0:nTargets-1).';
    z = 1-2*(k+0.5)/nTargets;
    goldenAngle = pi*(3-sqrt(5));
    phi = mod(k*goldenAngle,2*pi);
    rxy = sqrt(1-z.^2);
    x = rxy.*cos(phi);
    y = rxy.*sin(phi);
    S = [x y z];

    normError = max(abs(vecnorm(S,2,2)-1));
    if normError > 100*eps
        error('FINESST:FibonacciSphere:UnitVectorFailure', ...
            'Fibonacci target normalization error reached %.3e.',normError);
    end

    targets = struct();
    targets.sInertial = S;
    targets.longitudeDeg = mod(atan2d(y,x),360);
    targets.latitudeDeg = asind(z);
    targets.frame = [ ...
        'model-inertial axes assumed aligned with J2000 at phase zero; ' ...
        'no epoch-dependent ICRF orientation'];
    targets.method = 'fibonacci-midpoint';
    targets.nTargets = nTargets;
end
