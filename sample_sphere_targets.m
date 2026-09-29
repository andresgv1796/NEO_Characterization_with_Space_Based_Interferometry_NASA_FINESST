function targets = sample_sphere_targets(nTargets, method, rngSeed)
%SAMPLE_SPHERE_TARGETS Sample target directions on the unit sphere.
%
% targets = sample_sphere_targets(nTargets)
% targets = sample_sphere_targets(nTargets, method)
% targets = sample_sphere_targets(nTargets, method, rngSeed)
%
% Methods:
%   'fibonacci'  : deterministic quasi-uniform equal-area sampling.
%                  Best default for suitability scans.
%   'random'     : random uniform sampling on the sphere.
%   'equalarea'  : deterministic RA/z grid with z = sin(dec) uniformly spaced.
%   'latlon'     : regular RA/Dec grid. Useful for debugging, not area-uniform.
%
% Output fields:
%   sInertial : nTargets x 3 unit vectors
%   raDeg     : right-ascension-like angle in the chosen inertial frame [deg]
%   decDeg    : declination-like angle in the chosen inertial frame [deg]
%
% Notes:
%   RA/Dec here are coordinates in the same inertial basis used by the CR3BP
%   UVW pipeline. For real catalog targets, first transform the catalog
%   direction into this inertial frame.

    if nargin < 1 || isempty(nTargets)
        nTargets = 300;
    end

    if nargin < 2 || isempty(method)
        method = 'fibonacci';
    end

    if nargin < 3
        rngSeed = [];
    end

    nTargets = round(nTargets);

    if nTargets < 1
        error('nTargets must be positive.');
    end

    methodKey = lower(strrep(strtrim(method), '-', ''));
    methodKey = strrep(methodKey, '_', '');

    switch methodKey

        case {'fibonacci', 'fib', 'golden', 'goldenangle'}
            % Quasi-uniform equal-area sphere sampling.
            %
            % The key point is that z = sin(dec) is uniformly spaced in
            % [-1, 1], while longitude advances by the golden angle. This
            % avoids pole clustering and avoids a regular longitude seam.
            k = (0:nTargets-1).';
            goldenAngle = pi*(3.0 - sqrt(5.0));

            z = 1.0 - 2.0*(k + 0.5)/nTargets;
            r = sqrt(max(0.0, 1.0 - z.^2));
            theta = goldenAngle*k;

            x = r.*cos(theta);
            y = r.*sin(theta);

            methodName = 'fibonacci';
            note = ['Quasi-uniform equal-area sampling: z = sin(dec) is ', ...
                    'uniformly spaced and RA advances by the golden angle.'];

        case {'random', 'uniformrandom', 'montecarlo'}
            % Random uniform sphere sampling.
            %
            % Uniform dec is wrong. Uniform z = sin(dec) is correct.
            if ~isempty(rngSeed)
                rng(rngSeed);
            end

            z = 2.0*rand(nTargets, 1) - 1.0;
            theta = 2.0*pi*rand(nTargets, 1);
            r = sqrt(max(0.0, 1.0 - z.^2));

            x = r.*cos(theta);
            y = r.*sin(theta);

            methodName = 'random';
            note = 'Random uniform sphere sampling with z = sin(dec) uniform in [-1, 1].';

        case {'equalarea', 'equalareagrid', 'radecequalarea', 'zra'}
            % Deterministic equal-area RA/z grid.
            %
            % This is more grid-like than Fibonacci but still area-aware
            % because z = sin(dec), not dec, is uniformly spaced.
            nDec = max(2, round(sqrt(nTargets/2)));
            nRa = ceil(nTargets/nDec);

            zGrid = linspace(1.0 - 1.0/nDec, -1.0 + 1.0/nDec, nDec).';
            raGrid = linspace(0.0, 2.0*pi, nRa + 1).';
            raGrid(end) = [];

            [thetaMat, zMat] = meshgrid(raGrid, zGrid);
            theta = thetaMat(:);
            z = zMat(:);

            theta = theta(1:nTargets);
            z = z(1:nTargets);

            r = sqrt(max(0.0, 1.0 - z.^2));
            x = r.*cos(theta);
            y = r.*sin(theta);

            methodName = 'equalarea';
            note = 'Deterministic RA/z grid. z = sin(dec) is uniformly spaced.';

        case {'latlon', 'radec', 'radecgrid', 'grid'}
            % Regular RA/Dec grid.
            %
            % This is NOT area-uniform because equal dec spacing clusters
            % samples toward the poles. It is useful for visual debugging.
            nDec = max(2, round(sqrt(nTargets/2)));
            nRa = ceil(nTargets/nDec);

            decGrid = linspace(-90.0, 90.0, nDec).';
            raGrid = linspace(0.0, 360.0, nRa + 1).';
            raGrid(end) = [];

            [raMat, decMat] = meshgrid(raGrid, decGrid);
            raDegRaw = raMat(:);
            decDegRaw = decMat(:);

            raDegRaw = raDegRaw(1:nTargets);
            decDegRaw = decDegRaw(1:nTargets);

            a = deg2rad(raDegRaw);
            d = deg2rad(decDegRaw);

            x = cos(d).*cos(a);
            y = cos(d).*sin(a);
            z = sin(d);

            methodName = 'latlon';
            note = 'Regular RA/Dec grid. Useful for debugging, but not area-uniform.';

        otherwise
            error(['Unknown sphere sampling method "%s". Use ', ...
                   '''fibonacci'', ''random'', ''equalarea'', or ''latlon''.'], method);
    end

    s = [x, y, z];
    s = s ./ vecnorm(s, 2, 2);

    raDeg = mod(atan2d(s(:,2), s(:,1)), 360.0);
    decDeg = asind(max(-1.0, min(1.0, s(:,3))));

    targets = struct();
    targets.sInertial = s;
    targets.raDeg = raDeg;
    targets.decDeg = decDeg;
    targets.nTargets = size(s, 1);
    targets.requestedNTargets = nTargets;
    targets.method = methodName;
    targets.note = note;
end
