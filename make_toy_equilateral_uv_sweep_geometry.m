function geom = make_toy_equilateral_uv_sweep_geometry( ...
    sys, ...
    lambda_m, ...
    rhoMin_Glambda, ...
    rhoMax_Glambda, ...
    nAngles, ...
    nRadii)
%MAKE_TOY_EQUILATERAL_UV_SWEEP_GEOMETRY
%
% Construct a controlled three-collector toy trajectory designed to sweep
% the UV plane.
%
% The collectors form an equilateral triangle centered at the origin.
%
% The triangle:
%
%   1. changes scale to sweep baseline magnitude;
%   2. rotates to sweep baseline orientation.
%
% The geometry lies entirely in the x-y plane, so a target along +z is
% exactly normal to the collector plane at every epoch.
%
% The requested baseline magnitude is specified directly in Glambda.
%
% This is a geometric positive-control case, not a dynamical CR3BP
% trajectory.


%% Validate inputs

if rhoMin_Glambda <= 0
    error('rhoMin_Glambda must be positive.');
end

if rhoMax_Glambda <= rhoMin_Glambda
    error('rhoMax_Glambda must exceed rhoMin_Glambda.');
end

if nAngles < 2 || nAngles ~= floor(nAngles)
    error('nAngles must be an integer >= 2.');
end

if nRadii < 2 || nRadii ~= floor(nRadii)
    error('nRadii must be an integer >= 2.');
end


%% Nondimensional baseline-to-wavelength scale

lambda_km = ...
    lambda_m*1e-3;

scaleToWavelengths = ...
    sys.Lstar_km/lambda_km;


%% Requested B12 orientations
%
% For an equilateral triangle, the three unoriented baseline directions
% are separated by 60 deg.
%
% Therefore rotating one baseline through only [0,60) deg is sufficient
% for the complete three-baseline set to span [0,180) deg.

thetaB12DegList = ...
    (0:nAngles-1).' * ...
    (60/nAngles);


%% Requested baseline radii
%
% Use a simple linear radial sweep so the coverage is transparent and easy
% to audit.

rhoList_Glambda = ...
    linspace( ...
        rhoMin_Glambda, ...
        rhoMax_Glambda, ...
        nRadii);


%% Allocate trajectory

nEpoch = ...
    nAngles*nRadii;

r = zeros(nEpoch,3,3);

rhoRequested_Glambda = ...
    zeros(nEpoch,1);

thetaRequestedDeg = ...
    zeros(nEpoch,1);

angleIndex = ...
    zeros(nEpoch,1);

radiusIndex = ...
    zeros(nEpoch,1);


%% Construct expanding/contracting rotating triangle

k = 0;

for ia = 1:nAngles

    % Alternate radial direction so consecutive angular sweeps join at
    % either rhoMax or rhoMin instead of jumping across the whole range.

    if mod(ia,2) == 1

        rhoThis = ...
            rhoList_Glambda;

    else

        rhoThis = ...
            fliplr(rhoList_Glambda);

    end


    thetaB12Deg = ...
        thetaB12DegList(ia);


    % For collectors located at
    %
    %   phi
    %   phi + 120 deg
    %   phi + 240 deg
    %
    % the B12 direction is phi + 150 deg.
    %
    % Therefore choose phi so B12 has the requested direction.

    phi = ...
        deg2rad(thetaB12Deg - 150);


    for ir = 1:nRadii

        k = k + 1;


        %% Desired physical baseline length in nondimensional CR3BP units

        rhoWavelengths = ...
            rhoThis(ir)*1e9;

        baselineLengthND = ...
            rhoWavelengths / ...
            scaleToWavelengths;


        % Equilateral triangle:
        %
        %   side length = sqrt(3)*circumradius

        collectorRadiusND = ...
            baselineLengthND/sqrt(3);


        %% Collector angular locations

        collectorAngle = ...
            phi + ...
            [0, 2*pi/3, 4*pi/3];


        for c = 1:3

            r(k,1,c) = ...
                collectorRadiusND * ...
                cos(collectorAngle(c));

            r(k,2,c) = ...
                collectorRadiusND * ...
                sin(collectorAngle(c));

            r(k,3,c) = 0;

        end


        rhoRequested_Glambda(k) = ...
            rhoThis(ir);

        thetaRequestedDeg(k) = ...
            thetaB12Deg;

        angleIndex(k) = ia;
        radiusIndex(k) = ir;

    end

end


%% Toy trajectory parameter
%
% This is not CR3BP time. It is simply an ordered trajectory parameter.

t = ...
    linspace(0,1,nEpoch).';


%% Baselines

b12 = ...
    r(:,:,2) - r(:,:,1);

b13 = ...
    r(:,:,3) - r(:,:,1);

b23 = ...
    r(:,:,3) - r(:,:,2);


L12 = vecnorm(b12,2,2);
L13 = vecnorm(b13,2,2);
L23 = vecnorm(b23,2,2);


%% Collector-plane normal

normalCross = ...
    cross(b12,b13,2);

twoArea = ...
    vecnorm(normalCross,2,2);

degenerate = ...
    twoArea <= 100*eps;

if any(degenerate)
    error('Toy geometry unexpectedly contains degenerate triangles.');
end

normalSynodic = ...
    normalCross ./ twoArea;


% +z remains +z under the nominal synodic/inertial rotation, but construct
% this consistently with the rest of the pipeline.

normalInertial = ...
    rotate_synodic_to_inertial_rows( ...
        normalSynodic, ...
        t);


%% Output structure

geom = struct();

geom.t = t;
geom.r = r;

geom.baselines.b12 = b12;
geom.baselines.b13 = b13;
geom.baselines.b23 = b23;

geom.baselines.L12 = L12;
geom.baselines.L13 = L13;
geom.baselines.L23 = L23;

geom.normalSynodic = ...
    normalSynodic;

geom.normalInertial = ...
    normalInertial;

geom.degenerate = ...
    degenerate;


%% Toy metadata

geom.toy.description = ...
    'equilateral rotating/radially sweeping UV positive control';

geom.toy.rhoRequested_Glambda = ...
    rhoRequested_Glambda;

geom.toy.thetaB12RequestedDeg = ...
    thetaRequestedDeg;

geom.toy.angleIndex = ...
    angleIndex;

geom.toy.radiusIndex = ...
    radiusIndex;

geom.toy.nAngles = ...
    nAngles;

geom.toy.nRadii = ...
    nRadii;

geom.toy.rhoMin_Glambda = ...
    rhoMin_Glambda;

geom.toy.rhoMax_Glambda = ...
    rhoMax_Glambda;

end