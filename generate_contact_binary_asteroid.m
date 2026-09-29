function shape = generate_contact_binary_asteroid( ...
    axes1, axes2, offset2, N, R1, R2)
%GENERATE_CONTACT_BINARY_ASTEROID
% Build a composite bilobed/contact-binary asteroid.
%
%   shape = generate_contact_binary_asteroid( ...
%       axes1, axes2, offset2, N, R1, R2)
%
% The model consists of two independently oriented triaxial ellipsoids.
% Lobe 1 is initially centered at the origin. Lobe 2 is translated by
% offset2. The complete composite body is then recentered at an
% approximate center of mass based on the analytical ellipsoid volumes.
%
% Inputs
% ------
% axes1 :
%       [a1 b1 c1] semi-axes of lobe 1 [km].
%
% axes2 :
%       [a2 b2 c2] semi-axes of lobe 2 [km].
%
% offset2 :
%       [dx dy dz] center of lobe 2 relative to lobe 1 [km].
%
% N :
%       Mesh resolution passed to generate_triaxial_ellipsoid.
%
% R1, R2 :
%       3 x 3 active right-handed rotations defining the orientation of
%       each lobe within the asteroid body geometry.
%
% Output
% ------
% shape.vertices :
%       Composite mesh vertices expressed relative to the approximate
%       center of mass [km].
%
% shape.faces :
%       Composite triangular connectivity.
%
% shape.type :
%       Shape-family identifier.
%
% shape.params :
%       Geometry and reference-center metadata.
%
% Important limitation
% --------------------
% The two ellipsoids are geometrically superposed but are not Boolean-
% merged into one watertight surface. Therefore the COM calculated here
% treats both complete ellipsoid volumes independently and double-counts
% any overlap volume.
%
% This is an APPROXIMATE COM suitable for the current synthetic imaging
% model. A future watertight contact-binary model should use the centroid
% of the actual union volume.

arguments
    axes1 (1,3) double {mustBePositive}
    axes2 (1,3) double {mustBePositive}
    offset2 (1,3) double
    N (1,1) double {mustBeInteger,mustBePositive} = 60
    R1 (3,3) double = eye(3)
    R2 (3,3) double = eye(3)
end


%% Validate lobe rotations

tol = 1e-10;

if norm(R1.'*R1 - eye(3),'fro') > tol || ...
   abs(det(R1) - 1) > tol
    error('R1 must be a proper rotation matrix.');
end

if norm(R2.'*R2 - eye(3),'fro') > tol || ...
   abs(det(R2) - 1) > tol
    error('R2 must be a proper rotation matrix.');
end


%% Generate the two ellipsoid meshes

shape1 = generate_triaxial_ellipsoid( ...
    axes1(1), axes1(2), axes1(3), N);

shape2 = generate_triaxial_ellipsoid( ...
    axes2(1), axes2(2), axes2(3), N);

V1 = shape1.vertices;
F1 = shape1.faces;

V2 = shape2.vertices;
F2 = shape2.faces;


%% Orient each lobe about its own center

V1 = (R1 * V1.').';
V2 = (R2 * V2.').';


%% Define lobe centers before composite recentering

center1 = [0, 0, 0];
center2 = offset2;


%% Translate lobe 2

V2 = V2 + center2;


%% Approximate lobe masses from analytical ellipsoid volumes

volume1 = (4/3)*pi*prod(axes1);
volume2 = (4/3)*pi*prod(axes2);

% Equal-density assumption:
%
%       mass_i = density * volume_i
%
% The common density cancels in the COM calculation.
massWeight1 = volume1;
massWeight2 = volume2;


%% Approximate composite center of mass

centerCOM = ...
    (massWeight1*center1 + massWeight2*center2) / ...
    (massWeight1 + massWeight2);


%% Recenter both lobes about the approximate COM

V1 = V1 - centerCOM;
V2 = V2 - centerCOM;

center1COM = center1 - centerCOM;
center2COM = center2 - centerCOM;


%% Concatenate meshes

F2 = F2 + size(V1,1);

V = [V1; V2];
F = [F1; F2];


%% Package output

shape.vertices = V;
shape.faces = F;

shape.type = 'contact binary';


%% Store geometric parameters

shape.params.axes1 = axes1;
shape.params.axes2 = axes2;

shape.params.R1 = R1;
shape.params.R2 = R2;

shape.params.offset2 = offset2;
shape.params.N = N;


%% Store reference-center information

shape.params.volume1_km3 = volume1;
shape.params.volume2_km3 = volume2;

shape.params.massWeight1 = massWeight1;
shape.params.massWeight2 = massWeight2;

% COM location in the original lobe-1-centered coordinates.
shape.params.centerCOMOriginal_km = centerCOM;

% Lobe centers after moving the asteroid body-frame origin to the COM.
shape.params.center1Body_km = center1COM;
shape.params.center2Body_km = center2COM;

shape.params.centerDefinition = ...
    'equal-density volume-weighted lobe COM; overlap volume double-counted';

end