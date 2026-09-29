function proj = project_shape_orthographic(shape, R_IB, target)
%PROJECT_SHAPE_ORTHOGRAPHIC
% Rotate a body-frame asteroid mesh into the inertial frame and express
% the vertices in the target-centered interferometric sky frame.
%
% Frames
% ------
%   B : asteroid body frame
%   I : inertial CR3BP frame
%   S : target-centered sky / uvw frame
%
% For a fixed inertial target, the sky frame is
%
%       +x_S = +eU  -> +l image coordinate
%       +y_S = +eV  -> +m image coordinate
%       +z_S = +eW  -> observer-to-target / phase-center direction
%
% where
%
%       eU x eV = eW.
%
% Therefore, expressed in the sky frame, the direction FROM the target
% TOWARD the observer is
%
%       observerDirection_S = [0 0 -1].
%
% Inputs
% ------
% shape.vertices : Nv x 3 asteroid vertices in body frame [km]
% shape.faces    : Nf x 3 triangular connectivity
%
% R_IB : 3 x 3 active body-to-inertial attitude matrix
%
%       r_I = R_IB * r_B
%
% target : target structure returned by make_target_definition.
%
%          For the currently supported fixed-inertial case it must contain
%
%              target.eUInertial0
%              target.eVInertial0
%              target.eWInertial0
%
% Outputs
% -------
% proj.vertices3D : Nv x 3 vertices expressed in sky frame [km]
% proj.x          : coordinate along +eU [km]
% proj.y          : coordinate along +eV [km]
% proj.z          : coordinate along +eW [km]
% proj.faces      : triangle connectivity
%
% proj.observerDirection :
%       Direction from target toward observer, expressed in sky frame.

arguments
    shape struct
    R_IB (3,3) double
    target struct
end

% -------------------------------------------------------------------------
% Current implementation is for a fixed inertial target.
%
% A synodic target requires transforming its sky basis into the inertial
% frame at the current epoch before projecting the asteroid.
% -------------------------------------------------------------------------
if ~strcmpi(target.frame, 'inertial')
    error(['project_shape_orthographic currently requires an inertial ', ...
           'target definition.']);
end

% -------------------------------------------------------------------------
% Read the exact same sky basis used by the interferometric uvw model.
% -------------------------------------------------------------------------
eU = target.eUInertial0(:);
eV = target.eVInertial0(:);
eW = target.eWInertial0(:);

tol = 1e-10;

% -------------------------------------------------------------------------
% Validate body-to-inertial attitude.
% -------------------------------------------------------------------------
if norm(R_IB.'*R_IB - eye(3), 'fro') > tol
    error('R_IB must be an orthogonal rotation matrix.');
end

if abs(det(R_IB) - 1) > tol
    error('R_IB must be a proper rotation matrix with determinant +1.');
end

% -------------------------------------------------------------------------
% Validate target sky basis.
% -------------------------------------------------------------------------
if abs(norm(eU)-1) > tol || ...
   abs(norm(eV)-1) > tol || ...
   abs(norm(eW)-1) > tol

    error('Target sky-basis vectors must be unit vectors.');
end

if abs(dot(eU,eV)) > tol || ...
   abs(dot(eU,eW)) > tol || ...
   abs(dot(eV,eW)) > tol

    error('Target sky-basis vectors must be mutually orthogonal.');
end

if norm(cross(eU,eV) - eW) > tol
    error('Target sky basis must satisfy cross(eU,eV) = eW.');
end

% -------------------------------------------------------------------------
% Body frame -> inertial frame.
%
% Each row of shape.vertices contains a body-frame position vector.
% -------------------------------------------------------------------------
V_B = shape.vertices;

V_I = (R_IB * V_B.').';

% -------------------------------------------------------------------------
% Inertial frame -> target-centered sky frame.
%
% Because eU, eV and eW are expressed in inertial coordinates,
%
%       x_S = eU^T r_I
%       y_S = eV^T r_I
%       z_S = eW^T r_I
%
% Therefore
%
%       r_S = C_SI r_I
%
% with the basis vectors forming the rows of C_SI.
% -------------------------------------------------------------------------
C_SI = [ ...
    eU.'
    eV.'
    eW.' ];

V_S = (C_SI * V_I.').';

% -------------------------------------------------------------------------
% Store projected geometry.
%
% These coordinates are now explicitly conjugate to the same u,v,w axes
% used by the interferometer:
%
%       x -> l -> u
%       y -> m -> v
% -------------------------------------------------------------------------
proj.vertices3D = V_S;

proj.x = V_S(:,1);
proj.y = V_S(:,2);
proj.z = V_S(:,3);

proj.faces = shape.faces;

% +eW points observer -> target.
% Therefore target -> observer is -eW.
proj.observerDirection = [0, 0, -1];

% Metadata useful for diagnostics.
proj.R_IB = R_IB;
proj.C_SI = C_SI;

proj.eUInertial = eU.';
proj.eVInertial = eV.';
proj.eWInertial = eW.';

end