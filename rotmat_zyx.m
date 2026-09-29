function R = rotmat_zyx(yawDeg, pitchDeg, rollDeg)
%ROTMAT_ZYX Construct a right-handed active Z-Y-X rotation matrix.
%
%   R = rotmat_zyx(yawDeg, pitchDeg, rollDeg)
%
%   Builds the rotation matrix
%
%       R = Rz(yaw) * Ry(pitch) * Rx(roll)
%
%   using standard right-handed ACTIVE rotations.
%
%   Inputs:
%       yawDeg     - Rotation about the +z axis [deg]
%       pitchDeg   - Rotation about the +y axis [deg]
%       rollDeg    - Rotation about the +x axis [deg]
%
%   Output:
%       R          - 3 x 3 active rotation matrix
%
%   Convention:
%
%       For a vector expressed as a column vector,
%
%           r_rot = R * r
%
%       physically rotates the vector/object while keeping the coordinate
%       frame fixed.
%
%       Positive rotations follow the right-hand rule.
%
%       Because
%
%           R = Rz * Ry * Rx,
%
%       the rotations are applied to a column vector in the order
%
%           1. roll  about +x
%           2. pitch about +y
%           3. yaw   about +z
%
%   Example:
%
%       R = rotmat_zyx(90, 0, 0);
%       rRot = R * [1; 0; 0];
%
%       gives approximately
%
%           rRot = [0; 1; 0]
%
%       corresponding to a +90 deg active rotation about +z.

arguments
    yawDeg   (1,1) double
    pitchDeg (1,1) double
    rollDeg  (1,1) double
end

% Convert input angles to radians.
yaw   = deg2rad(yawDeg);
pitch = deg2rad(pitchDeg);
roll  = deg2rad(rollDeg);

% -------------------------------------------------------------------------
% Active right-handed rotation about +z.
% -------------------------------------------------------------------------
Rz = [ ...
     cos(yaw), -sin(yaw), 0
     sin(yaw),  cos(yaw), 0
            0,         0, 1 ];

% -------------------------------------------------------------------------
% Active right-handed rotation about +y.
% -------------------------------------------------------------------------
Ry = [ ...
     cos(pitch),  0, sin(pitch)
              0,  1,          0
    -sin(pitch),  0, cos(pitch) ];

% -------------------------------------------------------------------------
% Active right-handed rotation about +x.
% -------------------------------------------------------------------------
Rx = [ ...
    1,          0,         0
    0,  cos(roll), -sin(roll)
    0,  sin(roll),  cos(roll) ];

% Combined active rotation.
R = Rz * Ry * Rx;

end