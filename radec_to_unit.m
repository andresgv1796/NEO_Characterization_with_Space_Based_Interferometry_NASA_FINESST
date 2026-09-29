function s = radec_to_unit(raDeg, decDeg)
%RADEC_TO_UNIT Convert RA/Dec in degrees to an inertial unit vector.
%
% The inertial CR3BP frame is treated as the celestial reference frame:
% RA is measured about +z from +x, and Dec is measured from the x-y plane.

    [~, ~, s] = radec_to_uvw_basis(raDeg, decDeg);
end
