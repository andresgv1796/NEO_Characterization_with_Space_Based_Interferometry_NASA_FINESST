function [raDeg, decDeg] = unit_to_radec(sHat)
%UNIT_TO_RADEC Convert a unit direction vector to RA/Dec in the same frame.
%
% [raDeg, decDeg] = unit_to_radec(sHat)
%
% The input vector must be expressed in the inertial frame whose +z axis is
% the celestial north pole and whose +x axis defines RA = 0. This is a
% coordinate conversion only; it does not transform between physical frames.

    sHat = normalize_row(sHat);

    raDeg = atan2d(sHat(2), sHat(1));
    if raDeg < 0
        raDeg = raDeg + 360;
    end

    decDeg = asind(max(-1, min(1, sHat(3))));
end
