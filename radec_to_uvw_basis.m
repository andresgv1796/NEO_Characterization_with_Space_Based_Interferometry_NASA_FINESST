function [eU, eV, eW] = radec_to_uvw_basis(raDeg, decDeg)
%RADEC_TO_UVW_BASIS Standard source-centered interferometry basis.
%
% [eU, eV, eW] = radec_to_uvw_basis(raDeg, decDeg)
%
% Inputs:
%   raDeg  : right ascension angle in degrees, measured about +z from +x
%   decDeg : declination angle in degrees, measured from the x-y plane
%
% Outputs:
%   eU : east direction on the sky
%   eV : north direction on the sky
%   eW : source / phase-center direction
%
% Convention:
%   eW points to the target phase center, eU points east, eV points north,
%   and eU x eV = eW. Baseline components along this basis divided by
%   wavelength are the conventional (u, v, w) coordinates.
%
% Note:
%   The inertial CR3BP frame is treated as the celestial reference frame for
%   this calculation. For an actual catalog target, first express the ICRF
%   source direction in the CR3BP inertial basis.

    ra = deg2rad(raDeg);
    dec = deg2rad(decDeg);

    eW = [cos(dec)*cos(ra), cos(dec)*sin(ra), sin(dec)];

    eU = [-sin(ra), cos(ra), 0];

    eV = [-sin(dec)*cos(ra), -sin(dec)*sin(ra), cos(dec)];

    eU = normalize_row(eU);
    eV = normalize_row(eV);
    eW = normalize_row(eW);
end
