function target = make_target_definition( ...
    targetFrame, useRaDec, raDeg, decDeg, sTargetInertial, sTargetSynodic, ...
    skyNorthReference)
%MAKE_TARGET_DEFINITION Build target direction and standard uvw sky basis.
%
% Standard fixed-target interferometry convention:
%   eW points to the target / phase center.
%   eU points east on the sky, i.e., increasing RA / longitude.
%   eV points north on the sky, i.e., increasing Dec / latitude.
%   eU x eV = eW.
%
% If useRaDec = true, the basis is built directly from RA/Dec.
%
% If useRaDec = false and targetFrame = 'inertial', the input vector is
% converted to RA/Dec in the same inertial coordinate system, then the same
% RA/Dec basis is used. Therefore vector-defined and RA/Dec-defined fixed
% inertial targets are equivalent, provided the vector is expressed in the
% same inertial frame where +z is the sky north pole and +x is RA = 0.
%
% For a real catalog target, first transform the ICRF/J2000 target direction
% into the CR3BP inertial basis, then pass either the resulting vector or its
% corresponding RA/Dec in that same basis.

    targetFrame = lower(strtrim(targetFrame));
    if ~strcmp(targetFrame, 'inertial') && ~strcmp(targetFrame, 'synodic')
        error('targetFrame must be ''inertial'' or ''synodic''.');
    end

    target = struct();
    target.frame = char(targetFrame);
    target.useRaDec = useRaDec;
    target.skyNorthReference = normalize_row(skyNorthReference);

    if useRaDec
        if ~strcmp(targetFrame, 'inertial')
            error('RA/Dec defines a fixed inertial target. Use targetFrame = ''inertial''.');
        end

        [eU, eV, eW] = radec_to_uvw_basis(raDeg, decDeg);

        target.raDeg = raDeg;
        target.decDeg = decDeg;
        target.sInertial0 = eW;
        target.eUInertial0 = eU;
        target.eVInertial0 = eV;
        target.eWInertial0 = eW;
        target.basisConvention = 'standard fixed-target uvw: u=east, v=north, w=phase center';
        target.description = sprintf( ...
            'fixed inertial phase center from RA %.6f deg, Dec %.6f deg', ...
            raDeg, decDeg);

    elseif strcmp(targetFrame, 'inertial')
        target.sInertial0 = normalize_row(sTargetInertial);

        % Convert the vector to RA/Dec in the same inertial frame, then use
        % the same standard east/north/source basis as the RA/Dec branch.
        [raFromVecDeg, decFromVecDeg] = unit_to_radec(target.sInertial0);
        [eU, eV, eW] = radec_to_uvw_basis(raFromVecDeg, decFromVecDeg);

        target.raDeg = raFromVecDeg;
        target.decDeg = decFromVecDeg;
        target.eUInertial0 = eU;
        target.eVInertial0 = eV;
        target.eWInertial0 = eW;
        target.basisConvention = 'standard fixed-target uvw from vector converted to RA/Dec';
        target.description = sprintf( ...
            'fixed inertial phase center from vector, equivalent RA %.6f deg, Dec %.6f deg', ...
            raFromVecDeg, decFromVecDeg);

    else
        % A synodic target is useful for CR3BP-relative diagnostics, but it
        % is not a fixed celestial source. Use this only for dynamics tests.
        target.sSynodic0 = normalize_row(sTargetSynodic);
        [eU, eV, eW] = build_sky_basis_from_target( ...
            target.sSynodic0, target.skyNorthReference);

        target.raDeg = NaN;
        target.decDeg = NaN;
        target.eUSynodic0 = eU;
        target.eVSynodic0 = eV;
        target.eWSynodic0 = eW;
        target.basisConvention = 'synodic diagnostic uvw: v from projected reference, u=v x w';
        target.description = 'fixed synodic CR3BP-relative direction; not a fixed celestial target';
    end
end
