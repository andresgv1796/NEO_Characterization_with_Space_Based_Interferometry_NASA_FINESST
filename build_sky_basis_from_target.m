function [eU, eV, eW] = build_sky_basis_from_target(sHat, northReference)
%BUILD_SKY_BASIS_FROM_TARGET Build a right-handed target-centered basis.
%
% This helper is retained for non-RA/Dec diagnostic cases. For standard
% fixed celestial targets, prefer radec_to_uvw_basis or convert the target
% vector to RA/Dec first with unit_to_radec.
%
% eW points to the target. eV is the projection of northReference into the
% target sky plane. eU is chosen so that eU x eV = eW.

    eW = normalize_row(sHat);
    ref1 = normalize_row(northReference);
    ref2 = [1, 0, 0];
    ref3 = [0, 1, 0];

    eV = ref1 - dot(ref1, eW) * eW;
    if norm(eV) < 1e-10
        eV = ref2 - dot(ref2, eW) * eW;
    end
    if norm(eV) < 1e-10
        eV = ref3 - dot(ref3, eW) * eW;
    end
    eV = normalize_row(eV);

    % Right-handed interferometry convention: east x north = source.
    eU = cross(eV, eW);
    eU = normalize_row(eU);

    % Re-orthogonalize eV to limit roundoff while preserving handedness.
    eV = cross(eW, eU);
    eV = normalize_row(eV);
end
