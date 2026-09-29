function uv = compute_fixed_target_uv(geom, target, scaleOrSystem, lambda)
%COMPUTE_FIXED_TARGET_UV Compute conventional fixed-target UV coordinates.
%
% Preferred call:
%   uv = compute_fixed_target_uv(geom, target, sys, lambda_m)
%
% Backward-compatible call:
%   uv = compute_fixed_target_uv(geom, target, lengthScale, lambda)

    if nargin < 4
        error('compute_fixed_target_uv requires geom, target, scaleOrSystem, and lambda.');
    end

    if isstruct(scaleOrSystem)
        sys = scaleOrSystem;
        lambda_m = lambda;
        if lambda_m <= 0
            error('lambda_m must be positive.');
        end
        lengthScale_km = sys.Lstar_km;
        lambda_km = lambda_m * 1e-3;
        scale = lengthScale_km / lambda_km;
        unitDescription = 'wavelengths';
    else
        lengthScale = scaleOrSystem;
        if lambda <= 0
            error('lambda must be positive.');
        end
        scale = lengthScale / lambda;
        lengthScale_km = NaN;
        lambda_m = NaN;
        lambda_km = NaN;
        unitDescription = 'baseline/lambda';
    end

    t = geom.t;
    nGrid = numel(t);

    b12S = geom.baselines.b12;
    b13S = geom.baselines.b13;
    b23S = geom.baselines.b23;

    switch lower(target.frame)
        case 'inertial'
            eU = repmat(target.eUInertial0, nGrid, 1);
            eV = repmat(target.eVInertial0, nGrid, 1);
            eW = repmat(target.eWInertial0, nGrid, 1);

            b12 = rotate_synodic_to_inertial_rows(b12S, t);
            b13 = rotate_synodic_to_inertial_rows(b13S, t);
            b23 = rotate_synodic_to_inertial_rows(b23S, t);
            frameDescription = 'fixed inertial target sky basis';

        case 'synodic'
            eU = repmat(target.eUSynodic0, nGrid, 1);
            eV = repmat(target.eVSynodic0, nGrid, 1);
            eW = repmat(target.eWSynodic0, nGrid, 1);

            b12 = b12S;
            b13 = b13S;
            b23 = b23S;
            frameDescription = 'fixed synodic CR3BP target basis';

        otherwise
            error('Unknown target frame.');
    end

    uvw12 = [sum(b12 .* eU, 2), sum(b12 .* eV, 2), sum(b12 .* eW, 2)] * scale;
    uvw13 = [sum(b13 .* eU, 2), sum(b13 .* eV, 2), sum(b13 .* eW, 2)] * scale;
    uvw23 = [sum(b23 .* eU, 2), sum(b23 .* eV, 2), sum(b23 .* eW, 2)] * scale;

    uvAll = [uvw12(:,1:2); uvw13(:,1:2); uvw23(:,1:2); ...
            -uvw12(:,1:2); -uvw13(:,1:2); -uvw23(:,1:2)];

    uv = struct();
    uv.uvw12 = uvw12;
    uv.uvw13 = uvw13;
    uv.uvw23 = uvw23;
    uv.uvAll = uvAll;
    uv.rho12 = vecnorm(uvw12(:,1:2), 2, 2);
    uv.rho13 = vecnorm(uvw13(:,1:2), 2, 2);
    uv.rho23 = vecnorm(uvw23(:,1:2), 2, 2);
    uv.w12 = uvw12(:,3);
    uv.w13 = uvw13(:,3);
    uv.w23 = uvw23(:,3);
    uv.eU = eU;
    uv.eV = eV;
    uv.eW = eW;
    uv.scale = scale;
    uv.lengthScale_km = lengthScale_km;
    uv.lambda_m = lambda_m;
    uv.lambda_km = lambda_km;
    uv.units = unitDescription;
    uv.frameDescription = frameDescription;
    if isfield(target, 'basisConvention')
        uv.basisConvention = target.basisConvention;
    else
        uv.basisConvention = '';
    end
end
