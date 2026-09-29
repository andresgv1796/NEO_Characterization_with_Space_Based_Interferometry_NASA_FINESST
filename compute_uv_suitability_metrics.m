function metrics = compute_uv_suitability_metrics(uv, accessibleMask)
%COMPUTE_UV_SUITABILITY_METRICS Scalar UV-quality metrics for one target.
%
% metrics = compute_uv_suitability_metrics(uv, accessibleMask)
%
% The UV coordinates are assumed to be in wavelengths, as produced by
% compute_fixed_target_uv. Metrics are computed using accessible samples
% only. Hermitian points are included for coverage-area/anisotropy metrics.

    if nargin < 2 || isempty(accessibleMask)
        accessibleMask = true(size(uv.uvw12, 1), 1);
    end
    accessibleMask = logical(accessibleMask(:));

    pts = [uv.uvw12(accessibleMask, 1:2); ...
           uv.uvw13(accessibleMask, 1:2); ...
           uv.uvw23(accessibleMask, 1:2)];

    pts = pts(all(isfinite(pts), 2), :);
    pts = [pts; -pts];

    metrics = struct();
    metrics.nUvSamples = size(pts, 1);
    metrics.uvMax_wavelengths = NaN;
    metrics.uvRms_wavelengths = NaN;
    metrics.uvMax_Glambda = NaN;
    metrics.uvRms_Glambda = NaN;
    metrics.uvAnisotropy = NaN;
    metrics.uvHullArea_Glambda2 = NaN;

    if isempty(pts)
        return;
    end

    rho = sqrt(sum(pts.^2, 2));
    metrics.uvMax_wavelengths = max(rho);
    metrics.uvRms_wavelengths = sqrt(mean(rho.^2));
    metrics.uvMax_Glambda = metrics.uvMax_wavelengths/1e9;
    metrics.uvRms_Glambda = metrics.uvRms_wavelengths/1e9;

    if size(pts, 1) >= 3
        centered = pts - mean(pts, 1);
        C = cov(centered);
        eigVals = sort(real(eig(C)), 'descend');
        if eigVals(1) > 0
            metrics.uvAnisotropy = sqrt(max(eigVals(2), 0)/eigVals(1));
        end

        if rank(centered, 1e-10*max(1, norm(centered, 'fro'))) >= 2
            try
                idx = convhull(pts(:,1), pts(:,2));
                areaWavelength2 = polyarea(pts(idx,1), pts(idx,2));
                metrics.uvHullArea_Glambda2 = areaWavelength2/1e18;
            catch
                metrics.uvHullArea_Glambda2 = 0.0;
            end
        else
            metrics.uvHullArea_Glambda2 = 0.0;
        end
    end
end
