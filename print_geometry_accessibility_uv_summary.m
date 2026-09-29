function print_geometry_accessibility_uv_summary(geom, access, uv)
%PRINT_GEOMETRY_ACCESSIBILITY_UV_SUMMARY Print compact diagnostics.

    fprintf('\nThree-collector QPO geometry summary\n');
    fprintf('  time span = %.3f primary periods\n', geom.tau(end));

    if isfield(geom, 'dim')
        fprintf('  dimensional time span = %.6g %s\n', geom.dim.timePlot(end), geom.dim.timePlotUnit);
        fprintf('  system = %s, L* = %.6g km, t* = %.6g days\n', ...
            geom.dim.system.name, geom.dim.system.Lstar_km, geom.dim.system.Tstar_days);
    end

    fprintf('\nCollector-collector baseline ranges [nondim.]\n');
    fprintf('  L12: [%.6e, %.6e]\n', min(geom.baselines.L12), max(geom.baselines.L12));
    fprintf('  L13: [%.6e, %.6e]\n', min(geom.baselines.L13), max(geom.baselines.L13));
    fprintf('  L23: [%.6e, %.6e]\n', min(geom.baselines.L23), max(geom.baselines.L23));

    if isfield(geom, 'dim')
        fprintf('Collector-collector baseline ranges [km]\n');
        fprintf('  L12: [%.6e, %.6e]\n', min(geom.dim.baselines.L12_km), max(geom.dim.baselines.L12_km));
        fprintf('  L13: [%.6e, %.6e]\n', min(geom.dim.baselines.L13_km), max(geom.dim.baselines.L13_km));
        fprintf('  L23: [%.6e, %.6e]\n', min(geom.dim.baselines.L23_km), max(geom.dim.baselines.L23_km));
    end

    fprintf('\nCombiner-centroid offset [nondim.]\n');
    fprintf('  ||centroid - combiner|| min/mean/max = %.6e / %.6e / %.6e\n', ...
        min(geom.centroidOffset), mean(geom.centroidOffset), max(geom.centroidOffset));
    fprintf('  normalized offset min/mean/max = %.6e / %.6e / %.6e\n', ...
        min(geom.normalizedCentroidOffset), mean(geom.normalizedCentroidOffset), ...
        max(geom.normalizedCentroidOffset));

    if isfield(geom, 'dim')
        fprintf('Combiner-centroid offset [km]\n');
        fprintf('  min/mean/max = %.6e / %.6e / %.6e\n', ...
            min(geom.dim.centroidOffset_km), mean(geom.dim.centroidOffset_km), max(geom.dim.centroidOffset_km));
    end

    fprintf('\nTriangle area and normal conditioning\n');
    fprintf('  area [nondim.^2] min/mean/max = %.6e / %.6e / %.6e\n', ...
        min(geom.area), mean(geom.area), max(geom.area));
    if isfield(geom, 'dim')
        fprintf('  area [km^2]      min/mean/max = %.6e / %.6e / %.6e\n', ...
            min(geom.dim.area_km2), mean(geom.dim.area_km2), max(geom.dim.area_km2));
    end
    fprintf('  sin(angle b12,b13) min/mean/max = %.6e / %.6e / %.6e\n', ...
        min(geom.sinAngle12_13), mean(geom.sinAngle12_13), max(geom.sinAngle12_13));
    fprintf('  degenerate samples = %d of %d\n', nnz(geom.degenerate), numel(geom.degenerate));

    fprintf('\nTarget accessibility\n');
    fprintf('  max off-normal angle = %.3f deg\n', access.thetaMaxDeg);
    fprintf('  off-normal angle min/mean/max = %.6f / %.6f / %.6f deg\n', ...
        min(access.thetaOffDeg), mean(access.thetaOffDeg), max(access.thetaOffDeg));
    fprintf('  accessible fraction = %.2f %%\n', 100*access.accessibleFraction);

    fprintf('\nFixed-target UV coverage summary\n');
    fprintf('  frame = %s\n', uv.frameDescription);
    if isfield(uv, 'basisConvention')
        fprintf('  basis = %s\n', uv.basisConvention);
    end
    fprintf('  units = %s, scale = %.6e\n', uv.units, uv.scale);
    if isfield(uv, 'lambda_m') && isfinite(uv.lambda_m)
        fprintf('  lambda = %.6e m, L*/lambda = %.6e\n', uv.lambda_m, uv.scale);
    end
    fprintf('  rho12 min/mean/max = %.6e / %.6e / %.6e\n', ...
        min(uv.rho12), mean(uv.rho12), max(uv.rho12));
    fprintf('  rho13 min/mean/max = %.6e / %.6e / %.6e\n', ...
        min(uv.rho13), mean(uv.rho13), max(uv.rho13));
    fprintf('  rho23 min/mean/max = %.6e / %.6e / %.6e\n', ...
        min(uv.rho23), mean(uv.rho23), max(uv.rho23));
end
