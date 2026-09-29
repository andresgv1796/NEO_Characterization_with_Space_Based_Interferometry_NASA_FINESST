function geom = compute_three_collector_qpo_geometry( ...
    mu, xCombiner0, TPrimary, XCollectors0, ...
    nPrimaryPeriods, nSavePerPrimaryPeriod, odeOpts, normalReference)
%COMPUTE_THREE_COLLECTOR_QPO_GEOMETRY Propagate and compute triangle geometry.

    if size(XCollectors0, 1) ~= 6 || size(XCollectors0, 2) ~= 3
        error('XCollectors0 must be 6 x 3.');
    end

    tEnd = nPrimaryPeriods * TPrimary;
    nGrid = nPrimaryPeriods * nSavePerPrimaryPeriod + 1;
    t = linspace(0, tEnd, nGrid).';
    tau = t / TPrimary;

    YC = propagate_state_on_grid(xCombiner0, t, mu, odeOpts);

    YCollectors = zeros(nGrid, 6, 3);
    for i = 1:3
        YCollectors(:, :, i) = propagate_state_on_grid(XCollectors0(:, i), t, mu, odeOpts);
    end

    rC = YC(:, 1:3);
    vC = YC(:, 4:6);

    r = zeros(nGrid, 3, 3);
    v = zeros(nGrid, 3, 3);
    for i = 1:3
        r(:, :, i) = YCollectors(:, 1:3, i);
        v(:, :, i) = YCollectors(:, 4:6, i);
    end

    b12 = r(:, :, 2) - r(:, :, 1);
    b13 = r(:, :, 3) - r(:, :, 1);
    b23 = r(:, :, 3) - r(:, :, 2);

    L12 = vecnorm(b12, 2, 2);
    L13 = vecnorm(b13, 2, 2);
    L23 = vecnorm(b23, 2, 2);

    bC1 = r(:, :, 1) - rC;
    bC2 = r(:, :, 2) - rC;
    bC3 = r(:, :, 3) - rC;

    LC1 = vecnorm(bC1, 2, 2);
    LC2 = vecnorm(bC2, 2, 2);
    LC3 = vecnorm(bC3, 2, 2);

    centroid = mean(r, 3);
    centroidMinusCombiner = centroid - rC;
    centroidOffset = vecnorm(centroidMinusCombiner, 2, 2);

    collectorMinusCentroid = zeros(nGrid, 3, 3);
    collectorRadius = zeros(nGrid, 3);
    for i = 1:3
        collectorMinusCentroid(:, :, i) = r(:, :, i) - centroid;
        collectorRadius(:, i) = vecnorm(collectorMinusCentroid(:, :, i), 2, 2);
    end

    rmsCollectorRadius = sqrt(mean(collectorRadius.^2, 2));
    normalizedCentroidOffset = centroidOffset ./ max(rmsCollectorRadius, eps);

    LC = [LC1, LC2, LC3];
    combinerRangeMean = mean(LC, 2);
    combinerRangeSpread = max(LC, [], 2) - min(LC, [], 2);
    normalizedCombinerRangeSpread = combinerRangeSpread ./ max(combinerRangeMean, eps);

    normalCross = cross(b12, b13, 2);
    twoArea = vecnorm(normalCross, 2, 2);
    area = 0.5 * twoArea;

    normalRaw = normalCross ./ max(twoArea, eps);
    degenerate = twoArea < 1e-12;
    normalRaw(degenerate, :) = NaN;
    normal = enforce_normal_sign_continuity(normalRaw, normalReference);

    sinAngle12_13 = twoArea ./ max(L12 .* L13, eps);

    dCombinerFromCentroid = rC - centroid;
    signedPlaneOffset = sum(dCombinerFromCentroid .* normal, 2);
    absPlaneOffset = abs(signedPlaneOffset);
    inPlaneOffsetVec = dCombinerFromCentroid - signedPlaneOffset .* normal;
    inPlaneOffset = vecnorm(inPlaneOffsetVec, 2, 2);
    normalizedInPlaneOffset = inPlaneOffset ./ max(rmsCollectorRadius, eps);
    normalizedAbsPlaneOffset = absPlaneOffset ./ max(rmsCollectorRadius, eps);

    normalAzSynodic = unwrap(atan2(normal(:, 2), normal(:, 1)));
    normalElSynodic = asin(max(-1, min(1, normal(:, 3))));

    normalInertial = rotate_synodic_to_inertial_rows(normal, t);
    normalAzInertial = unwrap(atan2(normalInertial(:, 2), normalInertial(:, 1)));
    normalElInertial = asin(max(-1, min(1, normalInertial(:, 3))));

    geom = struct();
    geom.t = t;
    geom.tau = tau;
    geom.YC = YC;
    geom.YCollectors = YCollectors;
    geom.rC = rC;
    geom.vC = vC;
    geom.r = r;
    geom.v = v;

    geom.baselines.b12 = b12;
    geom.baselines.b13 = b13;
    geom.baselines.b23 = b23;
    geom.baselines.L12 = L12;
    geom.baselines.L13 = L13;
    geom.baselines.L23 = L23;
    geom.baselines.bC1 = bC1;
    geom.baselines.bC2 = bC2;
    geom.baselines.bC3 = bC3;
    geom.baselines.LC1 = LC1;
    geom.baselines.LC2 = LC2;
    geom.baselines.LC3 = LC3;

    geom.centroid = centroid;
    geom.centroidMinusCombiner = centroidMinusCombiner;
    geom.centroidOffset = centroidOffset;
    geom.normalizedCentroidOffset = normalizedCentroidOffset;
    geom.collectorMinusCentroid = collectorMinusCentroid;
    geom.collectorRadius = collectorRadius;
    geom.rmsCollectorRadius = rmsCollectorRadius;

    geom.dCombinerFromCentroid = dCombinerFromCentroid;
    geom.signedPlaneOffset = signedPlaneOffset;
    geom.absPlaneOffset = absPlaneOffset;
    geom.inPlaneOffsetVec = inPlaneOffsetVec;
    geom.inPlaneOffset = inPlaneOffset;
    geom.normalizedInPlaneOffset = normalizedInPlaneOffset;
    geom.normalizedAbsPlaneOffset = normalizedAbsPlaneOffset;
    geom.combinerRangeMean = combinerRangeMean;
    geom.combinerRangeSpread = combinerRangeSpread;
    geom.normalizedCombinerRangeSpread = normalizedCombinerRangeSpread;

    geom.area = area;
    geom.twoArea = twoArea;
    geom.normalCross = normalCross;
    geom.normalRaw = normalRaw;
    geom.normalSynodic = normal;
    geom.normalInertial = normalInertial;
    geom.normalAzSynodic = normalAzSynodic;
    geom.normalElSynodic = normalElSynodic;
    geom.normalAzInertial = normalAzInertial;
    geom.normalElInertial = normalElInertial;
    geom.sinAngle12_13 = sinAngle12_13;
    geom.degenerate = degenerate;
end
