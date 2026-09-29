function access = compute_triangle_accessibility( ...
    geom, target, maxOffNormalDeg, usePlaneNormalSignAmbiguity)
%COMPUTE_TRIANGLE_ACCESSIBILITY Compare triangle normal to target over time.

    t = geom.t;
    nS = geom.normalSynodic;
    nI = geom.normalInertial;
    nGrid = numel(t);

    switch lower(target.frame)
        case 'inertial'
            sI = repmat(target.sInertial0, nGrid, 1);
            sS = rotate_inertial_to_synodic_rows(sI, t);
        case 'synodic'
            sS = repmat(target.sSynodic0, nGrid, 1);
            sI = rotate_synodic_to_inertial_rows(sS, t);
        otherwise
            error('Unknown target frame.');
    end

    dotSyn = sum(nS .* sS, 2);
    dotInertial = sum(nI .* sI, 2);

    if usePlaneNormalSignAmbiguity
        cosOff = abs(dotSyn);
    else
        cosOff = dotSyn;
    end
    cosOff = max(-1, min(1, cosOff));
    offNormalDeg = acosd(cosOff);
    thetaOffDeg = offNormalDeg;
    accessible = thetaOffDeg <= maxOffNormalDeg;

    access = struct();
    access.sTargetSynodic = sS;
    access.sTargetInertial = sI;
    access.dotSynodic = dotSyn;
    access.dotInertial = dotInertial;
    access.offNormalDeg = offNormalDeg;        % backward-compatible field name
    access.thetaOffDeg = thetaOffDeg;      % preferred: not confused with RA alpha
    access.maxOffNormalDeg = maxOffNormalDeg;
    access.thetaMaxDeg = maxOffNormalDeg;
    access.accessible = accessible;
    access.accessibleFraction = nnz(accessible) / numel(accessible);
    access.usePlaneNormalSignAmbiguity = usePlaneNormalSignAmbiguity;
end
