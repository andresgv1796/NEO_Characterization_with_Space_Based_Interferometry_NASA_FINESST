function scan = evaluate_target_suitability_scan( ...
    geom, sys, targetSamples, lambda_m, maxOffNormalDeg, ...
    usePlaneNormalSignAmbiguity, useParallel)
%EVALUATE_TARGET_SUITABILITY_SCAN Evaluate accessibility and UV quality.
%
% scan = evaluate_target_suitability_scan(geom, sys, targetSamples, ...
%     lambda_m, maxOffNormalDeg, usePlaneNormalSignAmbiguity, useParallel)
%
% Inputs:
%   geom                         : geometry struct from QPO pipeline
%   sys                          : CR3BP system struct
%   targetSamples                : output of sample_sphere_targets
%   lambda_m                     : wavelength [m]
%   maxOffNormalDeg              : accessibility half-angle [deg]
%   usePlaneNormalSignAmbiguity  : true uses abs(dot(n,s))
%   useParallel                  : true attempts parfor
%
% Output:
%   scan : struct with one row per sampled target.

    if nargin < 7 || isempty(useParallel)
        useParallel = false;
    end

    if nargin < 6 || isempty(usePlaneNormalSignAmbiguity)
        usePlaneNormalSignAmbiguity = true;
    end

    if nargin < 5 || isempty(maxOffNormalDeg)
        maxOffNormalDeg = 30.0;
    end

    if nargin < 4 || isempty(lambda_m)
        error('lambda_m must be provided.');
    end

    if nargin < 3 || isempty(targetSamples)
        targetSamples = sample_sphere_targets(300, 'fibonacci');
    end

    if ~isstruct(targetSamples) || ~isfield(targetSamples, 'sInertial')
        error('targetSamples must be the struct returned by sample_sphere_targets.');
    end

    sList = targetSamples.sInertial;
    nTargets = size(sList, 1);

    accessibleFraction = zeros(nTargets, 1);
    thetaMinDeg = zeros(nTargets, 1);
    thetaMeanDeg = zeros(nTargets, 1);
    thetaMedianDeg = zeros(nTargets, 1);
    uvAnisotropy = NaN(nTargets, 1);
    uvMax_Glambda = NaN(nTargets, 1);
    uvRms_Glambda = NaN(nTargets, 1);
    uvHullArea_Glambda2 = NaN(nTargets, 1);
    nUvSamples = zeros(nTargets, 1);
    suitabilityScore = zeros(nTargets, 1);

    runParallel = false;
    if useParallel && exist('gcp', 'file') == 2
        try
            if isempty(gcp('nocreate'))
                parpool;
            end
            runParallel = true;
        catch
            warning('Parallel pool could not be started. Falling back to serial loop.');
            runParallel = false;
        end
    end

    if runParallel
        parfor k = 1:nTargets
            [accessibleFraction(k), thetaMinDeg(k), thetaMeanDeg(k), ...
                thetaMedianDeg(k), uvAnisotropy(k), uvMax_Glambda(k), ...
                uvRms_Glambda(k), uvHullArea_Glambda2(k), nUvSamples(k), ...
                suitabilityScore(k)] = evaluate_one_target( ...
                geom, sys, lambda_m, maxOffNormalDeg, ...
                usePlaneNormalSignAmbiguity, sList(k,:));
        end
    else
        for k = 1:nTargets
            [accessibleFraction(k), thetaMinDeg(k), thetaMeanDeg(k), ...
                thetaMedianDeg(k), uvAnisotropy(k), uvMax_Glambda(k), ...
                uvRms_Glambda(k), uvHullArea_Glambda2(k), nUvSamples(k), ...
                suitabilityScore(k)] = evaluate_one_target( ...
                geom, sys, lambda_m, maxOffNormalDeg, ...
                usePlaneNormalSignAmbiguity, sList(k,:));
        end
    end

    scan = struct();
    scan.targets = targetSamples;
    scan.raDeg = targetSamples.raDeg(:);
    scan.decDeg = targetSamples.decDeg(:);
    scan.sInertial = sList;

    scan.accessibleFraction = accessibleFraction;
    scan.thetaMinDeg = thetaMinDeg;
    scan.thetaMeanDeg = thetaMeanDeg;
    scan.thetaMedianDeg = thetaMedianDeg;

    scan.uvAnisotropy = uvAnisotropy;
    scan.uvMax_Glambda = uvMax_Glambda;
    scan.uvRms_Glambda = uvRms_Glambda;
    scan.uvHullArea_Glambda2 = uvHullArea_Glambda2;
    scan.nUvSamples = nUvSamples;

    scan.suitabilityScore = suitabilityScore;
    scan.score = suitabilityScore; % short alias for convenience

    scan.maxOffNormalDeg = maxOffNormalDeg;
    scan.usePlaneNormalSignAmbiguity = usePlaneNormalSignAmbiguity;
    scan.lambda_m = lambda_m;
    scan.usedParallel = runParallel;

    [~, bestIdx] = max(suitabilityScore);
    scan.bestIndex = bestIdx;
    scan.bestRaDeg = scan.raDeg(bestIdx);
    scan.bestDecDeg = scan.decDeg(bestIdx);
    scan.bestTargetVector = scan.sInertial(bestIdx,:);
end

function [fAcc, thetaMin, thetaMean, thetaMedian, etaUv, rhoMaxG, rhoRmsG, areaG2, nUv, score] = ...
    evaluate_one_target(geom, sys, lambda_m, maxOffNormalDeg, useAmbiguity, sTarget)

    target = make_target_definition( ...
        'inertial', false, NaN, NaN, sTarget, [1, 0, 0], [0, 0, 1]);

    access = compute_triangle_accessibility( ...
        geom, target, maxOffNormalDeg, useAmbiguity);

    uv = compute_fixed_target_uv(geom, target, sys, lambda_m);
    m = compute_uv_suitability_metrics(uv, access.accessible);

    theta = access.thetaOffDeg(:);
    fAcc = access.accessibleFraction;
    thetaMin = min(theta, [], 'omitnan');
    thetaMean = mean(theta, 'omitnan');
    thetaMedian = median(theta, 'omitnan');

    etaUv = m.uvAnisotropy;
    rhoMaxG = m.uvMax_Glambda;
    rhoRmsG = m.uvRms_Glambda;
    areaG2 = m.uvHullArea_Glambda2;
    nUv = m.nUvSamples;

    if ~isfinite(etaUv)
        etaForScore = 0.0;
    else
        etaForScore = etaUv;
    end

    score = fAcc*etaForScore;
end
