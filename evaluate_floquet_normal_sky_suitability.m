function scan = evaluate_floquet_normal_sky_suitability( ...
    nInertial,sTargetsInertial,alphaThresholdsDeg,alphaUsefulDeg)
%EVALUATE_FLOQUET_NORMAL_SKY_SUITABILITY Sign-invariant normal suitability.
%
%   alpha(t) = acos(|s_target^T n_F(t)|)
%
% The continuous kernel is a raised cosine on [0,alphaUsefulDeg] and zero
% outside. The reported score is the equal-weight average over the supplied
% time/quadrature samples.

    narginchk(4,4);

    validateattributes(nInertial,{'numeric'},{'real','finite','2d'}, ...
        mfilename,'nInertial',1);
    validateattributes(sTargetsInertial,{'numeric'},{'real','finite','2d'}, ...
        mfilename,'sTargetsInertial',2);
    validateattributes(alphaThresholdsDeg,{'numeric'}, ...
        {'real','finite','vector','positive'},mfilename,'alphaThresholdsDeg',3);
    validateattributes(alphaUsefulDeg,{'numeric'}, ...
        {'real','finite','scalar','positive','<=',90},mfilename,'alphaUsefulDeg',4);

    if size(nInertial,1) ~= 3
        error('FINESST:NormalSuitability:NormalShape','nInertial must be 3xNt.');
    end
    if size(sTargetsInertial,2) ~= 3
        error('FINESST:NormalSuitability:TargetShape', ...
            'sTargetsInertial must be Ntarget x 3.');
    end

    alphaThresholdsDeg = alphaThresholdsDeg(:).';
    if any(diff(alphaThresholdsDeg) <= 0) || any(alphaThresholdsDeg > 90)
        error('FINESST:NormalSuitability:InvalidThresholds', ...
            'alphaThresholdsDeg must be strictly increasing and <= 90 deg.');
    end

    normalNormError = max(abs(vecnorm(nInertial,2,1)-1));
    targetNormError = max(abs(vecnorm(sTargetsInertial,2,2)-1));
    if normalNormError > 1e-10
        error('FINESST:NormalSuitability:NonunitNormal', ...
            'Input normals have maximum norm error %.3e.',normalNormError);
    end
    if targetNormError > 1e-10
        error('FINESST:NormalSuitability:NonunitTarget', ...
            'Input targets have maximum norm error %.3e.',targetNormError);
    end

    nTargets = size(sTargetsInertial,1);
    nThresholds = numel(alphaThresholdsDeg);
    nTime = size(nInertial,2);
    if nTime < 2
        error('FINESST:NormalSuitability:InsufficientTimeSamples', ...
            'At least two Floquet-normal samples are required.');
    end

    continuousScore = zeros(nTargets,1);
    fractionWithin = zeros(nTargets,nThresholds);
    alphaMinDeg = zeros(nTargets,1);
    alphaMeanDeg = zeros(nTargets,1);
    alphaRMSDeg = zeros(nTargets,1);
    targetChunkSize = 1000;

    for i0 = 1:targetChunkSize:nTargets
        i1 = min(i0+targetChunkSize-1,nTargets);
        S = sTargetsInertial(i0:i1,:);
        C = min(1,max(0,abs(S*nInertial)));
        alphaDeg = acosd(C);

        g = zeros(size(alphaDeg));
        inside = alphaDeg <= alphaUsefulDeg;
        g(inside) = 0.5*(1+cos(pi*alphaDeg(inside)/alphaUsefulDeg));

        continuousScore(i0:i1) = mean(g,2);
        alphaMinDeg(i0:i1) = min(alphaDeg,[],2);
        alphaMeanDeg(i0:i1) = mean(alphaDeg,2);
        alphaRMSDeg(i0:i1) = sqrt(mean(alphaDeg.^2,2));
        for j = 1:nThresholds
            fractionWithin(i0:i1,j) = mean(alphaDeg <= alphaThresholdsDeg(j),2);
        end
    end

    if any(~isfinite(continuousScore)) || any(~isfinite(fractionWithin(:))) || ...
            any(~isfinite(alphaMinDeg)) || any(~isfinite(alphaMeanDeg)) || ...
            any(~isfinite(alphaRMSDeg))
        error('FINESST:NormalSuitability:NonfiniteOutput', ...
            'Nonfinite values were produced by the suitability evaluation.');
    end

    scan = struct();
    scan.continuousScore = continuousScore;
    scan.fractionWithin = fractionWithin;
    scan.alphaThresholdsDeg = alphaThresholdsDeg;
    scan.alphaUsefulDeg = alphaUsefulDeg;
    scan.alphaMinDeg = alphaMinDeg;
    scan.alphaMeanDeg = alphaMeanDeg;
    scan.alphaRMSDeg = alphaRMSDeg;
    scan.metric = 'raised-cosine normal alignment';
    scan.signAmbiguity = true;
end
