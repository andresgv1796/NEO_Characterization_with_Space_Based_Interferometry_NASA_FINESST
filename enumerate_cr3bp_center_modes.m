function spectrum = enumerate_cr3bp_center_modes(PhiT,TPrimary,xPO0,mu,opts)
%ENUMERATE_CR3BP_CENTER_MODES Enumerate every nontrivial center eigenpair.
%
% spectrum = enumerate_cr3bp_center_modes(PhiT,TPrimary,xPO0,mu,opts)
%
% The routine makes no planar/vertical branch request. Every positive-
% imaginary, nontrivial, near-unit-circle multiplier is returned once, with
% its conjugate understood. Modes are ordered only by Floquet angle.
%
% Required opts fields
%   nTrivialMultipliers       number of autonomous +1 multipliers (2 in CR3BP)
%   trivialPairMaxDistance    quality threshold for the identified trivial pair
%   trivialFlowAlignmentMin   quality threshold for the flow eigenvector
%   trivialMultiplierTol     exclusion radius around lambda = +1
%   imagTol                   positive-imaginary threshold
%   centerRadialLogTol        admissible abs(log(abs(lambda)))
%   centerConjugacyTol        admissible conjugate-pair mismatch
%   stabilityRadialLogTol     radial tolerance used for stability counts
%   complexQuartetAngleTol    equal-angle test for a complex quartet
%   complexQuartetRadialTol   minimum radial splitting of a complex quartet
%   eigenpairResidualRelTol  maximum normalized eigenpair residual

% Selected outputs
%   spectrum.lambdaAll
%   spectrum.centerModes(k).lambda
%   spectrum.centerModes(k).rho
%   spectrum.centerModes(k).eigenvector
%   spectrum.centerModes(k).basis6
%   spectrum.centerModes(k).planarProjection
%   spectrum.centerModes(k).verticalProjection
%   spectrum.featureVector
%   spectrum.stabilitySignature

% No unit-circle projection, eigenvector substitution, or requested-mode
% fallback is performed.

    narginchk(5,5);

    validateattributes(PhiT,{'numeric'}, ...
        {'finite','2d','size',[6 6]},mfilename,'PhiT',1);
    validateattributes(TPrimary,{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,'TPrimary',2);
    validateattributes(xPO0,{'numeric'}, ...
        {'real','finite','vector','numel',6},mfilename,'xPO0',3);
    validateattributes(mu,{'numeric'}, ...
        {'real','finite','scalar','>',0,'<',0.5},mfilename,'mu',4);
    xPO0 = xPO0(:);

    requiredFields = { ...
        'nTrivialMultipliers', ...
        'trivialPairMaxDistance', ...
        'trivialFlowAlignmentMin', ...
        'trivialMultiplierTol', ...
        'imagTol', ...
        'centerRadialLogTol', ...
        'centerConjugacyTol', ...
        'stabilityRadialLogTol', ...
        'complexQuartetAngleTol', ...
        'complexQuartetRadialTol', ...
        'eigenpairResidualRelTol'};

    if ~isstruct(opts) || ~all(isfield(opts,requiredFields))
        error('FINESST:CenterModes:InvalidOptions', ...
            'opts must explicitly contain: %s.',strjoin(requiredFields,', '));
    end

    validateattributes(opts.nTrivialMultipliers,{'numeric'}, ...
        {'real','finite','scalar','integer','positive','<',6}, ...
        mfilename,'opts.nTrivialMultipliers');
    if opts.nTrivialMultipliers ~= 2
        error('FINESST:CenterModes:InvalidTrivialDimension', ...
            'The autonomous CR3BP monodromy requires exactly two trivial multipliers.');
    end
    validateattributes(opts.trivialPairMaxDistance,{'numeric'}, ...
        {'real','finite','scalar','positive'}, ...
        mfilename,'opts.trivialPairMaxDistance');
    validateattributes(opts.trivialFlowAlignmentMin,{'numeric'}, ...
        {'real','finite','scalar','>=',0,'<=',1}, ...
        mfilename,'opts.trivialFlowAlignmentMin');
    validateattributes(opts.trivialMultiplierTol,{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,'opts.trivialMultiplierTol');
    validateattributes(opts.imagTol,{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,'opts.imagTol');
    validateattributes(opts.centerRadialLogTol,{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,'opts.centerRadialLogTol');
    validateattributes(opts.centerConjugacyTol,{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,'opts.centerConjugacyTol');
    validateattributes(opts.stabilityRadialLogTol,{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,'opts.stabilityRadialLogTol');
    validateattributes(opts.complexQuartetAngleTol,{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,'opts.complexQuartetAngleTol');
    validateattributes(opts.complexQuartetRadialTol,{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,'opts.complexQuartetRadialTol');
    validateattributes(opts.eigenpairResidualRelTol,{'numeric'}, ...
        {'real','finite','scalar','positive'}, ...
        mfilename,'opts.eigenpairResidualRelTol');

    [V,D] = eig(PhiT);
    lambdaAll = diag(D);

    if any(~isfinite(lambdaAll)) || any(~isfinite(V(:)))
        error('FINESST:CenterModes:NonfiniteEigendecomposition', ...
            'The monodromy eigendecomposition contains nonfinite values.');
    end

    modulus = abs(lambdaAll);
    phase = angle(lambdaAll);

    % Identify the autonomous flow multiplier from its eigenvector, then pair
    % it with the remaining multiplier closest to +1. This prevents a genuine
    % nontrivial mode near +1 from being removed merely because eig() ordering
    % or radial noise changes along the family.
    flowTangent = local_cr3bp_vector_field(xPO0,mu);
    flowNorm = norm(flowTangent);
    if ~(isfinite(flowNorm) && flowNorm > 0)
        error('FINESST:CenterModes:InvalidFlowTangent', ...
            'The periodic-orbit flow tangent is invalid.');
    end

    flowAlignment = zeros(6,1);
    for i = 1:6
        flowAlignment(i) = abs(V(:,i)'*flowTangent)/(norm(V(:,i))*flowNorm);
    end
    [trivialFlowAlignment,iFlow] = max(flowAlignment);

    remaining = setdiff((1:6).',iFlow,'stable');
    [~,iPartnerLocal] = min(abs(lambdaAll(remaining)-1));
    iPartner = remaining(iPartnerLocal);
    trivialIdx = [iFlow;iPartner];
    trivialDistances = abs(lambdaAll(trivialIdx)-1);
    trivialPairResolved = ...
        max(trivialDistances) <= opts.trivialPairMaxDistance && ...
        trivialFlowAlignment >= opts.trivialFlowAlignmentMin;

    trivialPairMask = false(6,1);
    trivialPairMask(trivialIdx) = true;
    nearTrivialMask = abs(lambdaAll-1) <= opts.trivialMultiplierTol;
    nontrivialMask = ~trivialPairMask;

    radialLog = log(max(modulus,realmin));
    positiveComplexIdx = find( ...
        nontrivialMask & ~nearTrivialMask & imag(lambdaAll) > opts.imagTol);
    negativeComplexIdx = find( ...
        nontrivialMask & ~nearTrivialMask & imag(lambdaAll) < -opts.imagTol);

    candidateIdx = zeros(0,1);
    conjugacyDefect = zeros(0,1);
    for k = 1:numel(positiveComplexIdx)
        idx = positiveComplexIdx(k);
        if isempty(negativeComplexIdx)
            defect = inf;
        else
            defect = min(abs(lambdaAll(negativeComplexIdx)-conj(lambdaAll(idx)))) / ...
                max(1,abs(lambdaAll(idx)));
        end

        if abs(radialLog(idx)) <= opts.centerRadialLogTol && ...
                defect <= opts.centerConjugacyTol
            candidateIdx(end+1,1) = idx; %#ok<AGROW>
            conjugacyDefect(end+1,1) = defect; %#ok<AGROW>
        end
    end

    % Two positive-imaginary multipliers with the same angle and reciprocal
    % radii represent a complex-unstable quartet, not two center pairs.
    complexQuartetDetected = false;
    if numel(candidateIdx) == 2
        rhoCandidate = phase(candidateIdx);
        sigmaCandidate = radialLog(candidateIdx);
        sameAngle = abs(rhoCandidate(1)-rhoCandidate(2)) <= ...
            opts.complexQuartetAngleTol;
        reciprocalRadii = abs(sum(sigmaCandidate)) <= opts.centerConjugacyTol;
        resolvedRadialSplit = max(abs(sigmaCandidate)) >= ...
            opts.complexQuartetRadialTol;
        if sameAngle && reciprocalRadii && resolvedRadialSplit
            candidateIdx = zeros(0,1);
            conjugacyDefect = zeros(0,1);
            complexQuartetDetected = true;
        end
    end

    [~,iOrder] = sort(phase(candidateIdx),'ascend');
    candidateIdx = candidateIdx(iOrder);
    conjugacyDefect = conjugacyDefect(iOrder);

    emptyMode = struct( ...
        'localIndex',{}, ...
        'eigenvalueIndex',{}, ...
        'lambda',{}, ...
        'rho',{}, ...
        'floquetExponent',{}, ...
        'unitCircleDefect',{}, ...
        'radialLogDefect',{}, ...
        'conjugacyDefect',{}, ...
        'eigenpairResidualRel',{}, ...
        'eigenvector',{}, ...
        'basis6',{}, ...
        'planarProjection',{}, ...
        'verticalProjection',{}, ...
        'positionRankRatioAtSection',{}, ...
        'continuationSegmentId',{}, ...
        'segmentBranchId',{}, ...
        'branchId',{});

    centerModes = emptyMode;

    for k = 1:numel(candidateIdx)
        idx = candidateIdx(k);
        lambda = lambdaAll(idx);
        w = V(:,idx);

        nw = norm(w);
        if ~(isfinite(nw) && nw > 0)
            error('FINESST:CenterModes:InvalidEigenvector', ...
                'Center candidate %d has an invalid eigenvector norm.',k);
        end
        w = w/nw;

        % Deterministic phase gauge for saved diagnostics. The physical real
        % two-plane is unchanged by this choice.
        [~,iPivot] = max(abs(w));
        w = w*exp(-1i*angle(w(iPivot)));
        if real(w(iPivot)) < 0
            w = -w;
        end

        residualRel = norm(PhiT*w-lambda*w)/(norm(PhiT)*norm(w));
        if ~(isfinite(residualRel) && ...
                residualRel <= opts.eigenpairResidualRelTol)
            error('FINESST:CenterModes:PoorEigenpairResidual', ...
                ['Center candidate %d has normalized eigenpair residual ' ...
                 '%.6e (allowed %.6e).'], ...
                k,residualRel,opts.eigenpairResidualRelTol);
        end

        B6 = [real(w),imag(w)];
        [Q6,R6] = qr(B6,0);
        if size(Q6,2) ~= 2 || rank(R6) ~= 2
            error('FINESST:CenterModes:DefectiveRealCenterPlane', ...
                'Center candidate %d does not define a real 2-D state subspace.',k);
        end

        Apos0 = B6(1:3,:);
        sPos0 = svd(Apos0,'econ');
        if numel(sPos0) == 2 && sPos0(1) > 0
            positionRankRatioAtSection = sPos0(2)/sPos0(1);
        else
            positionRankRatioAtSection = 0;
        end

        mode = struct();
        mode.localIndex = k;
        mode.eigenvalueIndex = idx;
        mode.lambda = lambda;
        mode.rho = angle(lambda);
        mode.floquetExponent = log(lambda)/TPrimary;
        mode.unitCircleDefect = abs(abs(lambda)-1);
        mode.radialLogDefect = abs(log(abs(lambda)));
        mode.conjugacyDefect = conjugacyDefect(k);
        mode.eigenpairResidualRel = residualRel;
        mode.eigenvector = w;
        mode.basis6 = Q6;
        mode.planarProjection = norm(w([1 2 4 5]));
        mode.verticalProjection = norm(w([3 6]));
        mode.positionRankRatioAtSection = positionRankRatioAtSection;
        mode.continuationSegmentId = NaN;
        mode.segmentBranchId = NaN;
        mode.branchId = NaN;

        centerModes(end+1) = mode; %#ok<AGROW>
    end

    nontrivialModulus = modulus(nontrivialMask);
    nontrivialPhase = phase(nontrivialMask);
    nontrivialRadialLog = log(max(nontrivialModulus,realmin));
    nUnstableMultipliers = nnz(nontrivialRadialLog > opts.stabilityRadialLogTol);
    nStableMultipliers = nnz(nontrivialRadialLog < -opts.stabilityRadialLogTol);
    nUnitMultipliers = nnz(abs(nontrivialRadialLog) <= opts.stabilityRadialLogTol);
    spectralRadius = max(nontrivialModulus);

    % Permutation-invariant spectrum descriptor used only for family sampling.
    % The three sorted blocks retain radial and angular information without
    % relying on eig() column ordering.
    logModFeature = sort(asinh(abs(log(max(nontrivialModulus,realmin)))),'ascend');
    cosPhaseFeature = sort(cos(nontrivialPhase),'ascend');
    absSinPhaseFeature = sort(abs(sin(nontrivialPhase)),'ascend');
    featureVector = [ ...
        logModFeature(:).', ...
        cosPhaseFeature(:).', ...
        absSinPhaseFeature(:).'];

    spectrum = struct();
    spectrum.lambdaAll = lambdaAll;
    spectrum.modulusAll = modulus;
    spectrum.phaseAll = phase;
    spectrum.trivialMultiplierIndices = trivialIdx;
    spectrum.trivialMultiplierDistance = trivialDistances;
    spectrum.trivialPairResolved = trivialPairResolved;
    spectrum.trivialFlowAlignment = trivialFlowAlignment;
    spectrum.nontrivialMultiplierIndices = find(nontrivialMask);
    spectrum.centerModes = centerModes;
    spectrum.nCenterModes = numel(centerModes);
    spectrum.nUnstableMultipliers = nUnstableMultipliers;
    spectrum.nStableMultipliers = nStableMultipliers;
    spectrum.nUnitMultipliers = nUnitMultipliers;
    spectrum.complexQuartetDetected = complexQuartetDetected;
    spectrum.spectralRadius = spectralRadius;
    spectrum.featureVector = featureVector;
    spectrum.stabilitySignature = sprintf( ...
        'C%d_U%d_S%d_Q%d',numel(centerModes), ...
        nUnstableMultipliers,nStableMultipliers,complexQuartetDetected);
    spectrum.options = opts;
end


function f = local_cr3bp_vector_field(x,mu)
    X = x(1);
    Y = x(2);
    Z = x(3);
    VX = x(4);
    VY = x(5);
    VZ = x(6);

    r1 = sqrt((X+mu)^2+Y^2+Z^2);
    r2 = sqrt((X-1+mu)^2+Y^2+Z^2);
    if ~(isfinite(r1) && isfinite(r2) && r1 > 0 && r2 > 0)
        error('FINESST:CenterModes:PrimaryCollision', ...
            'The state lies at a CR3BP primary singularity.');
    end

    AX = 2*VY+X-(1-mu)*(X+mu)/r1^3-mu*(X-1+mu)/r2^3;
    AY = -2*VX+Y-(1-mu)*Y/r1^3-mu*Y/r2^3;
    AZ = -(1-mu)*Z/r1^3-mu*Z/r2^3;
    f = [VX;VY;VZ;AX;AY;AZ];
end
