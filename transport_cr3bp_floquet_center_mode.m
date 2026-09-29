function out = transport_cr3bp_floquet_center_mode( ...
    xPO0,TPrimary,mu,modeReference,nPrimaryPeriods, ...
    nSavePerPrimaryPeriod,synodicPhase0_rad,odeOpts,opts)
%TRANSPORT_CR3BP_FLOQUET_CENTER_MODE Transport one explicitly identified mode.
%
% out = transport_cr3bp_floquet_center_mode(...)
%
% modeReference must come from enumerate_cr3bp_center_modes, normally through
% analyze_cr3bp_family_floquet_structure. The mode is reidentified against
% the dense-grid monodromy by Floquet angle and real-subspace overlap. No
% planar/vertical request and no alternate-mode substitution are used.
%
% The position-plane normal is evaluated at uniform midpoint quadrature nodes.
% If the selected state-space center mode projects to a rank-deficient spatial
% plane at any node, the mode is explicitly inadmissible for this normal-based
% surrogate.

    narginchk(9,9);

    validateattributes(xPO0,{'numeric'}, ...
        {'real','finite','vector','numel',6},mfilename,'xPO0',1);
    xPO0 = xPO0(:);
    validateattributes(TPrimary,{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,'TPrimary',2);
    validateattributes(mu,{'numeric'}, ...
        {'real','finite','scalar','>',0,'<',0.5},mfilename,'mu',3);
    validateattributes(nPrimaryPeriods,{'numeric'}, ...
        {'real','finite','scalar','integer','positive'}, ...
        mfilename,'nPrimaryPeriods',5);
    validateattributes(nSavePerPrimaryPeriod,{'numeric'}, ...
        {'real','finite','scalar','integer','>=',20}, ...
        mfilename,'nSavePerPrimaryPeriod',6);
    validateattributes(synodicPhase0_rad,{'numeric'}, ...
        {'real','finite','scalar'},mfilename,'synodicPhase0_rad',7);

    requiredReference = {'lambda','rho','basis6','branchId'};
    if ~isstruct(modeReference) || ~all(isfield(modeReference,requiredReference))
        error('FINESST:ModeTransport:InvalidModeReference', ...
            'modeReference must explicitly contain: %s.', ...
            strjoin(requiredReference,', '));
    end
    if ~isequal(size(modeReference.basis6),[6 2]) || ...
            any(~isfinite(modeReference.basis6(:)))
        error('FINESST:ModeTransport:InvalidModeReferenceBasis', ...
            'modeReference.basis6 must be a finite 6x2 real basis.');
    end
    if exist('ode89','file') ~= 2
        error('FINESST:ModeTransport:MissingOde89', ...
            'ode89 is required; no alternate solver is substituted.');
    end
    if ~isstruct(odeOpts)
        error('FINESST:ModeTransport:InvalidOdeOptions', ...
            'odeOpts must be an explicit ODE options structure.');
    end

    requiredOpts = { ...
        'nTrivialMultipliers','trivialPairMaxDistance', ...
        'trivialFlowAlignmentMin','trivialMultiplierTol','imagTol', ...
        'centerRadialLogTol','centerConjugacyTol','stabilityRadialLogTol', ...
        'complexQuartetAngleTol','complexQuartetRadialTol', ...
        'eigenpairResidualRelTol','planeRankRelTol', ...
        'continuityAbsDotTol','floquetClosureRelTol', ...
        'modeMatchMaxCost','modeMatchUniquenessTol','normalNormTol'};
    if ~isstruct(opts) || ~all(isfield(opts,requiredOpts))
        error('FINESST:ModeTransport:InvalidOptions', ...
            'opts must explicitly contain: %s.',strjoin(requiredOpts,', '));
    end

    N = nSavePerPrimaryPeriod;
    dt = TPrimary/N;
    tOne = ((0:N-1)+0.5)*dt;
    tIntegrate = [0,tOne,TPrimary];

    Y0 = zeros(42,1);
    Y0(1:36) = reshape(eye(6),36,1);
    Y0(37:42) = xPO0;

    [tSol,YSol] = ode89(@(t,Y) STM_vec(t,Y,mu),tIntegrate,Y0,odeOpts);
    if numel(tSol) ~= numel(tIntegrate) || ...
            size(YSol,1) ~= numel(tIntegrate) || size(YSol,2) ~= 42
        error('FINESST:ModeTransport:UnexpectedIntegratorOutput', ...
            'ode89 did not return the requested one-period STM grid.');
    end
    if any(~isfinite(YSol(:)))
        error('FINESST:ModeTransport:NonfinitePropagation', ...
            'Nonfinite state/STM values occurred during mode transport.');
    end

    PhiT = reshape(YSol(end,1:36),6,6);
    xT = YSol(end,37:42).';
    stateClosureNorm = norm(xT-xPO0);

    currentSpectrum = enumerate_cr3bp_center_modes( ...
        PhiT,TPrimary,xPO0,mu,opts);
    currentModes = currentSpectrum.centerModes;
    if isempty(currentModes)
        error('FINESST:ModeTransport:ReferenceModeMissing', ...
            'No nontrivial center mode exists in the dense-grid monodromy.');
    end

    matchCost = nan(numel(currentModes),1);
    matchOverlap = nan(numel(currentModes),1);
    matchRhoDistance = nan(numel(currentModes),1);

    for k = 1:numel(currentModes)
        overlap = norm(modeReference.basis6.'*currentModes(k).basis6,'fro')/sqrt(2);
        overlap = min(1,max(0,overlap));
        rhoDistance = abs(modeReference.rho-currentModes(k).rho)/pi;
        matchOverlap(k) = overlap;
        matchRhoDistance(k) = rhoDistance;
        matchCost(k) = rhoDistance + 1-overlap;
    end

    [sortedCost,sortedIndex] = sort(matchCost,'ascend');
    if ~(isfinite(sortedCost(1)) && sortedCost(1) <= opts.modeMatchMaxCost)
        error('FINESST:ModeTransport:ReferenceModeMismatch', ...
            ['The closest dense-grid center mode has match cost %.6e ' ...
             '(allowed %.6e).'],sortedCost(1),opts.modeMatchMaxCost);
    end
    if numel(sortedCost) > 1 && ...
            sortedCost(2)-sortedCost(1) <= opts.modeMatchUniquenessTol
        error('FINESST:ModeTransport:AmbiguousModeMatch', ...
            ['The two closest dense-grid modes have costs %.6e and %.6e; ' ...
             'the reference match is not unique within %.6e.'], ...
            sortedCost(1),sortedCost(2),opts.modeMatchUniquenessTol);
    end

    selectedMode = currentModes(sortedIndex(1));
    lambda = selectedMode.lambda;
    w = selectedMode.eigenvector;
    floquetExponent = log(lambda)/TPrimary;

    p0 = w;
    pT = PhiT*w*exp(-floquetExponent*TPrimary);
    floquetClosureRel = norm(pT-p0)/norm(p0);
    if ~(isfinite(floquetClosureRel) && ...
            floquetClosureRel <= opts.floquetClosureRelTol)
        error('FINESST:ModeTransport:PoorFloquetClosure', ...
            ['The full-state periodic Floquet factor has closure %.6e ' ...
             '(allowed %.6e).'], ...
            floquetClosureRel,opts.floquetClosureRelTol);
    end

    nOne = nan(3,N);
    chiOne = nan(1,N);

    for k = 1:N
        Phi = reshape(YSol(k+1,1:36),6,6);
        p = Phi*w*exp(-floquetExponent*tOne(k));
        Apos = [real(p(1:3)),imag(p(1:3))];

        if any(~isfinite(Apos(:)))
            error('FINESST:ModeTransport:NonfinitePositionPlane', ...
                'Nonfinite position-plane basis at midpoint sample %d.',k);
        end

        s = svd(Apos,'econ');
        if numel(s) ~= 2 || ~(s(1) > 0)
            error('FINESST:ModeTransport:PositionPlaneRankDeficient', ...
                'The position plane is rank deficient at midpoint sample %d.',k);
        end

        chiOne(k) = s(2)/s(1);
        if s(2) <= opts.planeRankRelTol*s(1)
            error('FINESST:ModeTransport:PositionPlaneRankDeficient', ...
                ['The position plane is numerically rank deficient at ' ...
                 'midpoint sample %d: sigma2/sigma1 = %.6e.'],k,chiOne(k));
        end

        c = cross(Apos(:,1),Apos(:,2));
        nc = norm(c);
        if ~(isfinite(nc) && nc > 0)
            error('FINESST:ModeTransport:PositionPlaneRankDeficient', ...
                'The position-plane normal is invalid at midpoint sample %d.',k);
        end
        nOne(:,k) = c/nc;
    end

    for k = 2:N
        d = dot(nOne(:,k-1),nOne(:,k));
        if abs(d) <= opts.continuityAbsDotTol
            error('FINESST:ModeTransport:NormalContinuityFailure', ...
                ['Adjacent normal lines are nearly orthogonal at samples ' ...
                 '%d-%d (|dot|=%.6e).'],k-1,k,abs(d));
        end
        if d < 0
            nOne(:,k) = -nOne(:,k);
        end
    end

    seamAbsDot = abs(dot(nOne(:,end),nOne(:,1)));
    if seamAbsDot <= opts.continuityAbsDotTol
        error('FINESST:ModeTransport:NormalContinuityFailure', ...
            'The normal line is nearly orthogonal across the periodic seam.');
    end
    periodicSeamStepDeg = acosd(min(1,max(0,seamAbsDot)));

    nTotal = nPrimaryPeriods*N;
    tFull = ((0:nTotal-1)+0.5)*dt;
    periodicIndex = mod(0:nTotal-1,N)+1;
    nSynodic = nOne(:,periodicIndex);
    chi = chiOne(periodicIndex);

    for k = 2:nTotal
        if dot(nSynodic(:,k-1),nSynodic(:,k)) < 0
            nSynodic(:,k) = -nSynodic(:,k);
        end
    end

    nInertial = zeros(size(nSynodic));
    for k = 1:nTotal
        theta = synodicPhase0_rad+tFull(k);
        R_IS = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
        nInertial(:,k) = R_IS*nSynodic(:,k);
    end

    normalNormError = max(abs(vecnorm(nInertial,2,1)-1));
    if normalNormError > opts.normalNormTol
        error('FINESST:ModeTransport:UnitNormalFailure', ...
            'Inertial normal norm error reached %.6e.',normalNormError);
    end

    out = struct();
    out.t = tFull;
    out.TPrimary = TPrimary;
    out.nPrimaryPeriods = nPrimaryPeriods;
    out.nSavePerPrimaryPeriod = N;
    out.dt = dt;
    out.sampleConvention = 'uniform midpoint quadrature';
    out.nSynodic = nSynodic;
    out.nInertial = nInertial;
    out.chi = chi;
    out.tOnePeriod = tOne;
    out.nSynodicOnePeriod = nOne;
    out.parentOrbitSynodic = [xPO0(1:3),YSol(2:end-1,37:39).',xT(1:3)];
    out.chiOnePeriod = chiOne;
    out.modeReference = modeReference;
    out.selectedMode = selectedMode;
    out.lambda = lambda;
    out.rho = selectedMode.rho;
    out.floquetExponent = floquetExponent;
    out.unitCircleDefect = selectedMode.unitCircleDefect;
    out.eigenpairResidualRel = selectedMode.eigenpairResidualRel;
    out.floquetClosureRel = floquetClosureRel;
    out.stateClosureNorm = stateClosureNorm;
    out.periodicSeamStepDeg = periodicSeamStepDeg;
    out.modeMatchCost = sortedCost(1);
    out.modeMatchOverlap = matchOverlap(sortedIndex(1));
    out.modeMatchRhoDistance = matchRhoDistance(sortedIndex(1));
end
