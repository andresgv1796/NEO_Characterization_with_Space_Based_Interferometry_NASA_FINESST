function out = transport_floquet_center_normal_CR3BP( ...
    xPO0,TPrimary,mu,mode,nPrimaryPeriods,nSavePerPrimaryPeriod, ...
    synodicPhase0_rad,odeOpts)
%TRANSPORT_FLOQUET_CENTER_NORMAL_CR3BP Transport a periodic Floquet plane normal.
%
% out = transport_floquet_center_normal_CR3BP(...)
%
% Computes a center Floquet eigendirection from the one-period STM of a
% CR3BP periodic orbit, constructs the periodic position-space Floquet plane,
% and returns its sign-continuous normal over an integer number of primary
% periods.
%
% IMPORTANT SAMPLING CONVENTION
%   The time-dependent plane normal is evaluated at midpoint quadrature nodes
%
%       t_j = (j-1/2) TPrimary/N,   j = 1,...,N,
%
%   rather than at periodic-orbit section endpoints. At an exact symmetry
%   section, the POSITION projection of an otherwise valid complex Floquet
%   mode can collapse instantaneously to rank one because of phase symmetry.
%   The plane normal is then genuinely undefined at that isolated instant.
%   Such measure-zero instants do not contribute to the time integral used
%   by the sky-suitability study and are therefore not sampled.
%
%   This is NOT a fallback: if the position-space Floquet plane is rank
%   deficient at any midpoint quadrature node, the function throws an error.
%
% Inputs
%   xPO0                    6x1 periodic-orbit initial state
%   TPrimary                primary-orbit period [TU]
%   mu                      CR3BP mass parameter
%   mode                    'planar' or 'vertical'
%   nPrimaryPeriods         positive integer
%   nSavePerPrimaryPeriod   positive integer >= 20; midpoint samples/period
%   synodicPhase0_rad       inertial rotation phase at t=0 [rad]
%   odeOpts                 explicit ODE options passed to ode89
%
% Outputs (selected fields)
%   out.t                    1xNt midpoint time history [TU]
%   out.nSynodic             3xNt sign-continuous Floquet plane normal
%   out.nInertial            3xNt inertial Floquet plane normal
%   out.chi                  1xNt periodic ellipse circularity sigma2/sigma1
%   out.tOnePeriod           1xN midpoint nodes in one primary period
%   out.nSynodicOnePeriod    3xN one-period normal line representative
%   out.chiOnePeriod         1xN one-period circularity
%   out.lambda               selected center multiplier
%   out.rho                  angle(out.lambda)
%   out.modeProjection       planar/vertical normalized projection scores
%   out.floquetClosureRel    relative full-state Floquet-factor closure error
%   out.stateClosureNorm     periodic-orbit state closure norm over one period
%
% Design choices
%   * State + STM are integrated for ONE primary period only.
%   * The periodic Floquet factor p(t)=Phi(t)w exp(-i*rho*t/T) is formed.
%   * Periodicity is validated in the full 6-D Floquet factor, not by forcing
%     a configuration-space normal to exist at a symmetry section.
%   * Position-plane normals are evaluated only where the 3x2 projected basis
%     is genuinely rank two.
%   * No fallback normal, interpolation, rank repair, or alternate mode is
%     substituted.
%   * The validated periodic normal line is repeated for nPrimaryPeriods.
%
% Existing dependencies
%   STM_vec, A_CR3BP, CR3BP_mu, ode89

    narginchk(8,8);

    validateattributes(xPO0,{'numeric'},{'real','finite','vector','numel',6}, ...
        mfilename,'xPO0',1);
    xPO0 = xPO0(:);

    validateattributes(TPrimary,{'numeric'},{'real','finite','scalar','positive'}, ...
        mfilename,'TPrimary',2);
    validateattributes(mu,{'numeric'},{'real','finite','scalar','>',0,'<',0.5}, ...
        mfilename,'mu',3);
    validateattributes(nPrimaryPeriods,{'numeric'}, ...
        {'real','finite','scalar','integer','positive'},mfilename,'nPrimaryPeriods',5);
    validateattributes(nSavePerPrimaryPeriod,{'numeric'}, ...
        {'real','finite','scalar','integer','>=',20},mfilename,'nSavePerPrimaryPeriod',6);
    validateattributes(synodicPhase0_rad,{'numeric'}, ...
        {'real','finite','scalar'},mfilename,'synodicPhase0_rad',7);

    mode = lower(string(mode));
    if ~isscalar(mode) || ~(mode == "planar" || mode == "vertical")
        error('FINESST:FloquetNormal:InvalidMode', ...
            'mode must be exactly ''planar'' or ''vertical''.');
    end

    if exist('ode89','file') ~= 2
        error('FINESST:FloquetNormal:MissingOde89', ...
            'ode89 is required. No alternate ODE solver is substituted.');
    end

    if ~isstruct(odeOpts)
        error('FINESST:FloquetNormal:InvalidOdeOptions', ...
            'odeOpts must be an explicit ODE options structure.');
    end

    unitCircleTol = 1e-6;
    imagTol = 1e-10;
    planeRankRelTol = 1e-10;
    continuityAbsDotTol = 1e-6;
    floquetClosureRelTol = 1e-8;

    % One-period integration grid: endpoints for monodromy/closure and
    % midpoint nodes for the configuration-space plane geometry.
    N = nSavePerPrimaryPeriod;
    dt = TPrimary/N;
    tOne = ((0:N-1) + 0.5)*dt;
    tIntegrate = [0, tOne, TPrimary];

    Y0 = zeros(42,1);
    Phi0 = eye(6);
    Y0(1:36) = Phi0(:);
    Y0(37:42) = xPO0;

    [tSol,YSol] = ode89(@(t,Y) STM_vec(t,Y,mu),tIntegrate,Y0,odeOpts);

    if numel(tSol) ~= numel(tIntegrate) || ...
       size(YSol,1) ~= numel(tIntegrate) || size(YSol,2) ~= 42
        error('FINESST:FloquetNormal:UnexpectedIntegratorOutput', ...
            'ode89 did not return the requested one-period STM grid.');
    end

    if any(~isfinite(YSol(:)))
        error('FINESST:FloquetNormal:NonfiniteSTM', ...
            'Nonfinite values occurred during one-period STM propagation.');
    end

    PhiT = reshape(YSol(end,1:36),6,6);
    xT = YSol(end,37:42).';
    stateClosureNorm = norm(xT-xPO0);

    [V,D] = eig(PhiT);
    lambdaAll = diag(D);

    candidateIdx = find( ...
        abs(abs(lambdaAll)-1) <= unitCircleTol & ...
        imag(lambdaAll) > imagTol);

    if isempty(candidateIdx)
        error('FINESST:FloquetNormal:NoCenterMode', ...
            ['No positive-imaginary complex unit-circle Floquet multiplier ' ...
             'was found within | |lambda|-1 | <= %.1e.'],unitCircleTol);
    end

    planarProjection = zeros(numel(candidateIdx),1);
    verticalProjection = zeros(numel(candidateIdx),1);

    for k = 1:numel(candidateIdx)
        wCandidate = V(:,candidateIdx(k));
        nw = norm(wCandidate);
        if ~(isfinite(nw) && nw > 0)
            error('FINESST:FloquetNormal:InvalidEigenvector', ...
                'A center-multiplier eigenvector has invalid norm.');
        end
        wCandidate = wCandidate/nw;
        planarProjection(k) = norm(wCandidate([1 2 4 5]));
        verticalProjection(k) = norm(wCandidate([3 6]));
    end

    if mode == "planar"
        requestedProjection = planarProjection;
    else
        requestedProjection = verticalProjection;
    end

    [~,kBest] = max(requestedProjection);
    iBest = candidateIdx(kBest);

    lambda = lambdaAll(iBest);
    w = V(:,iBest);
    w = w/norm(w);
    rho = angle(lambda);
    beta = rho/TPrimary;

    modeProjection = struct();
    modeProjection.planar = norm(w([1 2 4 5]));
    modeProjection.vertical = norm(w([3 6]));

    if mode == "vertical" && modeProjection.vertical < modeProjection.planar
        warning('FINESST:FloquetNormal:MixedMode', ...
            ['Requested vertical center branch is planar-dominant in the ' ...
             'selected eigenvector (planar %.3f, vertical %.3f). The ' ...
             'candidate is retained and explicitly reported; no alternate ' ...
             'mode is substituted.'], ...
            modeProjection.planar,modeProjection.vertical);
    elseif mode == "planar" && modeProjection.planar < modeProjection.vertical
        warning('FINESST:FloquetNormal:MixedMode', ...
            ['Requested planar center branch is vertical-dominant in the ' ...
             'selected eigenvector (planar %.3f, vertical %.3f). The ' ...
             'candidate is retained and explicitly reported; no alternate ' ...
             'mode is substituted.'], ...
            modeProjection.planar,modeProjection.vertical);
    end

    % Full 6-D Floquet closure remains well defined at exact symmetry sections.
    p0 = w;
    pT = PhiT*w*exp(-1i*beta*TPrimary);
    floquetClosureRel = norm(pT-p0)/norm(p0);

    if ~(isfinite(floquetClosureRel) && floquetClosureRel <= floquetClosureRelTol)
        error('FINESST:FloquetNormal:PoorFloquetClosure', ...
            ['Periodic full-state Floquet factor failed one-period closure: ' ...
             'relative error %.6e (allowed %.6e).'], ...
            floquetClosureRel,floquetClosureRelTol);
    end

    % Position-space Floquet plane on midpoint quadrature nodes.
    nOne = nan(3,N);
    chiOne = nan(1,N);

    for k = 1:N
        Phi = reshape(YSol(k+1,1:36),6,6);
        p = Phi*w*exp(-1i*beta*tOne(k));
        Apos = [real(p(1:3)), imag(p(1:3))];

        if any(~isfinite(Apos(:)))
            error('FINESST:FloquetNormal:NonfinitePlane', ...
                'Nonfinite Floquet position-plane basis at midpoint sample %d.',k);
        end

        s = svd(Apos,'econ');
        if numel(s) ~= 2 || ~(s(1) > 0)
            error('FINESST:FloquetNormal:PositionPlaneRankDeficient', ...
                'Floquet position plane is rank-deficient at midpoint sample %d.',k);
        end

        chiOne(k) = s(2)/s(1);
        if s(2) <= planeRankRelTol*s(1)
            error('FINESST:FloquetNormal:PositionPlaneRankDeficient', ...
                ['Floquet position plane is numerically rank-deficient at ' ...
                 'midpoint sample %d: sigma2/sigma1 = %.3e.'],k,chiOne(k));
        end

        c = cross(Apos(:,1),Apos(:,2));
        nc = norm(c);
        if ~(isfinite(nc) && nc > 0)
            error('FINESST:FloquetNormal:PositionPlaneRankDeficient', ...
                'Floquet plane normal is invalid at midpoint sample %d.',k);
        end

        nOne(:,k) = c/nc;
    end

    % Deterministic sign continuity within one period.
    for k = 2:N
        d = dot(nOne(:,k-1),nOne(:,k));
        if abs(d) <= continuityAbsDotTol
            error('FINESST:FloquetNormal:NormalContinuityFailure', ...
                ['Consecutive Floquet normal lines are nearly orthogonal at ' ...
                 'midpoint samples %d-%d (|dot|=%.3e). Increase temporal ' ...
                 'resolution or inspect the selected center geometry.'], ...
                k-1,k,abs(d));
        end
        if d < 0
            nOne(:,k) = -nOne(:,k);
        end
    end

    % Neighboring line-angle across the periodic seam; separated by one dt.
    seamAbsDot = abs(dot(nOne(:,end),nOne(:,1)));
    if seamAbsDot <= continuityAbsDotTol
        error('FINESST:FloquetNormal:NormalContinuityFailure', ...
            ['Floquet normal line is nearly orthogonal across the periodic ' ...
             'seam (|dot|=%.3e). Increase temporal resolution or inspect ' ...
             'the selected center geometry.'],seamAbsDot);
    end
    periodicSeamStepDeg = acosd(min(1,max(0,seamAbsDot)));

    % Repeat the periodic normal line over exactly nPrimaryPeriods.
    nTotal = nPrimaryPeriods*N;
    tFull = ((0:nTotal-1) + 0.5)*dt;
    periodicIndex = mod(0:nTotal-1,N) + 1;
    nSynodic = nOne(:,periodicIndex);
    chi = chiOne(periodicIndex);

    % Make the displayed representative sign-continuous across all periods.
    for k = 2:nTotal
        d = dot(nSynodic(:,k-1),nSynodic(:,k));
        if abs(d) <= continuityAbsDotTol
            error('FINESST:FloquetNormal:NormalContinuityFailure', ...
                ['Repeated Floquet normal lines are nearly orthogonal at ' ...
                 'full-history samples %d-%d (|dot|=%.3e).'],k-1,k,abs(d));
        end
        if d < 0
            nSynodic(:,k) = -nSynodic(:,k);
        end
    end

    % Synodic -> inertial. Standard CR3BP nondimensional angular rate = 1.
    nInertial = zeros(size(nSynodic));
    for k = 1:nTotal
        th = synodicPhase0_rad + tFull(k);
        cth = cos(th);
        sth = sin(th);
        R_IS = [cth -sth 0; sth cth 0; 0 0 1];
        nInertial(:,k) = R_IS*nSynodic(:,k);
    end

    normError = max(abs(vecnorm(nInertial,2,1)-1));
    if normError > 1e-10
        error('FINESST:FloquetNormal:UnitNormalFailure', ...
            'Inertial normal normalization error reached %.3e.',normError);
    end

    out = struct();
    out.t = tFull;
    out.TPrimary = TPrimary;
    out.nPrimaryPeriods = nPrimaryPeriods;
    out.nSavePerPrimaryPeriod = nSavePerPrimaryPeriod;
    out.sampleConvention = 'uniform midpoint quadrature';
    out.dt = dt;
    out.synodicPhase0_rad = synodicPhase0_rad;

    out.nSynodic = nSynodic;
    out.nInertial = nInertial;
    out.chi = chi;

    out.tOnePeriod = tOne;
    out.nSynodicOnePeriod = nOne;
    out.chiOnePeriod = chiOne;

    out.lambda = lambda;
    out.rho = rho;
    out.beta = beta;
    out.eigenvector = w;
    out.requestedMode = char(mode);
    out.modeProjection = modeProjection;

    out.floquetClosureRel = floquetClosureRel;
    out.stateClosureNorm = stateClosureNorm;
    out.periodicSeamStepDeg = periodicSeamStepDeg;

    out.centerCandidateEigenvalues = lambdaAll(candidateIdx);
    out.centerCandidatePlanarProjection = planarProjection;
    out.centerCandidateVerticalProjection = verticalProjection;
end
