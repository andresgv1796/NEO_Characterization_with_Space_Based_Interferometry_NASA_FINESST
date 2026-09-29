function fg = transport_floquet_center_geometry_CR3BP( ...
    tspan, xPO0, mu, w0, nSamples, ode_opts, nPhaseSamples)
%TRANSPORT_FLOQUET_CENTER_GEOMETRY_CR3BP
% Transport a selected complex Floquet center eigenvector around one CR3BP
% periodic orbit and extract scale-free formation-geometry predictors.
%
% Inputs
%   tspan         : [t0 tf], normally one parent-orbit period
%   xPO0          : 6x1 periodic-orbit initial state
%   mu            : CR3BP mass ratio
%   w0            : 6x1 complex Floquet eigenvector at t0
%   nSamples      : number of output points over the parent period
%   ode_opts      : ODE options; defaults to RelTol=AbsTol=1e-12
%   nPhaseSamples : common torus phases used for linear triangle-shape
%                   sensitivity; default 61 over [0,2pi/3)
%
% Theory
% ------
% The transported Floquet perturbation is
%
%       w(t) = Phi(t,t0) w0.
%
% Let the position part be
%
%       w_r(t) = a(t) + i b(t).
%
% Then a first-order invariant-curve section can be written as
%
%       delta r(t,theta)
%       = Re[w_r(t) exp(i theta)]
%       = a(t) cos(theta) - b(t) sin(theta).
%
% Define
%
%       A(t) = [a(t), -b(t)].
%
% The local configuration-space center ellipse is the image of the unit
% circle under A.  Therefore:
%
%   plane normal:
%       n_F(t) ~ a x b
%
%   ellipse singular values:
%       sigma_1 >= sigma_2
%
%   circularity:
%       chi = sigma_2/sigma_1
%
%   position RMS radius:
%       R_F = ||A||_F/sqrt(2)
%
% For three collectors separated by 120 degrees in center phase, the
% linearized RMS baseline obeys the exact identity
%
%       B_RMS = sqrt(3) R_F,
%
% independent of common collector phase.  Common phase can still change
% the individual side lengths when chi < 1.
%
% All principal predictors returned here are invariant to an arbitrary
% complex scaling w0 -> c exp(i alpha) w0, except the absolute R_F scale.
% The survey therefore uses only ratios / normalized quantities from R_F.
%
% Frames
% ------
% nSynodic is the plane normal in CR3BP rotating coordinates.
%
% nInertial uses the same convention as the FINESST fixed-target UV model:
%
%       r_I = Rz(+t) r_syn.
%
% This lets us distinguish intrinsic QPO-plane morphology from the inertial
% normal evolution relevant to fixed celestial targets.

    if nargin < 5 || isempty(nSamples)
        nSamples = 401;
    end

    if nargin < 6 || isempty(ode_opts)
        ode_opts = odeset( ...
            'RelTol',1e-12, ...
            'AbsTol',1e-12);
    end

    if nargin < 7 || isempty(nPhaseSamples)
        nPhaseSamples = 61;
    end

    assert(numel(tspan)==2,'tspan must contain [t0 tf].');
    assert(numel(xPO0)==6,'xPO0 must be 6x1.');
    assert(numel(w0)==6,'w0 must be 6x1 complex Floquet eigenvector.');
    assert(nSamples >= 2,'nSamples must be at least 2.');
    assert(nPhaseSamples >= 3,'nPhaseSamples must be at least 3.');

    xPO0 = xPO0(:);
    w0 = w0(:);

    tEval = linspace(tspan(1),tspan(2),nSamples);

    Phi0 = eye(6);

    Y0 = [ ...
        Phi0(:); ...
        xPO0];

    sol = ode89( ...
        @(t,Y) STM_vec(t,Y,mu), ...
        tEval, ...
        Y0, ...
        ode_opts);

    t = sol.x(:);
    Y = sol.y;

    nT = numel(t);

    xPO = Y(37:42,:).';

    wTrack = complex(zeros(nT,6));

    nSynodic = nan(nT,3);
    nInertial = nan(nT,3);

    stateVerticalFraction = nan(nT,1);
    positionVerticalFraction = nan(nT,1);
    velocityVerticalFraction = nan(nT,1);

    sigma1 = nan(nT,1);
    sigma2 = nan(nT,1);

    circularity = nan(nT,1);
    ellipseAspect = nan(nT,1);

    positionRMSRadius = nan(nT,1);
    linearBaselineRMS = nan(nT,1);

    phaseAspectMedian = nan(nT,1);
    phaseAspectMax = nan(nT,1);

    baselineIdentityRelError = nan(nT,1);

    phaseGrid = linspace(0,2*pi/3,nPhaseSamples+1);
    phaseGrid(end) = [];

    for k = 1:nT

        Phi = reshape(Y(1:36,k),[6,6]);

        wk = Phi*w0;

        wTrack(k,:) = wk.';

        wr = wk(1:3);
        wv = wk(4:6);

        a = real(wr);
        b = imag(wr);

        A = [a,-b];

        s = svd(A,'econ');

        if numel(s) >= 2

            sigma1(k) = s(1);
            sigma2(k) = s(2);

            if sigma1(k) > 0

                circularity(k) = ...
                    sigma2(k)/sigma1(k);

            end

            if sigma2(k) > 0

                ellipseAspect(k) = ...
                    sigma1(k)/sigma2(k);

            end

        end


        %% Plane normal

        nRaw = cross(a,b);

        nNorm = norm(nRaw);

        if nNorm > 1e-14

            nSynodic(k,:) = ...
                (nRaw/nNorm).';

        end


        %% Vertical fractions

        denomState = ...
            sum(abs(wk).^2);

        if denomState > 0

            stateVerticalFraction(k) = ...
                (abs(wk(3))^2 + abs(wk(6))^2)/ ...
                denomState;

        end


        denomPos = ...
            sum(abs(wr).^2);

        if denomPos > 0

            positionVerticalFraction(k) = ...
                abs(wr(3))^2/ ...
                denomPos;

        end


        denomVel = ...
            sum(abs(wv).^2);

        if denomVel > 0

            velocityVerticalFraction(k) = ...
                abs(wv(3))^2/ ...
                denomVel;

        end


        %% Linear ellipse scale

        positionRMSRadius(k) = ...
            norm(A,'fro')/sqrt(2);

        linearBaselineRMS(k) = ...
            sqrt(3)*positionRMSRadius(k);


        %% Common-phase sensitivity of the three-collector triangle
        %
        % The exact RMS baseline is common-phase invariant.  Individual
        % side lengths can vary with phase when the ellipse is non-circular.

        aspectByPhase = nan(nPhaseSamples,1);

        maxBrmsRelErrorThisTime = 0;

        for jp = 1:nPhaseSamples

            phi = phaseGrid(jp);

            theta = phi + ...
                2*pi*(0:2)/3;

            Q = A*[ ...
                cos(theta); ...
                sin(theta)];

            b12 = norm(Q(:,2)-Q(:,1));
            b13 = norm(Q(:,3)-Q(:,1));
            b23 = norm(Q(:,3)-Q(:,2));

            side = [b12,b13,b23];

            bMin = min(side);
            bMax = max(side);

            if bMin > 0

                aspectByPhase(jp) = ...
                    bMax/bMin;

            end


            bRMSDirect = ...
                sqrt(mean(side.^2));

            denom = ...
                max(linearBaselineRMS(k),eps);

            maxBrmsRelErrorThisTime = ...
                max( ...
                    maxBrmsRelErrorThisTime, ...
                    abs(bRMSDirect-linearBaselineRMS(k))/denom);

        end

        phaseAspectMedian(k) = ...
            median(aspectByPhase,'omitnan');

        phaseAspectMax(k) = ...
            max(aspectByPhase,[],'omitnan');

        baselineIdentityRelError(k) = ...
            maxBrmsRelErrorThisTime;

    end


    %% Enforce normal-sign continuity before transforming frames

    nSynodic = ...
        enforce_plane_normal_continuity(nSynodic);


    %% Synodic -> inertial normal

    for k = 1:nT

        if all(isfinite(nSynodic(k,:)))

            c = cos(t(k));
            s = sin(t(k));

            Rz = [ ...
                 c, -s, 0; ...
                 s,  c, 0; ...
                 0,  0, 1];

            nInertial(k,:) = ...
                (Rz*nSynodic(k,:).').';

        end

    end

    nInertial = ...
        enforce_plane_normal_continuity(nInertial);


    %% Sign-invariant normal-track metrics

    synMetrics = ...
        unoriented_normal_track_metrics(nSynodic);

    inertialMetrics = ...
        unoriented_normal_track_metrics(nInertial);


    %% Floquet closure checks

    validNormalSyn = ...
        all(isfinite(nSynodic),2);

    if nnz(validNormalSyn) >= 2

        iFirst = find(validNormalSyn,1,'first');
        iLast = find(validNormalSyn,1,'last');

        normalClosure_deg = ...
            acosd( ...
                clamp_unit( ...
                    abs(dot( ...
                        nSynodic(iFirst,:), ...
                        nSynodic(iLast,:)))));

    else

        normalClosure_deg = ...
            NaN;

    end


    lambdaEstimated = ...
        dot(w0,wTrack(end,:).')/ ...
        dot(w0,w0);

    floquetClosureResidual = ...
        norm( ...
            wTrack(end,:).' - ...
            lambdaEstimated*w0)/ ...
        max(norm(w0),eps);


    %% Scale-free breathing

    validRadius = ...
        positionRMSRadius > 0 & ...
        isfinite(positionRMSRadius);

    if any(validRadius)

        posRadiusMin = min(positionRMSRadius(validRadius));
        posRadiusMax = max(positionRMSRadius(validRadius));

        positionBreathingRatio = ...
            posRadiusMax/posRadiusMin;

    else

        positionBreathingRatio = NaN;

    end


    %% Package

    fg = struct();

    fg.t = t;
    fg.tau = (t-t(1))/(t(end)-t(1));

    fg.xPO = xPO;
    fg.wTrack = wTrack;

    fg.nSynodic = nSynodic;
    fg.nInertial = nInertial;

    fg.stateVerticalFraction = stateVerticalFraction;
    fg.positionVerticalFraction = positionVerticalFraction;
    fg.velocityVerticalFraction = velocityVerticalFraction;

    fg.sigma1 = sigma1;
    fg.sigma2 = sigma2;

    fg.circularity = circularity;
    fg.ellipseAspect = ellipseAspect;

    fg.positionRMSRadius = positionRMSRadius;
    fg.linearBaselineRMS = linearBaselineRMS;

    fg.phaseAspectMedian = phaseAspectMedian;
    fg.phaseAspectMax = phaseAspectMax;

    fg.baselineIdentityRelError = baselineIdentityRelError;

    fg.lambdaEstimated = lambdaEstimated;
    fg.floquetClosureResidual = floquetClosureResidual;
    fg.normalClosure_deg = normalClosure_deg;

    fg.metrics = struct();

    fg.metrics.synodic = synMetrics;
    fg.metrics.inertial = inertialMetrics;

    fg.metrics.meanStateVerticalFraction = ...
        mean(stateVerticalFraction,'omitnan');

    fg.metrics.meanPositionVerticalFraction = ...
        mean(positionVerticalFraction,'omitnan');

    fg.metrics.maxPositionVerticalFraction = ...
        max(positionVerticalFraction,[],'omitnan');

    fg.metrics.meanVelocityVerticalFraction = ...
        mean(velocityVerticalFraction,'omitnan');

    fg.metrics.medianCircularity = ...
        median(circularity,'omitnan');

    fg.metrics.minCircularity = ...
        min(circularity,[],'omitnan');

    fg.metrics.medianEllipseAspect = ...
        median(ellipseAspect,'omitnan');

    fg.metrics.maxEllipseAspect = ...
        max(ellipseAspect,[],'omitnan');

    fg.metrics.positionBreathingRatio = ...
        positionBreathingRatio;

    fg.metrics.medianPhaseAspect = ...
        median(phaseAspectMedian,'omitnan');

    fg.metrics.maxPhaseAspect = ...
        max(phaseAspectMax,[],'omitnan');

    fg.metrics.maxBaselineIdentityRelError = ...
        max(baselineIdentityRelError,[],'omitnan');

end


function nTrack = enforce_plane_normal_continuity(nTrack)
%ENFORCE_PLANE_NORMAL_CONTINUITY
% Flip n -> -n when needed so the displayed track is locally continuous.
% This does not alter the physical plane, for which n and -n are equivalent.

    valid = all(isfinite(nTrack),2);

    idx = find(valid);

    if numel(idx) < 2
        return
    end

    for j = 2:numel(idx)

        kPrev = idx(j-1);
        k = idx(j);

        if dot(nTrack(kPrev,:),nTrack(k,:)) < 0

            nTrack(k,:) = ...
                -nTrack(k,:);

        end

    end

end


function M = unoriented_normal_track_metrics(nTrack)
%UNORIENTED_NORMAL_TRACK_METRICS
% Sign-invariant summary of a plane-normal track.
%
% A representative mean plane axis is obtained from the dominant eigenvector
% of sum(n n^T), which is appropriate when n and -n represent the same plane.

    valid = all(isfinite(nTrack),2);

    N = nTrack(valid,:);

    M = struct();

    if isempty(N)

        M.meanAxis = [NaN NaN NaN];
        M.RMSExcursion_deg = NaN;
        M.maxExcursion_deg = NaN;
        M.pathLength_deg = NaN;
        return

    end

    Q = N.'*N;

    [V,D] = eig(Q);

    [~,iMax] = max(diag(D));

    nMean = V(:,iMax);
    nMean = nMean/norm(nMean);

    d = abs(N*nMean);

    d = arrayfun(@clamp_unit,d);

    angleFromMean = acosd(d);

    pathLength_deg = 0;

    for k = 2:size(N,1)

        c = clamp_unit(abs(dot(N(k-1,:),N(k,:))));

        pathLength_deg = ...
            pathLength_deg + ...
            acosd(c);

    end

    M.meanAxis = nMean(:).';
    M.RMSExcursion_deg = sqrt(mean(angleFromMean.^2));
    M.maxExcursion_deg = max(angleFromMean);
    M.pathLength_deg = pathLength_deg;

end


function x = clamp_unit(x)
%CLAMP_UNIT

    x = min(1,max(-1,x));

end
