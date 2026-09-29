function out = evaluate_finite_qpo_vs_floquet( ...
    XCurve, TQPO, rhoQPO, po, sys, floquetGeom, ...
    nCollectors, phaseStartIndices, ...
    nPrimaryPeriods, nSavePerPrimaryPeriod, ...
    odeOpts, normalReferenceSynodic)
%EVALUATE_FINITE_QPO_VS_FLOQUET
% Compare one corrected finite-amplitude QPO torus with the linear Floquet
% plane predicted from its parent periodic orbit.
%
% For every physically distinct equal-phase 3-collector starting phase:
%
%   1. propagate the finite formation;
%   2. compare its plane normal against the periodic Floquet-plane normal;
%   3. measure one-period non-closure of the finite formation plane;
%   4. compare finite RMS baseline breathing with the linearly predicted
%      breathing after calibrating the linear scale using the actual
%      corrected torus-section RMS radius.
%
% The phase ensemble also provides a direct nonlinear phase-sensitivity
% metric.  In the linear Floquet limit that spread must vanish.

    N = size(XCurve,2);

    if isempty(phaseStartIndices)

        phaseStartIndices = ...
            unique_equal_phase_start_indices( ...
                N, ...
                nCollectors);

    end

    nPhases = numel(phaseStartIndices);

    TPrimary = po.T;

    xCombiner0 = po.x0(:);

    %% Corrected section size

    rCenter = mean(XCurve(1:3,:),2);

    dR = XCurve(1:3,:) - rCenter;

    rSamples = vecnorm(dR,2,1);

    sectionRadiusRMS_ND = ...
        sqrt(mean(rSamples.^2));

    sectionRadiusRMS_km = ...
        sectionRadiusRMS_ND*sys.Lstar_km;


    %% Linear baseline scale calibrated to the finite section RMS radius

    RF0 = floquetGeom.positionRMSRadius(1);

    if ~(isfinite(RF0) && RF0 > 0)
        error('Invalid Floquet position RMS radius at t=0.');
    end

    linearScale = ...
        sectionRadiusRMS_ND/RF0;


    %% Per-phase storage

    floquetErrorFirstRMS_deg = nan(nPhases,1);
    floquetErrorFirstMax_deg = nan(nPhases,1);

    floquetErrorAllRMS_deg = nan(nPhases,1);
    floquetErrorAllMax_deg = nan(nPhases,1);

    periodClosureRMS_deg = nan(nPhases,1);
    periodClosureMax_deg = nan(nPhases,1);

    baselineRelErrorRMS = nan(nPhases,1);
    baselineRelErrorMax = nan(nPhases,1);

    medianBaseline_km = nan(nPhases,1);

    % Instantaneous triangle-plane conditioning:
    %
    %              4 sqrt(3) A_triangle
    %   q_A = --------------------------------
    %         B12^2 + B13^2 + B23^2
    %
    % where A_triangle = 0.5*||B12 x B13||.
    %
    % This symmetric dimensionless area is independent of collector
    % labeling, equals 1 for an equilateral triangle, and tends to 0 as the
    % three collectors become collinear.  Thus small q_A means the
    % formation-plane normal is geometrically ill-conditioned.
    qAMinFirstByPhase = nan(nPhases,1);

    normalTracksFirst = cell(nPhases,1);

    referencePhaseTrack = struct();


    % The parent-orbit propagation grid is independent of common collector
    % phase.  Cache the exact-phase Floquet predictor after the first phase
    % case and reuse it for the remaining phase starts.
    floquetCache = struct();

    floquetCache.t = [];
    floquetCache.n = [];
    floquetCache.B_RMS_unit = [];


    for iph = 1:nPhases

        phaseStartIndex = ...
            phaseStartIndices(iph);


        jCollectors = ...
            qpo_equal_phase_indices( ...
                N, ...
                nCollectors, ...
                phaseStartIndex);


        XCollectors0 = ...
            XCurve(:,jCollectors);


        geom = ...
            compute_three_collector_qpo_geometry( ...
                sys.mu, ...
                xCombiner0, ...
                TPrimary, ...
                XCollectors0, ...
                nPrimaryPeriods, ...
                nSavePerPrimaryPeriod, ...
                odeOpts, ...
                normalReferenceSynodic);


        t = extract_time_track(geom);


        [B12,B13,B23] = ...
            extract_three_baselines_nd(geom);


        nQPO = ...
            normal_from_baselines( ...
                B12, ...
                B13);


        %% Floquet predictor at the exact finite-QPO sample phases
        %
        % Avoid interpolation of an already-normalized Floquet-normal track.
        % At very small formation size that interpolation can be comparable
        % to the finite-amplitude effect being measured.

        sameFloquetGrid = ...
            ~isempty(floquetCache.t) && ...
            numel(floquetCache.t) == numel(t) && ...
            max(abs(floquetCache.t(:)-t(:))) <= ...
                1e-12*max(1,TPrimary);


        if sameFloquetGrid

            nFloquet = ...
                floquetCache.n;

            B_RMS_Floquet_unit = ...
                floquetCache.B_RMS_unit;

        else

            [nFloquet,B_RMS_Floquet_unit] = ...
                direct_floquet_predictor_at_times( ...
                    po, ...
                    sys, ...
                    floquetGeom, ...
                    t, ...
                    odeOpts);

            floquetCache.t = ...
                t(:);

            floquetCache.n = ...
                nFloquet;

            floquetCache.B_RMS_unit = ...
                B_RMS_Floquet_unit;

        end


        angleAll = ...
            plane_angle_deg( ...
                nQPO, ...
                nFloquet);


        firstMask = ...
            t <= t(1)+TPrimary*(1+1e-11);


        angleFirst = ...
            angleAll(firstMask);


        % Symmetric triangle-area conditioning over the same first-period
        % interval.
        b12Now = vecnorm(B12,2,2);
        b13Now = vecnorm(B13,2,2);
        b23Now = vecnorm(B23,2,2);

        crossNow = cross(B12,B13,2);

        qA = ...
            2*sqrt(3)*vecnorm(crossNow,2,2) ./ ...
            max( ...
                b12Now.^2 + b13Now.^2 + b23Now.^2, ...
                realmin);

        qA = min(1,max(0,qA));

        qAMinFirstByPhase(iph) = ...
            min(qA(firstMask),[],'omitnan');


        floquetErrorFirstRMS_deg(iph) = ...
            sqrt(mean(angleFirst.^2,'omitnan'));

        floquetErrorFirstMax_deg(iph) = ...
            max(angleFirst,[],'omitnan');


        floquetErrorAllRMS_deg(iph) = ...
            sqrt(mean(angleAll.^2,'omitnan'));

        floquetErrorAllMax_deg(iph) = ...
            max(angleAll,[],'omitnan');


        normalTracksFirst{iph} = ...
            nQPO(firstMask,:);


        %% Period-to-period non-closure of finite-QPO formation plane

        if nPrimaryPeriods >= 2

            tFirst = ...
                t(firstMask);


            validFirst = ...
                tFirst + TPrimary <= t(end)+1e-10*TPrimary;


            tFirst = ...
                tFirst(validFirst);


            nFirst = ...
                nQPO(firstMask,:);

            nFirst = ...
                nFirst(validFirst,:);


            nSecond = ...
                shifted_normal_track( ...
                    t, ...
                    nQPO, ...
                    tFirst+TPrimary, ...
                    TPrimary);


            closureAngle = ...
                plane_angle_deg( ...
                    nFirst, ...
                    nSecond);


            periodClosureRMS_deg(iph) = ...
                sqrt(mean(closureAngle.^2,'omitnan'));

            periodClosureMax_deg(iph) = ...
                max(closureAngle,[],'omitnan');

        end


        %% Baseline RMS breathing versus calibrated linear prediction

        b12 = vecnorm(B12,2,2);
        b13 = vecnorm(B13,2,2);
        b23 = vecnorm(B23,2,2);

        B_RMS_QPO = ...
            sqrt( ...
                (b12.^2+b13.^2+b23.^2)/3);


        B_RMS_Floquet = ...
            linearScale * ...
            B_RMS_Floquet_unit;


        firstB = ...
            firstMask & ...
            isfinite(B_RMS_QPO) & ...
            isfinite(B_RMS_Floquet) & ...
            B_RMS_Floquet > 0;


        relB = ...
            (B_RMS_QPO(firstB)-B_RMS_Floquet(firstB))./ ...
            B_RMS_Floquet(firstB);


        baselineRelErrorRMS(iph) = ...
            sqrt(mean(relB.^2,'omitnan'));

        baselineRelErrorMax(iph) = ...
            max(abs(relB),[],'omitnan');


        medianBaseline_km(iph) = ...
            median( ...
                [b12;b13;b23], ...
                'omitnan') * ...
            sys.Lstar_km;


        if iph == 1

            nQPOFirst = ...
                nQPO(firstMask,:);

            nFloquetFirst = ...
                nFloquet(firstMask,:);

            % Plane normals are unoriented: n and -n represent the same
            % physical plane.  The comparison metric already uses
            % |n_QPO dot n_F|.  For visualization only, orient the finite
            % QPO normal to the nearest representative of the Floquet
            % normal so coincident planes plot on top of one another.
            nQPOFirstAligned = ...
                align_plane_normal_to_reference( ...
                    nQPOFirst, ...
                    nFloquetFirst);

            referencePhaseTrack.t = ...
                t(firstMask);

            referencePhaseTrack.nQPO = ...
                nQPOFirstAligned;

            referencePhaseTrack.nFloquet = ...
                nFloquetFirst;

            referencePhaseTrack.angle_deg = ...
                angleFirst;

        end

    end


    %% --------------------------------------------------------
    %  Phase spread of finite formation normals
    %  --------------------------------------------------------

    phaseSpreadRMS_deg = ...
        NaN;

    phaseSpreadMax_deg = ...
        NaN;


    if ~isempty(normalTracksFirst)

        nTime = ...
            size(normalTracksFirst{1},1);


        phaseAngles = ...
            nan(nTime,nPhases);


        for it = 1:nTime

            Nmat = ...
                nan(nPhases,3);


            for iph = 1:nPhases

                nTrack = ...
                    normalTracksFirst{iph};


                if size(nTrack,1) >= it

                    Nmat(iph,:) = ...
                        nTrack(it,:);

                end

            end


            valid = ...
                all(isfinite(Nmat),2);


            if nnz(valid) < 1
                continue
            end


            Nuse = ...
                Nmat(valid,:);


            Q = ...
                Nuse.'*Nuse;


            [V,D] = ...
                eig(Q);


            [~,iMax] = ...
                max(diag(D));


            nMean = ...
                V(:,iMax);


            nMean = ...
                nMean/norm(nMean);


            d = ...
                abs(Nuse*nMean);


            d = ...
                min(1,max(-1,d));


            phaseAngles(it,valid) = ...
                acosd(d);

        end


        phaseSpreadRMS_deg = ...
            sqrt( ...
                mean( ...
                    phaseAngles(:).^2, ...
                    'omitnan'));


        phaseSpreadMax_deg = ...
            max( ...
                phaseAngles(:), ...
                [], ...
                'omitnan');

    end


    %% --------------------------------------------------------
    %  Nonlinear rotation-number departure from linear Floquet value
    %  --------------------------------------------------------

    rho0 = ...
        infer_linear_rho_from_floquet_geometry( ...
            floquetGeom);


    rhoError_rad = ...
        min( ...
            abs(wrap_to_pi(rhoQPO-rho0)), ...
            abs(wrap_to_pi(rhoQPO+rho0)));


    %% Package

    out = struct();

    out.phaseStartIndices = phaseStartIndices(:);

    out.sectionRadiusRMS_ND = sectionRadiusRMS_ND;
    out.sectionRadiusRMS_km = sectionRadiusRMS_km;

    out.rhoQPO = rhoQPO;
    out.rhoLinear = rho0;
    out.rhoError_rad = rhoError_rad;

    out.floquetErrorFirstRMS_deg = floquetErrorFirstRMS_deg;
    out.floquetErrorFirstMax_deg = floquetErrorFirstMax_deg;

    out.floquetErrorAllRMS_deg = floquetErrorAllRMS_deg;
    out.floquetErrorAllMax_deg = floquetErrorAllMax_deg;

    out.periodClosureRMS_deg = periodClosureRMS_deg;
    out.periodClosureMax_deg = periodClosureMax_deg;

    out.baselineRelErrorRMS = baselineRelErrorRMS;
    out.baselineRelErrorMax = baselineRelErrorMax;

    out.medianBaseline_km = medianBaseline_km;

    out.qAMinFirstByPhase = qAMinFirstByPhase;

    out.phaseSpreadRMS_deg = phaseSpreadRMS_deg;
    out.phaseSpreadMax_deg = phaseSpreadMax_deg;

    out.referencePhaseTrack = referencePhaseTrack;

    out.summary = struct();

    % E50, E95, and Emax are all statistics of the SAME quantity:
    % the per-phase RMS normal-track error over the first parent period.
    %
    %   E50  : median over common phase starts
    %   E95  : 95th percentile over common phase starts
    %   Emax : worst common phase start
    %
    % This keeps the three metrics directly comparable and avoids mixing an
    % RMS phase statistic with a separate instantaneous maximum.
    out.summary.normalTrackE50_deg = ...
        median(floquetErrorFirstRMS_deg,'omitnan');

    out.summary.normalTrackE95_deg = ...
        percentile_finite(floquetErrorFirstRMS_deg,95);

    out.summary.normalTrackEmax_deg = ...
        max(floquetErrorFirstRMS_deg,[],'omitnan');

    out.summary.qAMinFirst = ...
        min(qAMinFirstByPhase,[],'omitnan');

    % Backward-compatible aliases used by older plotting / analysis code.
    out.summary.floquetErrorFirstMedianRMS_deg = ...
        out.summary.normalTrackE50_deg;

    out.summary.floquetErrorFirstWorstRMS_deg = ...
        out.summary.normalTrackEmax_deg;

    out.summary.floquetErrorFirstWorstMax_deg = ...
        max(floquetErrorFirstMax_deg,[],'omitnan');

    out.summary.floquetErrorAllMedianRMS_deg = ...
        median(floquetErrorAllRMS_deg,'omitnan');

    out.summary.periodClosureMedianRMS_deg = ...
        median(periodClosureRMS_deg,'omitnan');

    out.summary.periodClosureWorstMax_deg = ...
        max(periodClosureMax_deg,[],'omitnan');

    out.summary.baselineRelErrorMedianRMS = ...
        median(baselineRelErrorRMS,'omitnan');

    out.summary.baselineRelErrorWorstMax = ...
        max(baselineRelErrorMax,[],'omitnan');

    out.summary.medianBaseline_km = ...
        median(medianBaseline_km,'omitnan');

end


%% ========================================================================
%  LOCAL HELPERS
%  ========================================================================


function phaseStartIndices = unique_equal_phase_start_indices(N,nCollectors)

    canonicalSets = zeros(N,nCollectors);

    for k = 1:N

        j = qpo_equal_phase_indices(N,nCollectors,k);

        canonicalSets(k,:) = sort(j(:).');

    end

    [~,ia] = unique(canonicalSets,'rows','stable');

    phaseStartIndices = ia(:).';

end


function t = extract_time_track(geom)

    if isfield(geom,'t')

        t = geom.t(:);
        return

    end

    if isfield(geom,'tau')

        t = geom.tau(:);
        return

    end

    error('Could not identify the geometry time vector.');

end


function [B12,B13,B23] = extract_three_baselines_nd(geom)

    B12 = extract_baseline_nd(geom,'12');
    B13 = extract_baseline_nd(geom,'13');
    B23 = extract_baseline_nd(geom,'23');

end


function B = extract_baseline_nd(geom,pairLabel)

    exactNames = { ...
        ['B',pairLabel], ...
        ['b',pairLabel], ...
        ['baseline',pairLabel], ...
        ['B',pairLabel,'Synodic'], ...
        ['b',pairLabel,'Synodic'], ...
        ['baseline',pairLabel,'Synodic']};

    for k = 1:numel(exactNames)

        if isfield(geom,exactNames{k})

            X = geom.(exactNames{k});

            if is_xyz_history(X)

                B = orient_xyz_history(X);
                return

            end

        end

    end


    candidates = collect_xyz_fields(geom,'geom');

    bestScore = -Inf;
    bestValue = [];

    pLabel = lower(pairLabel);

    for k = 1:numel(candidates)

        pathLower = lower(candidates(k).path);

        if contains(pathLower,'km') || ...
           contains(pathLower,'uvw') || ...
           contains(pathLower,'velocity')
            continue
        end

        if ~(contains(pathLower,['b',pLabel]) || ...
             contains(pathLower,['baseline',pLabel]))
            continue
        end

        score = 0;

        if contains(pathLower,'baseline')
            score = score + 5;
        end

        if contains(pathLower,['b',pLabel])
            score = score + 5;
        end

        if contains(pathLower,'synodic')
            score = score + 2;
        end

        if score > bestScore

            bestScore = score;
            bestValue = orient_xyz_history(candidates(k).value);

        end

    end

    if isempty(bestValue)

        error( ...
            ['Could not identify baseline ',pairLabel, ...
             ' in geometry structure.']);

    end

    B = bestValue;

end


function n = normal_from_baselines(B12,B13)

    nRaw = cross(B12,B13,2);

    nNorm = vecnorm(nRaw,2,2);

    valid = nNorm > 0 & isfinite(nNorm);

    n = nan(size(nRaw));

    n(valid,:) = ...
        nRaw(valid,:)./nNorm(valid);

    n = enforce_sign_continuity(n);

end


function [n,B_RMS_unit] = direct_floquet_predictor_at_times( ...
    po,sys,fg,t,odeOpts)
%DIRECT_FLOQUET_PREDICTOR_AT_TIMES
% Integrate the parent STM at the exact parent-orbit phases used by the
% finite-QPO propagation.

    T = po.T;

    phaseTime = ...
        mod(t-t(1),T);

    values = ...
        [0;phaseTime(:)];

    [tEval,~,mapAll] = ...
        unique(values,'sorted');

    map = ...
        mapAll(2:end);


    Y0 = [ ...
        reshape(eye(6),36,1); ...
        po.x0(:)];


    [tReturned,YReturned] = ...
        ode89( ...
            @(tau,Y) STM_vec(tau,Y,sys.mu), ...
            tEval, ...
            Y0, ...
            odeOpts);


    tReturned = ...
        tReturned(:);


    % When tspan contains more than two entries, the two-output ODE
    % interface must return one state row at each requested tEval value.
    % Check this explicitly so a solver-interface mismatch can never turn
    % into a misleading array-bounds error.
    if numel(tReturned) ~= numel(tEval) || ...
            size(YReturned,1) ~= numel(tEval)

        error( ...
            ['Exact-phase Floquet propagation returned %d time samples ', ...
             'and %d state rows for %d requested times.'], ...
            numel(tReturned), ...
            size(YReturned,1), ...
            numel(tEval));

    end


    timeGridError = ...
        max(abs(tReturned-tEval(:)));


    if timeGridError > 1e-11*max(1,T)

        error( ...
            ['Exact-phase Floquet propagation did not return the requested ', ...
             'time grid.  max |tReturned-tEval| = %.3e.'], ...
            timeGridError);

    end


    w0 = ...
        fg.wTrack(1,:).';


    nEval = nan(numel(tEval),3);
    B_RMS_eval = nan(numel(tEval),1);


    for k = 1:numel(tEval)

        Phi = ...
            reshape(YReturned(k,1:36),[6,6]);

        wk = ...
            Phi*w0;

        a = real(wk(1:3));
        b = imag(wk(1:3));

        nRaw = cross(a,b);
        nNorm = norm(nRaw);

        if isfinite(nNorm) && nNorm > 0

            nEval(k,:) = ...
                (nRaw/nNorm).';

        end

        A = ...
            [a,-b];

        R_F = ...
            norm(A,'fro')/sqrt(2);

        B_RMS_eval(k) = ...
            sqrt(3)*R_F;

    end


    nEval = ...
        enforce_sign_continuity(nEval);


    n = ...
        nEval(map,:);

    B_RMS_unit = ...
        B_RMS_eval(map);

end


function nOut = shifted_normal_track(t,n,tq,T)
%SHIFTED_NORMAL_TRACK Use exact one-period paired samples when available.

    nOut = nan(numel(tq),3);

    timeTol = ...
        1e-10*max(1,T);


    for k = 1:numel(tq)

        [dMin,j] = ...
            min(abs(t-tq(k)));

        if isfinite(dMin) && dMin <= timeTol

            nOut(k,:) = ...
                n(j,:);

        end

    end


    missing = ...
        ~all(isfinite(nOut),2);


    if any(missing)

        nInterp = ...
            interpolate_normal_track( ...
                t, ...
                n, ...
                tq(missing));

        nOut(missing,:) = ...
            nInterp;

    end

end


function n = periodic_floquet_normal_at_times(fg,t,T)

    tau = mod(t-t(1),T)/T;

    % Make exact period endpoints map to tau=1 rather than zero only when
    % they are numerically at the first-period endpoint.  Since the plane is
    % periodic either representation is equivalent.
    n = nan(numel(t),3);

    for j = 1:3

        n(:,j) = ...
            interp1( ...
                fg.tau, ...
                fg.nSynodic(:,j), ...
                tau, ...
                'pchip');

    end

    n = normalize_rows(n);

    n = enforce_sign_continuity(n);

end


function y = periodic_scalar_at_times(tBase,yBase,t,T)

    tauBase = ...
        (tBase-tBase(1))/(tBase(end)-tBase(1));

    tau = ...
        mod(t-t(1),T)/T;

    y = ...
        interp1( ...
            tauBase, ...
            yBase, ...
            tau, ...
            'pchip');

end


function nOut = interpolate_normal_track(t,n,tq)

    nOut = nan(numel(tq),3);

    for j = 1:3

        nOut(:,j) = ...
            interp1( ...
                t, ...
                n(:,j), ...
                tq, ...
                'pchip', ...
                NaN);

    end

    nOut = normalize_rows(nOut);

end


function a = plane_angle_deg(n1,n2)

    d = sum(n1.*n2,2);

    d = abs(d);

    d = min(1,max(-1,d));

    a = acosd(d);

end


function n = normalize_rows(n)

    r = vecnorm(n,2,2);

    valid = r > 0 & isfinite(r);

    n(valid,:) = n(valid,:)./r(valid);

    n(~valid,:) = NaN;

end


function n = enforce_sign_continuity(n)

    validIdx = find(all(isfinite(n),2));

    for j = 2:numel(validIdx)

        iPrev = validIdx(j-1);
        iNow = validIdx(j);

        if dot(n(iPrev,:),n(iNow,:)) < 0

            n(iNow,:) = -n(iNow,:);

        end

    end

end


function rho0 = infer_linear_rho_from_floquet_geometry(fg)
% The transported helper estimates lambda from w(T)=lambda*w(0).

    lambda = fg.lambdaEstimated;

    rho0 = angle(lambda);

end


function nAligned = align_plane_normal_to_reference(n,nRef)
%ALIGN_PLANE_NORMAL_TO_REFERENCE
% Pointwise sign alignment for plotting unoriented plane normals.
%
% This does not change the physical plane or any sign-invariant metric.

    nAligned = n;

    valid = ...
        all(isfinite(n),2) & ...
        all(isfinite(nRef),2);

    idx = find(valid);

    for k = idx(:).'

        if dot(nAligned(k,:),nRef(k,:)) < 0

            nAligned(k,:) = ...
                -nAligned(k,:);

        end

    end

end


function x = wrap_to_pi(x)

    x = mod(x+pi,2*pi)-pi;

end


function tf = is_xyz_history(X)

    tf = ...
        isnumeric(X) && ...
        ismatrix(X) && ...
        ~isempty(X) && ...
        (size(X,2)==3 || size(X,1)==3);

end


function X = orient_xyz_history(X)

    if size(X,2)==3
        return
    end

    if size(X,1)==3

        X = X.';
        return

    end

    error('Expected N-by-3 or 3-by-N array.');

end


function candidates = collect_xyz_fields(S,rootName)

    candidates = struct('path',{},'value',{});

    if ~isstruct(S)
        return
    end

    f = fieldnames(S);

    for k = 1:numel(f)

        name = f{k};
        value = S.(name);

        path = [rootName,'.',name];

        if is_xyz_history(value)

            candidates(end+1).path = path; %#ok<AGROW>
            candidates(end).value = value;

        elseif isstruct(value) && isscalar(value)

            child = collect_xyz_fields(value,path);

            if ~isempty(child)

                candidates = [candidates,child]; %#ok<AGROW>

            end

        end

    end

end


function q = percentile_finite(x,p)
%PERCENTILE_FINITE Toolbox-free linear percentile of finite samples.

    x = x(isfinite(x));
    x = sort(x(:));

    if isempty(x)
        q = NaN;
        return
    end

    if numel(x)==1
        q = x;
        return
    end

    p = min(100,max(0,p));

    s = 1 + (numel(x)-1)*(p/100);

    i0 = floor(s);
    i1 = ceil(s);

    if i0==i1
        q = x(i0);
    else
        a = s-i0;
        q = (1-a)*x(i0) + a*x(i1);
    end

end
