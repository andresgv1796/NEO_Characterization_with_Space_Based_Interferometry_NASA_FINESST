function famOut = build_gmos_fixedT_amplitude_family( ...
    po, sys, N, KSeeds, mode, nContinuationSteps, ...
    continuationStepMultiplier, odeOpts, newtonOpts)
%BUILD_GMOS_FIXEDT_AMPLITUDE_FAMILY
% Build a finite-amplitude QPO family from one parent periodic orbit using
% two nearby corrected seeds followed by fixed-T pseudo-arclength
% continuation.
%
% This is specifically intended for finite-amplitude breakdown studies.
% The continuation family is preferable to solving many independent seed
% amplitudes because independent Newton solves can jump to unrelated
% branches.
%
% Inputs
%   po                         parent periodic-orbit structure
%   sys                        CR3BP system structure
%   N                          invariant-curve discretization
%   KSeeds                     two small seed amplitudes [K0 K1]
%   mode                       'planar' or 'vertical'
%   nContinuationSteps         steps after the two corrected seeds
%   continuationStepMultiplier ds / ||z1-z0||
%   odeOpts                    ODE options
%   newtonOpts                 Newton options
%
% Output
%   famOut.fam                 raw gmos_continue_fixT_family output
%   famOut.radiusRMS_ND        actual corrected section RMS radius
%   famOut.radiusRMS_km        same, dimensional
%   famOut.radiusMax_km
%   famOut.jacobiSpread
%   famOut.rho
%   famOut.T
%   famOut.ds
%   famOut.seedInfo
%
% Required existing functions
%   gmos_demo_invariance_only
%   gmos_pack_z
%   gmos_unpack_z
%   gmos_newton_fixT_analytic
%   gmos_continue_fixT_family
%   Jacobi_Constant

    validateattributes(KSeeds,{'numeric'}, ...
        {'vector','numel',2,'positive','finite'});

    if KSeeds(2) <= KSeeds(1)
        error('KSeeds must satisfy KSeeds(2) > KSeeds(1).');
    end

    mu = sys.mu;

    xPO0 = po.x0(:);
    Ttarget = po.T;

    %% --------------------------------------------------------
    %  Build and strictly correct the first startup member
    %  --------------------------------------------------------
    %
    % A universal K0 is not equally easy for every parent halo.  The first
    % corrected member is used only to initialize PALC, so if the requested
    % K0 stalls above tolInf we contract K0 toward the parent periodic orbit.
    %
    % IMPORTANT:
    %   - tolInf is NEVER relaxed;
    %   - a seed is accepted only when the existing fixed-T Newton corrector
    %     meets the same strict convergence criterion used elsewhere;
    %   - the requested and actually used amplitudes are recorded.

    startupNewtonOpts = newtonOpts;

    % Allow additional iterations for startup only.  This does not change
    % the Newton residual tolerance.
    startupNewtonOpts.maxIter = ...
        max(newtonOpts.maxIter,40);


    K0Requested = KSeeds(1);

    firstSeedScales = ...
        [1.0 0.75 0.50 0.35 0.25 0.15 0.10 0.05];

    firstSeedConverged = false;

    z0Sol = [];
    info0 = struct();

    X0Sol = [];
    T0Sol = NaN;
    rho0Sol = NaN;

    K0Used = NaN;

    firstSeedAttempts = struct( ...
        'K',{}, ...
        'finalNormInf',{}, ...
        'radiusInit_ND',{}, ...
        'radiusCorrected_ND',{}, ...
        'localityRatio',{}, ...
        'converged',{});


    for iTry = 1:numel(firstSeedScales)

        K0Try = ...
            K0Requested*firstSeedScales(iTry);


        demo0Try = ...
            gmos_demo_invariance_only( ...
                mu, ...
                xPO0, ...
                Ttarget, ...
                N, ...
                K0Try, ...
                false, ...
                mode);


        init0Try = ...
            demo0Try.init;


        z0Try = ...
            gmos_pack_z( ...
                init0Try.X0, ...
                init0Try.T0, ...
                init0Try.rho0);


        [z0TrySol,info0Try] = ...
            gmos_newton_fixT_analytic( ...
                z0Try, ...
                N, ...
                mu, ...
                init0Try.X0, ...
                init0Try.T0, ...
                init0Try.rho0, ...
                Ttarget, ...
                odeOpts, ...
                startupNewtonOpts);


        if isfield(info0Try,'normInf') && ...
                ~isempty(info0Try.normInf)

            finalNormInf = ...
                info0Try.normInf(end);

        else

            finalNormInf = ...
                Inf;

        end


        [X0TrySol,T0TrySol,rho0TrySol] = ...
            gmos_unpack_z( ...
                z0TrySol, ...
                N);

        rInitTry = ...
            section_rms_radius_nd(init0Try.X0);

        rCorrectedTry = ...
            section_rms_radius_nd(X0TrySol);

        localityRatio = ...
            symmetric_positive_ratio( ...
                rCorrectedTry, ...
                rInitTry);

        localTry = ...
            isfinite(localityRatio) && ...
            localityRatio <= 10;

        convergedTry = ...
            isfinite(finalNormInf) && ...
            finalNormInf < startupNewtonOpts.tolInf && ...
            localTry;


        firstSeedAttempts(end+1) = struct( ... %#ok<AGROW>
            'K',K0Try, ...
            'finalNormInf',finalNormInf, ...
            'radiusInit_ND',rInitTry, ...
            'radiusCorrected_ND',rCorrectedTry, ...
            'localityRatio',localityRatio, ...
            'converged',convergedTry);


        fprintf( ...
            ['First startup seed attempt: K = %.6e, ', ...
             'final ||F||inf = %.3e, locality ratio = %.3g%s\n'], ...
            K0Try, ...
            finalNormInf, ...
            localityRatio, ...
            ternary_text(convergedTry,'  [accepted]',''));


        if convergedTry

            z0Sol = ...
                z0TrySol;

            info0 = ...
                info0Try;

            K0Used = ...
                K0Try;

            X0Sol = ...
                X0TrySol;

            T0Sol = ...
                T0TrySol;

            rho0Sol = ...
                rho0TrySol;

            firstSeedConverged = ...
                true;

            break

        end

    end


    if ~firstSeedConverged

        attemptedK = ...
            [firstSeedAttempts.K];

        attemptedF = ...
            [firstSeedAttempts.finalNormInf];

        error( ...
            ['Unable to obtain a first strict fixed-T startup member. ', ...
             'Requested K0 = %.3e; tried K = [%s], with final ', ...
             '||F||inf = [%s].  PALC was not started.'], ...
            K0Requested, ...
            num2str(attemptedK,' %.3e'), ...
            num2str(attemptedF,' %.3e'));

    end


    %% --------------------------------------------------------
    %  Correct second seed using first corrected member as gauge reference
    %  --------------------------------------------------------
    %
    % A single universal KSeeds(2) is not equally easy for every parent
    % halo.  Some center modes are much more poorly conditioned, so a seed
    % that works for one parent can stall a few orders above tolInf for
    % another.  This is a STARTUP issue, not a failure of the QPO family.
    %
    % Keep the requested second seed first.  If it does not converge to the
    % SAME strict tolerance, contract only the seed separation toward K0.
    % We never relax tolInf and we never accept the 1e-8-ish residual merely
    % because it is "close".

    K0 = K0Used;
    K1Requested = KSeeds(2);

    seedFractions = [1.0 0.75 0.50 0.35 0.25 0.15 0.10 0.05];

    secondSeedConverged = false;

    z1Sol = [];
    info1 = struct();
    K1Used = NaN;

    secondSeedAttempts = struct( ...
        'K',{}, ...
        'finalNormInf',{}, ...
        'radiusInit_ND',{}, ...
        'radiusCorrected_ND',{}, ...
        'localityRatio',{}, ...
        'converged',{});

    for iTry = 1:numel(seedFractions)

        K1Try = ...
            K0 + seedFractions(iTry)*(K1Requested-K0);

        % Ensure the second startup member is genuinely distinct.
        if K1Try <= K0*(1+1e-8)
            continue
        end

        demo1Try = ...
            gmos_demo_invariance_only( ...
                mu, ...
                xPO0, ...
                Ttarget, ...
                N, ...
                K1Try, ...
                false, ...
                mode);

        init1Try = demo1Try.init;

        z1Try = ...
            gmos_pack_z( ...
                init1Try.X0, ...
                init1Try.T0, ...
                init1Try.rho0);

        [z1TrySol,info1Try] = ...
            gmos_newton_fixT_analytic( ...
                z1Try, ...
                N, ...
                mu, ...
                X0Sol, ...
                T0Sol, ...
                rho0Sol, ...
                Ttarget, ...
                odeOpts, ...
                startupNewtonOpts);

        if isfield(info1Try,'normInf') && ...
                ~isempty(info1Try.normInf)

            finalNormInf = ...
                info1Try.normInf(end);

        else

            finalNormInf = Inf;

        end

        [X1TrySol,~,~] = ...
            gmos_unpack_z( ...
                z1TrySol, ...
                N);

        rInitTry = ...
            section_rms_radius_nd(init1Try.X0);

        rCorrectedTry = ...
            section_rms_radius_nd(X1TrySol);

        localityRatio = ...
            symmetric_positive_ratio( ...
                rCorrectedTry, ...
                rInitTry);

        localTry = ...
            isfinite(localityRatio) && ...
            localityRatio <= 10;

        convergedTry = ...
            isfinite(finalNormInf) && ...
            finalNormInf < startupNewtonOpts.tolInf && ...
            localTry;

        secondSeedAttempts(end+1) = struct( ... %#ok<AGROW>
            'K',K1Try, ...
            'finalNormInf',finalNormInf, ...
            'radiusInit_ND',rInitTry, ...
            'radiusCorrected_ND',rCorrectedTry, ...
            'localityRatio',localityRatio, ...
            'converged',convergedTry);

        fprintf( ...
            ['Second startup seed attempt: K = %.6e, ', ...
             'final ||F||inf = %.3e, locality ratio = %.3g%s\n'], ...
            K1Try, ...
            finalNormInf, ...
            localityRatio, ...
            ternary_text(convergedTry,'  [accepted]',''));

        if convergedTry

            z1Sol = z1TrySol;
            info1 = info1Try;
            K1Used = K1Try;

            secondSeedConverged = true;
            break

        end

    end

    if ~secondSeedConverged

        attemptedK = [secondSeedAttempts.K];
        attemptedF = [secondSeedAttempts.finalNormInf];

        error( ...
            ['Unable to obtain a second strict fixed-T startup member. ', ...
             'Requested K1 = %.3e; tried K = [%s], with final ', ...
             '||F||inf = [%s].  The continuation itself was not started.'], ...
            K1Requested, ...
            num2str(attemptedK,' %.3e'), ...
            num2str(attemptedF,' %.3e'));

    end

    [X1Sol,T1Sol,rho1Sol] = ...
        gmos_unpack_z(z1Sol,N);


    %% --------------------------------------------------------
    %  Ensure the secant points toward increasing physical section size
    %  --------------------------------------------------------

    r0 = section_rms_radius_nd(X0Sol);
    r1 = section_rms_radius_nd(X1Sol);

    swappedSeeds = false;

    if r1 < r0

        tmp = z0Sol;
        z0Sol = z1Sol;
        z1Sol = tmp;

        tmp = r0;
        r0 = r1;
        r1 = tmp;

        swappedSeeds = true;

    end

    secantNorm = norm(z1Sol-z0Sol);

    if ~(isfinite(secantNorm) && secantNorm > 0)
        error('Corrected startup members are identical or non-finite.');
    end

    ds = continuationStepMultiplier*secantNorm;


    %% --------------------------------------------------------
    %  Pseudo-arclength continuation at fixed parent period
    %  --------------------------------------------------------

    fam = gmos_continue_fixT_family( ...
        z0Sol, ...
        z1Sol, ...
        nContinuationSteps, ...
        ds, ...
        N, ...
        mu, ...
        Ttarget, ...
        odeOpts, ...
        newtonOpts);


    %% --------------------------------------------------------
    %  Physical diagnostics for every accepted family member
    %  --------------------------------------------------------

    nFam = numel(fam.z);

    radiusRMS_ND = nan(nFam,1);
    radiusRMS_km = nan(nFam,1);
    radiusMax_km = nan(nFam,1);

    jacobiSpread = nan(nFam,1);

    rho = nan(nFam,1);
    T = nan(nFam,1);

    for k = 1:nFam

        [XCurve,Tk,rhok] = ...
            gmos_unpack_z(fam.z{k},N);

        [rRMS,rMax] = ...
            section_radius_nd(XCurve);

        radiusRMS_ND(k) = rRMS;
        radiusRMS_km(k) = rRMS*sys.Lstar_km;
        radiusMax_km(k) = rMax*sys.Lstar_km;

        C = Jacobi_Constant(XCurve.',mu);

        jacobiSpread(k) = ...
            max(C)-min(C);

        rho(k) = rhok;
        T(k) = Tk;

    end


    %% --------------------------------------------------------
    %  Resolution / continuation-quality mask
    %  --------------------------------------------------------
    %
    % A Newton-converged large-amplitude invariant curve can still be
    % under-resolved at fixed N.  Use the already-computed Jacobi spread as
    % the resolution diagnostic.

    nRef = min(5,nFam);

    jacobiReference = ...
        median( ...
            jacobiSpread(1:nRef), ...
            'omitnan');

    jacobiLimit = ...
        max( ...
            1e-10, ...
            1e4*jacobiReference);

    analysisValid = ...
        isfinite(radiusRMS_km) & ...
        radiusRMS_km > 0 & ...
        isfinite(jacobiSpread) & ...
        jacobiSpread <= jacobiLimit;


    %% --------------------------------------------------------
    %  Package
    %  --------------------------------------------------------

    famOut = struct();

    famOut.fam = fam;

    famOut.radiusRMS_ND = radiusRMS_ND;
    famOut.radiusRMS_km = radiusRMS_km;
    famOut.radiusMax_km = radiusMax_km;

    famOut.jacobiSpread = jacobiSpread;

    famOut.analysisValid = analysisValid;
    famOut.jacobiReference = jacobiReference;
    famOut.jacobiLimit = jacobiLimit;

    famOut.rho = rho;
    famOut.T = T;

    famOut.ds = ds;

    famOut.seedInfo = struct();

    famOut.seedInfo.KSeedsRequested = KSeeds(:).';
    famOut.seedInfo.KSeedsUsed = [K0Used,K1Used];

    famOut.seedInfo.firstSeedAttempts = firstSeedAttempts;
    famOut.seedInfo.secondSeedAttempts = secondSeedAttempts;

    famOut.seedInfo.swappedSeeds = swappedSeeds;

    famOut.seedInfo.initialCorrectedRadiusRMS_ND = ...
        [r0,r1];

    famOut.seedInfo.info0 = info0;
    famOut.seedInfo.info1 = info1;

end


function r = symmetric_positive_ratio(a,b)
%SYMMETRIC_POSITIVE_RATIO Return max(a/b,b/a) for positive finite scalars.

    if ~(isfinite(a) && isfinite(b) && a > 0 && b > 0)

        r = Inf;
        return

    end

    r = max(a/b,b/a);

end


function txt = ternary_text(condition,txtTrue,txtFalse)
%TERNARY_TEXT Small local utility for compact diagnostic printing.

    if condition
        txt = txtTrue;
    else
        txt = txtFalse;
    end

end


function r = section_rms_radius_nd(XCurve)
%SECTION_RMS_RADIUS_ND

    r = section_radius_nd(XCurve);

end


function [rRMS,rMax] = section_radius_nd(XCurve)
%SECTION_RADIUS_ND
% RMS and maximum configuration-space radius about the section centroid.

    rCenter = mean(XCurve(1:3,:),2);

    dR = XCurve(1:3,:) - rCenter;

    rSamples = vecnorm(dR,2,1);

    rRMS = sqrt(mean(rSamples.^2));
    rMax = max(rSamples);

end
