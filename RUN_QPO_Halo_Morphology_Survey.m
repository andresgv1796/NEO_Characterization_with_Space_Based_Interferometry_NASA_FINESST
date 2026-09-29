%% RUN_QPO_HALO_MORPHOLOGY_SURVEY
%
% FINESST QPO FORMATION-MORPHOLOGY SURVEY
%
% Purpose
% -------
% This script studies how the geometry and interferometric behavior of
% three-collector QPO formations vary with:
%
%   1. parent northern-halo orbit selected from the JPL catalog,
%   2. corrected QPO torus size,
%   3. common starting phase of the three equally phased collectors.
%
% The script is deliberately a MORPHOLOGY / SCREENING study rather than a
% full asteroid-reconstruction study.  It answers:
%
%   A. Does the array-plane normal history depend strongly on torus size?
%   B. Does the array-plane normal history depend strongly on common phase?
%   C. Does torus size mainly rescale the baselines, or does it change the
%      normalized formation morphology?
%   D. How strongly does common phase change the baseline triangle?
%   E. How do those changes map into target accessibility and UV coverage?
%
% Study organization
% ------------------
% For each representative halo parent:
%
%   - several corrected GMOS QPOs are generated at different seed amplitudes;
%   - the ACTUAL corrected section radius is measured from XCurve;
%   - all physically distinct equal-phase three-collector triangles are flown;
%   - target-independent formation quantities are measured:
%
%         normal-vector track
%         baseline lengths
%         triangle area
%         triangle aspect ratio
%
%   - one reference formation (middle successful torus size, first unique
%     common phase) is used to select ONE all-sky target for that parent;
%   - that same target is then held fixed across every torus-size / phase case
%     for that parent, so changes in UV coverage are caused by the formation
%     rather than by changing the target.
%
% Interferometric screening quantities
% ------------------------------------
% The script records:
%
%   accessible fraction
%   maximum and RMS projected spatial frequency
%   Fourier-plane covariance isotropy
%   angular-sector fill
%   radial-bin fill
%   normalized UV-cell occupancy
%   normalized convex-hull area
%   fringe-spacing scale
%   covariance-based beam-axis-ratio proxy
%
% It does NOT run Itokawa, dirty imaging, or CLEAN for every case.  The point
% is to identify the formation families and parameter regions that deserve
% the expensive end-to-end reconstruction study.
%
% Required existing FINESST helpers
% ---------------------------------
%   cr3bp_system_parameters
%   load_JPL_orbit_data
%   choose_periodic_orbit
%   generate_gmos_qpo_torus_from_po
%   qpo_equal_phase_indices
%   compute_three_collector_qpo_geometry
%   add_dimensional_geometry
%   sample_sphere_targets
%   evaluate_target_suitability_scan
%   make_target_definition
%   compute_triangle_accessibility
%   compute_fixed_target_uv
%
% The local helper routines at the bottom intentionally tolerate several
% plausible field names in geom so the study can be used with the current
% FINESST geometry structure without hard-coding one fragile spelling.

clear;
clc;
close all;


%% ============================================================
%  0. PLOT DEFAULTS
%  ============================================================

set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

set(groot,'DefaultAxesFontSize',15);
set(groot,'DefaultTextFontSize',17);
set(groot,'DefaultLineLineWidth',1.25);

set(groot,'DefaultFigureColor','w');
set(groot,'DefaultAxesColor','w');
set(groot,'DefaultAxesXColor','k');
set(groot,'DefaultAxesYColor','k');
set(groot,'DefaultTextColor','k');


%% ============================================================
%  1. USER CONFIGURATION
%  ============================================================

%% CR3BP system

systemName = 'EarthMoon';

sys = cr3bp_system_parameters(systemName);
mu = sys.mu;


%% JPL halo catalog

poFile = 'L1_northern_halo_orbits_JPL_IC_Earth_Moon.csv';

% Fractions of the physical z0-amplitude range used to choose representative
% parents.  These avoid the literal family endpoints.
haloFamilyFractions = [0.05 0.10 0.15 0.20 0.25 0.30];


%% GMOS QPO settings

qpoGenBase = struct();

qpoGenBase.N = 51;
qpoGenBase.mode = 'vertical';

qpoGenBase.corrector = 'fixT';
qpoGenBase.CtargetMode = 'parent';

qpoGenBase.maxIter = 15;
qpoGenBase.tolInf = 1e-10;

% The survey would otherwise print thousands of Newton lines.
qpoGenBase.verbose = false;
qpoGenBase.plotInitialGuess = false;
qpoGenBase.saveGeneratedQPO = false;


% Seed amplitudes used to generate tori.
%
% IMPORTANT:
% These are NOT treated as the physical torus-size coordinate in the plots.
% The script measures the actual corrected invariant-curve section radius.
%
% After confirming convergence, you can append 1e-3 to probe deeper into the
% nonlinear finite-amplitude regime.
torusSeedAmplitudes = [ ...
    1e-5, ...
    3e-5, ...
    1e-4, ...
    3e-4];


%% Three-collector formation

nCollectors = 3;

normalReferenceSynodic = [1 0 0];


%% Propagation horizon

nPrimaryPeriods = 3;
nSavePerPrimaryPeriod = 250;

timePlotUnit = 'auto';

odeOpts = odeset( ...
    'RelTol',1e-12, ...
    'AbsTol',1e-12);


%% Interferometer

lambda_m = 10e-6;


%% Accessibility / target study

maxOffNormalDeg = 10.0;
usePlaneNormalSignAmbiguity = true;

skyNorthReference = [0 0 1];

% One target scan is run for each representative halo parent, using that
% parent's reference formation.
nTargetsPerParentScan = 4000;
useParallelTargetScan = true;


%% UV morphology bins

nAngularBins = 18;
nRadialBins = 12;
nUVCellBins = 30;


%% Plot only a few representative phases in the detailed UV mosaics.
%
% All unique phases are still propagated and included in the numerical
% metrics and heat maps.
nPhasesToPlot = 4;


%% Saving

saveResults = true;
saveFigures = true;

outputDir = fullfile( ...
    pwd, ...
    'QPO_Halo_Morphology_Survey');

if ~exist(outputDir,'dir')
    mkdir(outputDir);
end


%% ============================================================
%  2. LOAD JPL CATALOG AND SELECT REPRESENTATIVE HALO PARENTS
%  ============================================================

[ICCatalog,TPrimaryCatalog] = ...
    load_JPL_orbit_data(poFile);


if size(ICCatalog,2) ~= 6

    error( ...
        'Expected JPL initial-condition matrix to have six state columns.');

end


nCatalog = ...
    size(ICCatalog,1);


catalogIndex = ...
    (1:nCatalog).';


% Northern-halo family coordinate used for representative selection.
%
% The JPL helper loads [x y z vx vy vz], so column 3 is z0.
z0Catalog = ...
    ICCatalog(:,3);


haloAmplitudeProxy = ...
    abs(z0Catalog);


validCatalog = ...
    all(isfinite(ICCatalog),2) & ...
    isfinite(TPrimaryCatalog) & ...
    isfinite(haloAmplitudeProxy);


validIdx = ...
    find(validCatalog);


if numel(validIdx) < numel(haloFamilyFractions)

    error( ...
        'Not enough valid halo entries in the JPL catalog.');

end


zValid = ...
    haloAmplitudeProxy(validIdx);


zMin = ...
    min(zValid);

zMax = ...
    max(zValid);


representativeOrbitIndices = ...
    zeros(numel(haloFamilyFractions),1);


representativeZ0 = ...
    zeros(numel(haloFamilyFractions),1);


for k = 1:numel(haloFamilyFractions)

    zTarget = ...
        zMin + ...
        haloFamilyFractions(k)*(zMax-zMin);


    [~,iLocal] = ...
        min(abs(zValid-zTarget));


    representativeOrbitIndices(k) = ...
        validIdx(iLocal);


    representativeZ0(k) = ...
        z0Catalog(representativeOrbitIndices(k));

end


% Remove accidental duplicate catalog rows if the family is strongly
% clustered in z0.
representativeOrbitIndices = ...
    unique( ...
        representativeOrbitIndices, ...
        'stable');


nParents = ...
    numel(representativeOrbitIndices);


representativeZ0 = ...
    z0Catalog(representativeOrbitIndices);


representativePeriods = ...
    TPrimaryCatalog(representativeOrbitIndices);


parentSelectionTable = ...
    table( ...
        (1:nParents).', ...
        representativeOrbitIndices, ...
        representativeZ0, ...
        abs(representativeZ0), ...
        representativePeriods, ...
        'VariableNames',{ ...
            'ParentCase', ...
            'CatalogRow', ...
            'z0', ...
            'AbsZ0', ...
            'PeriodND'});


fprintf('\n');
fprintf('============================================================\n');
fprintf('REPRESENTATIVE JPL HALO PARENTS\n');
fprintf('============================================================\n');

disp(parentSelectionTable);


%% Catalog-selection figure

figCatalog = ...
    figure( ...
        'Name','Representative JPL halo parents', ...
        'Color','w', ...
        'Position',[100 100 1100 500]);


plot( ...
    catalogIndex(validCatalog), ...
    haloAmplitudeProxy(validCatalog), ...
    'k.', ...
    'MarkerSize',7);

hold on;


scatter( ...
    representativeOrbitIndices, ...
    abs(representativeZ0), ...
    80, ...
    (1:nParents).', ...
    'filled');


grid on;

xlabel('JPL catalog row');
ylabel('$|z_0|$');

title( ...
    'Representative L1 northern-halo parents selected across the JPL family');


for ip = 1:nParents

    text( ...
        representativeOrbitIndices(ip), ...
        abs(representativeZ0(ip)), ...
        sprintf('  P%d',ip), ...
        'VerticalAlignment','bottom');

end


%% ============================================================
%  3. SURVEY STORAGE
%  ============================================================

nSizes = ...
    numel(torusSeedAmplitudes);


% Phase count is determined after the first successful QPO.
phaseStartIndices = [];
nPhaseCases = [];


%% QPO-level storage

qpoSuccess = ...
    false(nParents,nSizes);

qpoFailureMessage = ...
    cell(nParents,nSizes);

qpoStore = ...
    cell(nParents,nSizes);

actualSectionRadiusRMS_km = ...
    nan(nParents,nSizes);

actualSectionRadiusMax_km = ...
    nan(nParents,nSizes);

sectionCentroidOffset_km = ...
    nan(nParents,nSizes);

qpoRho = ...
    nan(nParents,nSizes);

qpoT = ...
    nan(nParents,nSizes);


%% Parent-level target storage

targetRA_deg = ...
    nan(nParents,1);

targetDec_deg = ...
    nan(nParents,1);

targetScore = ...
    nan(nParents,1);

targetAccessibleFractionReference = ...
    nan(nParents,1);

targetUVIsotropyReference = ...
    nan(nParents,1);

targetRhoMaxReference_Glambda = ...
    nan(nParents,1);

referenceSizeIndex = ...
    nan(nParents,1);

referencePhaseCase = ...
    ones(nParents,1);


%% These case arrays are allocated after the number of unique phases is known.

normalRMSDifference_deg = [];
normalMaxDifference_deg = [];

baselineMedian_km = [];
baselineRMS_km = [];
baselineMax_km = [];

triangleAreaMedian_km2 = [];
triangleAspectMedian = [];

baselineToTorusRatio = [];
areaToTorusRadius2 = [];

accessibleFraction = [];

rhoMax_Glambda = [];
rhoRMS_Glambda = [];

uvIsotropy = [];
uvAngularFill = [];
uvRadialFill = [];
uvCellOccupancy = [];
uvHullFraction = [];

fringeSpacing_mas = [];
beamAxisRatioProxy = [];

caseSuccess = [];
caseFailureMessage = [];

normalTrackStore = {};
timeTrackStore = {};
uvPointStore = {};


%% ============================================================
%  4. MAIN HALO-PARENT SURVEY
%  ============================================================

for ip = 1:nParents

    orbitIndex = ...
        representativeOrbitIndices(ip);


    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('PARENT %d / %d : JPL ROW %d\n', ...
        ip,nParents,orbitIndex);
    fprintf('============================================================\n');


    po = ...
        choose_periodic_orbit( ...
            poFile, ...
            orbitIndex, ...
            sys);


    xCombiner0 = ...
        po.x0(:);


    TPrimary = ...
        po.T;


    fprintf('  |z0|       = %.8e\n',abs(po.x0(3)));
    fprintf('  TPrimary   = %.8e nondim\n',TPrimary);
    fprintf('  TPrimary   = %.6f days\n',po.T_days);


    %% --------------------------------------------------------
    %  4.1 GENERATE THE TORUS-SIZE FAMILY FOR THIS PARENT
    %  --------------------------------------------------------

    for is = 1:nSizes

        qpoGen = ...
            qpoGenBase;


        qpoGen.initialAmplitude = ...
            torusSeedAmplitudes(is);


        fprintf('\n');
        fprintf('  Generating QPO size %d/%d: seed K = %.3e\n', ...
            is,nSizes,qpoGen.initialAmplitude);


        try

            qpo = ...
                generate_gmos_qpo_torus_from_po( ...
                    po, ...
                    sys, ...
                    qpoGen, ...
                    odeOpts);


            XCurve = ...
                qpo.XCurve;


            if size(XCurve,1) < 3

                error( ...
                    'Corrected QPO XCurve does not contain position states.');

            end


            % Actual corrected cross-section size, independent of the seed K.
            curveCenter = ...
                mean( ...
                    XCurve(1:3,:), ...
                    2);


            dR = ...
                XCurve(1:3,:) - ...
                curveCenter;


            radiusSamplesND = ...
                vecnorm( ...
                    dR, ...
                    2, ...
                    1);


            actualSectionRadiusRMS_km(ip,is) = ...
                sqrt(mean(radiusSamplesND.^2))* ...
                sys.Lstar_km;


            actualSectionRadiusMax_km(ip,is) = ...
                max(radiusSamplesND)* ...
                sys.Lstar_km;


            sectionCentroidOffset_km(ip,is) = ...
                norm( ...
                    curveCenter - ...
                    po.x0(1:3))* ...
                sys.Lstar_km;


            qpoRho(ip,is) = ...
                qpo.rho;


            qpoT(ip,is) = ...
                qpo.T;


            qpoStore{ip,is} = ...
                qpo;


            qpoSuccess(ip,is) = ...
                true;


            fprintf('    corrected RMS section radius = %.6f km\n', ...
                actualSectionRadiusRMS_km(ip,is));

            fprintf('    corrected max section radius = %.6f km\n', ...
                actualSectionRadiusMax_km(ip,is));

            fprintf('    section centroid offset      = %.6f km\n', ...
                sectionCentroidOffset_km(ip,is));


            %% Determine the physically distinct common phases once.

            if isempty(phaseStartIndices)

                phaseStartIndices = ...
                    unique_equal_phase_start_indices( ...
                        qpo.N, ...
                        nCollectors);


                nPhaseCases = ...
                    numel(phaseStartIndices);


                fprintf('\n');
                fprintf('  Unique equal-phase collector triangles = %d\n', ...
                    nPhaseCases);


                %% Allocate case arrays now that nPhaseCases is known.

                dims = ...
                    [nParents,nSizes,nPhaseCases];


                normalRMSDifference_deg = nan(dims);
                normalMaxDifference_deg = nan(dims);

                baselineMedian_km = nan(dims);
                baselineRMS_km = nan(dims);
                baselineMax_km = nan(dims);

                triangleAreaMedian_km2 = nan(dims);
                triangleAspectMedian = nan(dims);

                baselineToTorusRatio = nan(dims);
                areaToTorusRadius2 = nan(dims);

                accessibleFraction = nan(dims);

                rhoMax_Glambda = nan(dims);
                rhoRMS_Glambda = nan(dims);

                uvIsotropy = nan(dims);
                uvAngularFill = nan(dims);
                uvRadialFill = nan(dims);
                uvCellOccupancy = nan(dims);
                uvHullFraction = nan(dims);

                fringeSpacing_mas = nan(dims);
                beamAxisRatioProxy = nan(dims);

                caseSuccess = false(dims);
                caseFailureMessage = cell(dims);

                normalTrackStore = cell(dims);
                timeTrackStore = cell(dims);
                uvPointStore = cell(dims);

            end


        catch ME

            qpoSuccess(ip,is) = ...
                false;


            qpoFailureMessage{ip,is} = ...
                ME.message;


            fprintf('    QPO GENERATION FAILED: %s\n', ...
                ME.message);

        end

    end


    successfulSizes = ...
        find(qpoSuccess(ip,:));


    if isempty(successfulSizes)

        warning( ...
            'No successful QPOs for parent %d. Skipping this parent.', ...
            ip);

        continue

    end


    %% --------------------------------------------------------
    %  4.2 CHOOSE A REFERENCE FORMATION FOR TARGET SELECTION
    %  --------------------------------------------------------
    %
    % Use the successful torus whose ACTUAL corrected size is closest to the
    % median successful radius.  Phase case 1 is used only to define one
    % target for this parent.  The target is then frozen across all cases.

    successfulRadii = ...
        actualSectionRadiusRMS_km(ip,successfulSizes);


    medianRadius = ...
        median(successfulRadii);


    [~,iRefLocal] = ...
        min(abs(successfulRadii-medianRadius));


    isRef = ...
        successfulSizes(iRefLocal);


    referenceSizeIndex(ip) = ...
        isRef;


    iphRef = ...
        1;


    referencePhaseCase(ip) = ...
        iphRef;


    qpoRef = ...
        qpoStore{ip,isRef};


    qpoGenRef = ...
        qpoGenBase;

    qpoGenRef.initialAmplitude = ...
        torusSeedAmplitudes(isRef);


    jCollectorsRef = ...
        qpo_equal_phase_indices( ...
            qpoRef.N, ...
            nCollectors, ...
            phaseStartIndices(iphRef));


    XCollectors0Ref = ...
        qpoRef.XCurve(:,jCollectorsRef);


    fprintf('\n');
    fprintf('  Reference formation for target scan\n');
    fprintf('    size index                = %d\n',isRef);
    fprintf('    seed K                    = %.3e\n', ...
        torusSeedAmplitudes(isRef));
    fprintf('    actual RMS torus radius  = %.6f km\n', ...
        actualSectionRadiusRMS_km(ip,isRef));
    fprintf('    phaseStartIndex          = %d\n', ...
        phaseStartIndices(iphRef));


    geomRef = ...
        compute_three_collector_qpo_geometry( ...
            mu, ...
            xCombiner0, ...
            TPrimary, ...
            XCollectors0Ref, ...
            nPrimaryPeriods, ...
            nSavePerPrimaryPeriod, ...
            odeOpts, ...
            normalReferenceSynodic);


    % add_dimensional_geometry expects the same metadata block used by the
    % validated end-to-end pipeline.  compute_three_collector_qpo_geometry
    % returns only the propagated geometry, so attach the QPO/CR3BP metadata
    % before dimensionalization.
    geomRef.config = ...
        make_geometry_config( ...
            sys, ...
            qpoRef, ...
            qpoGenRef, ...
            poFile, ...
            orbitIndex, ...
            TPrimary, ...
            jCollectorsRef, ...
            phaseStartIndices(iphRef), ...
            lambda_m, ...
            normalReferenceSynodic);


    geomRef = ...
        add_dimensional_geometry( ...
            geomRef, ...
            sys, ...
            timePlotUnit);


    normalRef = ...
        extract_normal_track( ...
            geomRef);


    %% --------------------------------------------------------
    %  4.3 CHOOSE ONE TARGET FOR THIS PARENT
    %  --------------------------------------------------------

    scanTargets = ...
        sample_sphere_targets( ...
            nTargetsPerParentScan, ...
            'fibonacci');


    fprintf('\n');
    fprintf('  Running %d-direction target scan for parent %d...\n', ...
        nTargetsPerParentScan,ip);


    scan = ...
        evaluate_target_suitability_scan( ...
            geomRef, ...
            sys, ...
            scanTargets, ...
            lambda_m, ...
            maxOffNormalDeg, ...
            usePlaneNormalSignAmbiguity, ...
            useParallelTargetScan);


    validScan = ...
        isfinite(scan.score) & ...
        isfinite(scan.raDeg) & ...
        isfinite(scan.decDeg);


    if ~any(validScan)

        warning( ...
            'No finite target-scan candidate for parent %d.', ...
            ip);

        continue

    end


    validIdxScan = ...
        find(validScan);


    [~,iBestLocal] = ...
        max(scan.score(validIdxScan));


    bestIdx = ...
        validIdxScan(iBestLocal);


    targetRA_deg(ip) = ...
        scan.raDeg(bestIdx);


    targetDec_deg(ip) = ...
        scan.decDeg(bestIdx);


    targetScore(ip) = ...
        scan.score(bestIdx);


    targetAccessibleFractionReference(ip) = ...
        scan.accessibleFraction(bestIdx);


    targetUVIsotropyReference(ip) = ...
        scan.uvAnisotropy(bestIdx);


    targetRhoMaxReference_Glambda(ip) = ...
        scan.uvMax_Glambda(bestIdx);


    fprintf('    selected RA               = %.6f deg\n', ...
        targetRA_deg(ip));

    fprintf('    selected Dec              = %.6f deg\n', ...
        targetDec_deg(ip));

    fprintf('    reference target score    = %.6f\n', ...
        targetScore(ip));

    fprintf('    reference f_acc           = %.6f\n', ...
        targetAccessibleFractionReference(ip));

    fprintf('    reference eta_uv          = %.6f\n', ...
        targetUVIsotropyReference(ip));

    fprintf('    reference rho_max         = %.6f Glambda\n', ...
        targetRhoMaxReference_Glambda(ip));


    targetParent = ...
        make_target_definition( ...
            'inertial', ...
            true, ...
            targetRA_deg(ip), ...
            targetDec_deg(ip), ...
            [1 0 0], ...
            [1 0 0], ...
            skyNorthReference);


    %% --------------------------------------------------------
    %  4.4 FLY ALL TORUS SIZE x PHASE CASES
    %  --------------------------------------------------------

    for is = 1:nSizes

        if ~qpoSuccess(ip,is)
            continue
        end


        qpo = ...
            qpoStore{ip,is};


        qpoGenCase = ...
            qpoGenBase;

        qpoGenCase.initialAmplitude = ...
            torusSeedAmplitudes(is);


        for iph = 1:nPhaseCases

            phaseStartIndex = ...
                phaseStartIndices(iph);


            fprintf( ...
                '  parent %d/%d | size %d/%d | phase %d/%d (start %d)\n', ...
                ip,nParents,is,nSizes,iph,nPhaseCases,phaseStartIndex);


            try

                if is == isRef && iph == iphRef

                    geomCase = ...
                        geomRef;

                else

                    jCollectors = ...
                        qpo_equal_phase_indices( ...
                            qpo.N, ...
                            nCollectors, ...
                            phaseStartIndex);


                    XCollectors0 = ...
                        qpo.XCurve(:,jCollectors);


                    geomCase = ...
                        compute_three_collector_qpo_geometry( ...
                            mu, ...
                            xCombiner0, ...
                            TPrimary, ...
                            XCollectors0, ...
                            nPrimaryPeriods, ...
                            nSavePerPrimaryPeriod, ...
                            odeOpts, ...
                            normalReferenceSynodic);


                    geomCase.config = ...
                        make_geometry_config( ...
                            sys, ...
                            qpo, ...
                            qpoGenCase, ...
                            poFile, ...
                            orbitIndex, ...
                            TPrimary, ...
                            jCollectors, ...
                            phaseStartIndex, ...
                            lambda_m, ...
                            normalReferenceSynodic);


                    geomCase = ...
                        add_dimensional_geometry( ...
                            geomCase, ...
                            sys, ...
                            timePlotUnit);

                end


                %% --------------------------------------------
                %  TARGET-INDEPENDENT FORMATION MORPHOLOGY
                %  --------------------------------------------

                normalCase = ...
                    extract_normal_track( ...
                        geomCase);


                [dNRMS,dNMax] = ...
                    normal_track_difference_deg( ...
                        normalRef, ...
                        normalCase);


                normalRMSDifference_deg(ip,is,iph) = ...
                    dNRMS;


                normalMaxDifference_deg(ip,is,iph) = ...
                    dNMax;


                B12_km = ...
                    extract_baseline_track_km( ...
                        geomCase, ...
                        '12', ...
                        sys);


                B13_km = ...
                    extract_baseline_track_km( ...
                        geomCase, ...
                        '13', ...
                        sys);


                B23_km = ...
                    extract_baseline_track_km( ...
                        geomCase, ...
                        '23', ...
                        sys);


                b12 = ...
                    vecnorm(B12_km,2,2);

                b13 = ...
                    vecnorm(B13_km,2,2);

                b23 = ...
                    vecnorm(B23_km,2,2);


                bAll = [ ...
                    b12; ...
                    b13; ...
                    b23];


                baselineMedian_km(ip,is,iph) = ...
                    median(bAll);


                baselineRMS_km(ip,is,iph) = ...
                    sqrt(mean(bAll.^2));


                baselineMax_km(ip,is,iph) = ...
                    max(bAll);


                areaTrack = ...
                    0.5* ...
                    vecnorm( ...
                        cross( ...
                            B12_km, ...
                            B13_km, ...
                            2), ...
                        2, ...
                        2);


                triangleAreaMedian_km2(ip,is,iph) = ...
                    median(areaTrack);


                bMatrix = [ ...
                    b12, ...
                    b13, ...
                    b23];


                bMin = ...
                    min(bMatrix,[],2);


                bMax = ...
                    max(bMatrix,[],2);


                validAspect = ...
                    bMin > 0;


                if any(validAspect)

                    triangleAspectMedian(ip,is,iph) = ...
                        median( ...
                            bMax(validAspect)./ ...
                            bMin(validAspect));

                end


                Rtorus = ...
                    actualSectionRadiusRMS_km(ip,is);


                if isfinite(Rtorus) && Rtorus > 0

                    baselineToTorusRatio(ip,is,iph) = ...
                        baselineMedian_km(ip,is,iph)/ ...
                        Rtorus;


                    areaToTorusRadius2(ip,is,iph) = ...
                        triangleAreaMedian_km2(ip,is,iph)/ ...
                        Rtorus^2;

                end


                %% --------------------------------------------
                %  FIXED-TARGET ACCESSIBILITY + UV MORPHOLOGY
                %  --------------------------------------------

                accessCase = ...
                    compute_triangle_accessibility( ...
                        geomCase, ...
                        targetParent, ...
                        maxOffNormalDeg, ...
                        usePlaneNormalSignAmbiguity);


                uvCase = ...
                    compute_fixed_target_uv( ...
                        geomCase, ...
                        targetParent, ...
                        sys, ...
                        lambda_m);


                accessMask = ...
                    accessCase.accessible(:);


                accessibleFraction(ip,is,iph) = ...
                    nnz(accessMask)/ ...
                    numel(accessMask);


                uvMetrics = ...
                    compute_case_uv_metrics( ...
                        uvCase, ...
                        accessMask, ...
                        nAngularBins, ...
                        nRadialBins, ...
                        nUVCellBins);


                rhoMax_Glambda(ip,is,iph) = ...
                    uvMetrics.rhoMax_Glambda;


                rhoRMS_Glambda(ip,is,iph) = ...
                    uvMetrics.rhoRMS_Glambda;


                uvIsotropy(ip,is,iph) = ...
                    uvMetrics.isotropy;


                uvAngularFill(ip,is,iph) = ...
                    uvMetrics.angularFill;


                uvRadialFill(ip,is,iph) = ...
                    uvMetrics.radialFill;


                uvCellOccupancy(ip,is,iph) = ...
                    uvMetrics.cellOccupancy;


                uvHullFraction(ip,is,iph) = ...
                    uvMetrics.hullFraction;


                fringeSpacing_mas(ip,is,iph) = ...
                    uvMetrics.fringeSpacing_mas;


                beamAxisRatioProxy(ip,is,iph) = ...
                    uvMetrics.beamAxisRatioProxy;


                normalTrackStore{ip,is,iph} = ...
                    normalCase;


                timeTrackStore{ip,is,iph} = ...
                    extract_time_track( ...
                        geomCase);


                uvPointStore{ip,is,iph} = ...
                    uvMetrics.points;


                caseSuccess(ip,is,iph) = ...
                    true;


            catch ME

                caseSuccess(ip,is,iph) = ...
                    false;


                caseFailureMessage{ip,is,iph} = ...
                    ME.message;


                fprintf('    CASE FAILED: %s\n', ...
                    ME.message);

            end

        end

    end


    %% --------------------------------------------------------
    %  4.5 PARENT-SPECIFIC DIAGNOSTIC FIGURES
    %  --------------------------------------------------------

    make_parent_normal_diagnostic_figure( ...
        ip, ...
        representativeOrbitIndices(ip), ...
        torusSeedAmplitudes, ...
        actualSectionRadiusRMS_km(ip,:), ...
        phaseStartIndices, ...
        referenceSizeIndex(ip), ...
        normalTrackStore, ...
        normalRMSDifference_deg, ...
        baselineToTorusRatio, ...
        caseSuccess);


    make_parent_uv_phase_figure( ...
        ip, ...
        representativeOrbitIndices(ip), ...
        referenceSizeIndex(ip), ...
        phaseStartIndices, ...
        uvPointStore, ...
        rhoMax_Glambda, ...
        nPhasesToPlot, ...
        caseSuccess);


    make_parent_uv_metric_heatmaps( ...
        ip, ...
        representativeOrbitIndices(ip), ...
        actualSectionRadiusRMS_km(ip,:), ...
        phaseStartIndices, ...
        accessibleFraction, ...
        rhoMax_Glambda, ...
        uvIsotropy, ...
        uvAngularFill, ...
        uvRadialFill, ...
        uvCellOccupancy);


    if saveFigures

        parentFigDir = ...
            fullfile( ...
                outputDir, ...
                sprintf('Parent_%02d_row_%04d', ...
                    ip,representativeOrbitIndices(ip)));


        if ~exist(parentFigDir,'dir')
            mkdir(parentFigDir);
        end


        figList = ...
            findall(groot,'Type','figure');


        for jf = 1:numel(figList)

            fig = ...
                figList(jf);


            figName = ...
                get(fig,'Name');


            if isempty(figName)
                continue
            end


            if contains(figName,sprintf('Parent %d ',ip))

                safeName = ...
                    regexprep( ...
                        figName, ...
                        '[^a-zA-Z0-9_+-]', ...
                        '_');


                exportgraphics( ...
                    fig, ...
                    fullfile( ...
                        parentFigDir, ...
                        [safeName,'.png']), ...
                    'Resolution',200);

            end

        end

    end

end


%% ============================================================
%  5. BUILD COMPLETE CASE TABLE
%  ============================================================

caseRows = {};

rowCounter = 0;


for ip = 1:nParents

    for is = 1:nSizes

        if isempty(phaseStartIndices)
            continue
        end


        for iph = 1:nPhaseCases

            rowCounter = ...
                rowCounter + 1;


            caseRows(rowCounter,:) = { ...
                ip, ...
                representativeOrbitIndices(ip), ...
                torusSeedAmplitudes(is), ...
                actualSectionRadiusRMS_km(ip,is), ...
                actualSectionRadiusMax_km(ip,is), ...
                sectionCentroidOffset_km(ip,is), ...
                phaseStartIndices(iph), ...
                caseSuccess(ip,is,iph), ...
                normalRMSDifference_deg(ip,is,iph), ...
                normalMaxDifference_deg(ip,is,iph), ...
                baselineMedian_km(ip,is,iph), ...
                baselineRMS_km(ip,is,iph), ...
                baselineMax_km(ip,is,iph), ...
                triangleAreaMedian_km2(ip,is,iph), ...
                triangleAspectMedian(ip,is,iph), ...
                baselineToTorusRatio(ip,is,iph), ...
                areaToTorusRadius2(ip,is,iph), ...
                accessibleFraction(ip,is,iph), ...
                rhoMax_Glambda(ip,is,iph), ...
                rhoRMS_Glambda(ip,is,iph), ...
                uvIsotropy(ip,is,iph), ...
                uvAngularFill(ip,is,iph), ...
                uvRadialFill(ip,is,iph), ...
                uvCellOccupancy(ip,is,iph), ...
                uvHullFraction(ip,is,iph), ...
                fringeSpacing_mas(ip,is,iph), ...
                beamAxisRatioProxy(ip,is,iph)};

        end

    end

end


caseTable = ...
    cell2table( ...
        caseRows, ...
        'VariableNames',{ ...
            'ParentCase', ...
            'CatalogRow', ...
            'SeedAmplitude', ...
            'TorusRadiusRMS_km', ...
            'TorusRadiusMax_km', ...
            'TorusCentroidOffset_km', ...
            'PhaseStartIndex', ...
            'CaseSuccess', ...
            'NormalRMSDifference_deg', ...
            'NormalMaxDifference_deg', ...
            'BaselineMedian_km', ...
            'BaselineRMS_km', ...
            'BaselineMax_km', ...
            'TriangleAreaMedian_km2', ...
            'TriangleAspectMedian', ...
            'BaselineToTorusRatio', ...
            'AreaToTorusRadius2', ...
            'AccessibleFraction', ...
            'RhoMax_Glambda', ...
            'RhoRMS_Glambda', ...
            'UVIsotropy', ...
            'UVAngularFill', ...
            'UVRadialFill', ...
            'UVCellOccupancy', ...
            'UVHullFraction', ...
            'FringeSpacing_mas', ...
            'BeamAxisRatioProxy'});


fprintf('\n');
fprintf('============================================================\n');
fprintf('COMPLETE MORPHOLOGY SURVEY TABLE\n');
fprintf('============================================================\n');

disp(caseTable);


%% ============================================================
%  6. PARENT-LEVEL SUMMARY: WHAT CONTROLS WHAT?
%  ============================================================

parentSummaryRows = {};


for ip = 1:nParents

    isRef = ...
        referenceSizeIndex(ip);


    if ~isfinite(isRef)
        continue
    end


    isRef = ...
        round(isRef);


    iphRef = ...
        referencePhaseCase(ip);


    %% Normal sensitivity to torus size at fixed phase

    sizeNormal = ...
        squeeze( ...
            normalRMSDifference_deg(ip,:,iphRef));


    maxNormalSizeEffect_deg = ...
        max_finite(sizeNormal);


    %% Normal sensitivity to phase at fixed torus size

    phaseNormal = ...
        squeeze( ...
            normalRMSDifference_deg(ip,isRef,:));


    maxNormalPhaseEffect_deg = ...
        max_finite(phaseNormal);


    %% Overall normal sensitivity

    allNormal = ...
        squeeze( ...
            normalRMSDifference_deg(ip,:,:));


    maxNormalOverall_deg = ...
        max_finite(allNormal(:));


    %% Baseline scaling exponent
    %
    % Use phase-median baseline scale for each torus size, then fit:
    %
    %       B_med ~ R_torus^p

    R = ...
        actualSectionRadiusRMS_km(ip,:);


    Bsize = ...
        nan(1,nSizes);


    for is = 1:nSizes

        Bsize(is) = ...
            median_finite( ...
                squeeze( ...
                    baselineMedian_km(ip,is,:)));

    end


    validFit = ...
        isfinite(R) & ...
        R > 0 & ...
        isfinite(Bsize) & ...
        Bsize > 0;


    if nnz(validFit) >= 2

        pFit = ...
            polyfit( ...
                log(R(validFit)), ...
                log(Bsize(validFit)), ...
                1);


        baselineScalingExponent = ...
            pFit(1);

    else

        baselineScalingExponent = ...
            NaN;

    end


    %% Phase variability of interferometric morphology at reference size

    isoPhase = ...
        squeeze( ...
            uvIsotropy(ip,isRef,:));


    rhoPhase = ...
        squeeze( ...
            rhoMax_Glambda(ip,isRef,:));


    cellPhase = ...
        squeeze( ...
            uvCellOccupancy(ip,isRef,:));


    phaseCV_UVIsotropy = ...
        coefficient_of_variation_finite( ...
            isoPhase);


    phaseCV_RhoMax = ...
        coefficient_of_variation_finite( ...
            rhoPhase);


    phaseCV_UVCellOccupancy = ...
        coefficient_of_variation_finite( ...
            cellPhase);


    parentSummaryRows(end+1,:) = { ...
        ip, ...
        representativeOrbitIndices(ip), ...
        abs(representativeZ0(ip)), ...
        targetRA_deg(ip), ...
        targetDec_deg(ip), ...
        targetScore(ip), ...
        maxNormalSizeEffect_deg, ...
        maxNormalPhaseEffect_deg, ...
        maxNormalOverall_deg, ...
        baselineScalingExponent, ...
        phaseCV_RhoMax, ...
        phaseCV_UVIsotropy, ...
        phaseCV_UVCellOccupancy};

end


parentSummaryTable = ...
    cell2table( ...
        parentSummaryRows, ...
        'VariableNames',{ ...
            'ParentCase', ...
            'CatalogRow', ...
            'AbsZ0', ...
            'TargetRA_deg', ...
            'TargetDec_deg', ...
            'TargetScore', ...
            'MaxNormalSizeEffect_deg', ...
            'MaxNormalPhaseEffect_deg', ...
            'MaxNormalOverall_deg', ...
            'BaselineScalingExponent', ...
            'PhaseCV_RhoMax', ...
            'PhaseCV_UVIsotropy', ...
            'PhaseCV_UVCellOccupancy'});


fprintf('\n');
fprintf('============================================================\n');
fprintf('PARENT-LEVEL SCIENTIFIC SUMMARY\n');
fprintf('============================================================\n');

disp(parentSummaryTable);


%% ============================================================
%  7. GLOBAL BASELINE-SCALING FIGURE
%  ============================================================

figScaling = ...
    figure( ...
        'Name','Global baseline scaling', ...
        'Color','w', ...
        'Position',[100 100 1350 560]);


tl = ...
    tiledlayout( ...
        figScaling, ...
        1,2, ...
        'TileSpacing','compact', ...
        'Padding','compact');


%% Baseline scale versus actual torus radius

ax = ...
    nexttile(tl);


hold(ax,'on');
grid(ax,'on');


parentColors = ...
    lines(nParents);


for ip = 1:nParents

    R = ...
        actualSectionRadiusRMS_km(ip,:);


    Bsize = ...
        nan(1,nSizes);


    for is = 1:nSizes

        Bsize(is) = ...
            median_finite( ...
                squeeze( ...
                    baselineMedian_km(ip,is,:)));

    end


    valid = ...
        isfinite(R) & ...
        R > 0 & ...
        isfinite(Bsize) & ...
        Bsize > 0;


    if any(valid)

        loglog( ...
            ax, ...
            R(valid), ...
            Bsize(valid), ...
            'o-', ...
            'Color',parentColors(ip,:), ...
            'MarkerFaceColor',parentColors(ip,:), ...
            'DisplayName',sprintf('Parent %d',ip));

    end

end


xlabel(ax,'Corrected QPO section RMS radius [km]');
ylabel(ax,'Phase-median physical baseline [km]');

title(ax,'Baseline scale versus actual corrected torus size');

legend(ax,'Location','best');


%% Baseline / torus-radius ratio

ax = ...
    nexttile(tl);


hold(ax,'on');
grid(ax,'on');


for ip = 1:nParents

    R = ...
        actualSectionRadiusRMS_km(ip,:);


    ratioSize = ...
        nan(1,nSizes);


    for is = 1:nSizes

        ratioSize(is) = ...
            median_finite( ...
                squeeze( ...
                    baselineToTorusRatio(ip,is,:)));

    end


    valid = ...
        isfinite(R) & ...
        R > 0 & ...
        isfinite(ratioSize);


    if any(valid)

        semilogx( ...
            ax, ...
            R(valid), ...
            ratioSize(valid), ...
            'o-', ...
            'Color',parentColors(ip,:), ...
            'MarkerFaceColor',parentColors(ip,:), ...
            'DisplayName',sprintf('Parent %d',ip));

    end

end


xlabel(ax,'Corrected QPO section RMS radius [km]');
ylabel('$\mathrm{median}(B)/R_{\rm torus}$');

title(ax,'Does torus size only rescale the formation?');

legend(ax,'Location','best');


title( ...
    tl, ...
    'QPO torus size and physical interferometer scale');


%% ============================================================
%  8. GLOBAL NORMAL-SENSITIVITY FIGURE
%  ============================================================

figNormalSummary = ...
    figure( ...
        'Name','Global normal sensitivity', ...
        'Color','w', ...
        'Position',[120 120 1200 520]);


tl = ...
    tiledlayout( ...
        figNormalSummary, ...
        1,2, ...
        'TileSpacing','compact', ...
        'Padding','compact');


ax = ...
    nexttile(tl);


bar( ...
    ax, ...
    parentSummaryTable.ParentCase, ...
    [ ...
        parentSummaryTable.MaxNormalSizeEffect_deg, ...
        parentSummaryTable.MaxNormalPhaseEffect_deg]);


grid(ax,'on');

xlabel(ax,'Representative halo parent');
ylabel('Maximum RMS normal-track difference [deg]');

legend( ...
    ax, ...
    {'vary torus size','vary starting phase'}, ...
    'Location','best');

title(ax,'Sensitivity of array-plane normal pattern');


ax = ...
    nexttile(tl);


bar( ...
    ax, ...
    parentSummaryTable.ParentCase, ...
    parentSummaryTable.BaselineScalingExponent);


grid(ax,'on');

xlabel(ax,'Representative halo parent');
ylabel('Power-law exponent $p$');

yline(ax,1,'k--','$p=1$');

title(ax,'Fit of $B_{\rm med}\propto R_{\rm torus}^{p}$');


title( ...
    tl, ...
    'Formation-level hypotheses: normal invariance and baseline scaling');


%% ============================================================
%  9. SAVE RESULTS
%  ============================================================

if saveResults

    matFile = ...
        fullfile( ...
            outputDir, ...
            'QPO_Halo_Morphology_Survey.mat');


    save( ...
        matFile, ...
        'systemName', ...
        'sys', ...
        'poFile', ...
        'parentSelectionTable', ...
        'representativeOrbitIndices', ...
        'torusSeedAmplitudes', ...
        'phaseStartIndices', ...
        'qpoSuccess', ...
        'qpoFailureMessage', ...
        'actualSectionRadiusRMS_km', ...
        'actualSectionRadiusMax_km', ...
        'sectionCentroidOffset_km', ...
        'qpoRho', ...
        'qpoT', ...
        'referenceSizeIndex', ...
        'targetRA_deg', ...
        'targetDec_deg', ...
        'targetScore', ...
        'targetAccessibleFractionReference', ...
        'targetUVIsotropyReference', ...
        'targetRhoMaxReference_Glambda', ...
        'caseSuccess', ...
        'caseFailureMessage', ...
        'normalRMSDifference_deg', ...
        'normalMaxDifference_deg', ...
        'baselineMedian_km', ...
        'baselineRMS_km', ...
        'baselineMax_km', ...
        'triangleAreaMedian_km2', ...
        'triangleAspectMedian', ...
        'baselineToTorusRatio', ...
        'areaToTorusRadius2', ...
        'accessibleFraction', ...
        'rhoMax_Glambda', ...
        'rhoRMS_Glambda', ...
        'uvIsotropy', ...
        'uvAngularFill', ...
        'uvRadialFill', ...
        'uvCellOccupancy', ...
        'uvHullFraction', ...
        'fringeSpacing_mas', ...
        'beamAxisRatioProxy', ...
        'caseTable', ...
        'parentSummaryTable', ...
        '-v7.3');


    writetable( ...
        caseTable, ...
        fullfile( ...
            outputDir, ...
            'QPO_Halo_Morphology_Cases.csv'));


    writetable( ...
        parentSummaryTable, ...
        fullfile( ...
            outputDir, ...
            'QPO_Halo_Morphology_ParentSummary.csv'));


    fprintf('\nSaved survey results to:\n  %s\n',matFile);

end


%% ============================================================
%  10. INTERPRETATION GUIDE
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('HOW TO INTERPRET THE SURVEY\n');
fprintf('============================================================\n');

fprintf(['1. If MaxNormalSizeEffect_deg is tiny while baseline scale changes\n', ...
         '   strongly, torus size mainly rescales the interferometer.\n']);

fprintf(['2. If MaxNormalPhaseEffect_deg is tiny but UV metrics vary with phase,\n', ...
         '   starting phase changes in-plane baselines without changing the\n', ...
         '   array-plane pointing pattern substantially.\n']);

fprintf(['3. A BaselineScalingExponent near 1 supports approximately linear\n', ...
         '   baseline scaling with corrected torus radius.\n']);

fprintf(['4. UVIsotropy near 1, large UVAngularFill / UVRadialFill, and large\n', ...
         '   UVCellOccupancy indicate better two-dimensional Fourier diversity.\n']);

fprintf(['5. BeamAxisRatioProxy near 1 is favorable; large values indicate a\n', ...
         '   strongly anisotropic synthesized-resolution proxy.\n']);

fprintf(['6. FringeSpacing_mas reports the angular scale 1/rho_max, not a fitted\n', ...
         '   CLEAN beam FWHM.  Use the end-to-end reconstruction pipeline for\n', ...
         '   final imaging claims.\n']);

fprintf('============================================================\n');


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================


function config = make_geometry_config( ...
    sys, ...
    qpo, ...
    qpoGen, ...
    poFile, ...
    orbitIndex, ...
    TPrimary, ...
    jCollectors, ...
    phaseStartIndex, ...
    lambda_m, ...
    normalReferenceSynodic)
%MAKE_GEOMETRY_CONFIG
% Attach the metadata expected by add_dimensional_geometry.
%
% This mirrors the metadata block used by the validated end-to-end FINESST
% QPO pipeline.  Keeping it in one helper prevents the reference formation
% and the hundreds of survey cases from drifting apart.

    config = struct();

    config.systemName = sys.name;
    config.mu = sys.mu;
    config.N = qpo.N;

    config.poFile = poFile;
    config.orbitIndex = orbitIndex;

    config.TPrimary = TPrimary;
    config.TQPO = qpo.T;
    config.rhoQPO = qpo.rho;

    config.qpoSource = 'generate';

    % Store the actual seed amplitude used for this corrected torus when it
    % is available.  The rest of the QPO settings come from qpoGenBase.
    config.qpoGen = qpoGen;

    config.jCollectors = jCollectors;
    config.phaseStartIndex = phaseStartIndex;

    config.lambda_m = lambda_m;
    config.normalReferenceSynodic = normalReferenceSynodic;

end


function phaseStartIndices = unique_equal_phase_start_indices(N,nCollectors)
%UNIQUE_EQUAL_PHASE_START_INDICES
% Return only physically distinct equal-phase collector triangles.
%
% Collector labels are ignored when identifying duplicate triangles.  For
% N = 51 and three equally phased collectors this should normally return
% 17 distinct common starting phases.

    canonicalSets = zeros(N,nCollectors);

    for k = 1:N

        j = qpo_equal_phase_indices(N,nCollectors,k);

        canonicalSets(k,:) = sort(j(:).');

    end


    [~,ia] = unique( ...
        canonicalSets, ...
        'rows', ...
        'stable');


    phaseStartIndices = ia(:).';

end


function nHat = extract_normal_track(geom)
%EXTRACT_NORMAL_TRACK
% Robustly identify an N-by-3 unit-normal history inside the FINESST geom
% structure.  Exact common names are checked first, followed by a scored
% recursive search.

    exactPaths = { ...
        'normalHatSynodic', ...
        'nHatSynodic', ...
        'normalSynodic', ...
        'normalHat', ...
        'nHat', ...
        'normal', ...
        'geometry.normalHatSynodic', ...
        'geometry.nHatSynodic', ...
        'formation.normalHatSynodic'};


    for k = 1:numel(exactPaths)

        [ok,value] = try_get_nested_field(geom,exactPaths{k});

        if ok && is_numeric_xyz_history(value)

            nHat = orient_xyz_history(value);
            nHat = normalize_rows_safe(nHat);
            nHat = enforce_row_sign_continuity(nHat);
            return

        end

    end


    candidates = collect_numeric_xyz_fields(geom,'geom');


    bestScore = -Inf;
    bestValue = [];


    for k = 1:numel(candidates)

        pathLower = lower(candidates(k).path);

        if contains(pathLower,'reference') || ...
           contains(pathLower,'config')
            continue
        end


        if ~(contains(pathLower,'normal') || contains(pathLower,'nhat'))
            continue
        end


        X = orient_xyz_history(candidates(k).value);

        rowNorm = vecnorm(X,2,2);

        finiteNorm = rowNorm(isfinite(rowNorm));

        if isempty(finiteNorm)
            continue
        end


        score = 0;

        if contains(pathLower,'normal')
            score = score + 8;
        end

        if contains(pathLower,'nhat')
            score = score + 8;
        end

        if contains(pathLower,'synodic')
            score = score + 3;
        end

        if contains(pathLower,'unit') || contains(pathLower,'hat')
            score = score + 2;
        end


        medianNorm = median(finiteNorm);

        if abs(medianNorm-1) < 0.1
            score = score + 8;
        elseif abs(medianNorm-1) < 0.5
            score = score + 3;
        end


        if score > bestScore

            bestScore = score;
            bestValue = X;

        end

    end


    if isempty(bestValue)

        error( ...
            ['Could not identify the array-normal history in geom. ', ...
             'Inspect fieldnames(geom) and add the current field name to ', ...
             'extract_normal_track.']);

    end


    nHat = normalize_rows_safe(bestValue);
    nHat = enforce_row_sign_continuity(nHat);

end


function B_km = extract_baseline_track_km(geom,pairLabel,sys)
%EXTRACT_BASELINE_TRACK_KM
% Locate baseline B12/B13/B23 in geom.  Dimensional km fields are preferred.
% If only a dimensionless CR3BP baseline is found, multiply by L*.

    pairLabel = char(pairLabel);


    exactPaths = { ...
        ['B',pairLabel,'_km'], ...
        ['b',pairLabel,'_km'], ...
        ['baseline',pairLabel,'_km'], ...
        ['B',pairLabel], ...
        ['b',pairLabel], ...
        ['baseline',pairLabel], ...
        ['geometry.B',pairLabel,'_km'], ...
        ['geometry.B',pairLabel], ...
        ['formation.B',pairLabel,'_km'], ...
        ['formation.B',pairLabel]};


    for k = 1:numel(exactPaths)

        [ok,value] = try_get_nested_field(geom,exactPaths{k});

        if ok && is_numeric_xyz_history(value)

            B = orient_xyz_history(value);

            if contains(lower(exactPaths{k}),'km')

                B_km = B;

            else

                B_km = B*sys.Lstar_km;

            end

            return

        end

    end


    candidates = collect_numeric_xyz_fields(geom,'geom');


    bestScore = -Inf;
    bestValue = [];
    bestPath = '';


    for k = 1:numel(candidates)

        p = lower(candidates(k).path);


        if contains(p,'uvw') || contains(p,'velocity')
            continue
        end


        pairHit = ...
            contains(p,['b',pairLabel]) || ...
            contains(p,['baseline',pairLabel]) || ...
            contains(p,['baseline_',pairLabel]);


        if ~pairHit
            continue
        end


        score = 0;

        if contains(p,'baseline')
            score = score + 7;
        end

        if contains(p,['b',pairLabel])
            score = score + 7;
        end

        if contains(p,'km')
            score = score + 12;
        end

        if contains(p,'synodic')
            score = score + 2;
        end


        if score > bestScore

            bestScore = score;
            bestValue = orient_xyz_history(candidates(k).value);
            bestPath = candidates(k).path;

        end

    end


    if isempty(bestValue)

        error( ...
            ['Could not identify baseline ',pairLabel,' in geom. ', ...
             'Inspect fieldnames(geom) and add its current name to ', ...
             'extract_baseline_track_km.']);

    end


    if contains(lower(bestPath),'km')

        B_km = bestValue;

    else

        B_km = bestValue*sys.Lstar_km;

    end

end


function t = extract_time_track(geom)
%EXTRACT_TIME_TRACK

    if isfield(geom,'t')

        t = geom.t(:);
        return

    end


    if isfield(geom,'tau')

        t = geom.tau(:);
        return

    end


    error('Could not identify geom time vector.');

end


function [rmsDeg,maxDeg] = normal_track_difference_deg(nRef,nCase)
%NORMAL_TRACK_DIFFERENCE_DEG
% Compare plane normals with n == -n ambiguity.

    n = min(size(nRef,1),size(nCase,1));

    if n < 1
        rmsDeg = NaN;
        maxDeg = NaN;
        return
    end


    a = nRef(1:n,:);
    b = nCase(1:n,:);


    d = ...
        sum(a.*b,2);


    d = ...
        min( ...
            1, ...
            max( ...
                -1, ...
                abs(d)));


    angDeg = ...
        acosd(d);


    rmsDeg = ...
        sqrt(mean(angDeg.^2,'omitnan'));


    maxDeg = ...
        max(angDeg,[],'omitnan');

end


function metrics = compute_case_uv_metrics( ...
    uvCase, ...
    accessMask, ...
    nAngularBins, ...
    nRadialBins, ...
    nUVCellBins)
%COMPUTE_CASE_UV_METRICS
% Compute target-dependent Fourier-plane morphology metrics using physical
% baseline samples only.  Hermitian counterparts are added only where a
% symmetric Fourier-plane geometric metric requires them.

    baselineFields = { ...
        'uvw12', ...
        'uvw13', ...
        'uvw23'};


    u = [];
    v = [];
    baselineID = [];


    for b = 1:3

        uvw = ...
            uvCase.(baselineFields{b});


        u = [ ...
            u; ...
            uvw(accessMask,1)/1e9];


        v = [ ...
            v; ...
            uvw(accessMask,2)/1e9];


        baselineID = [ ...
            baselineID; ...
            b*ones(nnz(accessMask),1)];

    end


    rho = ...
        hypot(u,v);


    metrics.points = struct( ...
        'u_Glambda',u, ...
        'v_Glambda',v, ...
        'baselineID',baselineID);


    if isempty(rho) || all(rho == 0)

        metrics.rhoMax_Glambda = 0;
        metrics.rhoRMS_Glambda = 0;
        metrics.isotropy = 0;
        metrics.angularFill = 0;
        metrics.radialFill = 0;
        metrics.cellOccupancy = 0;
        metrics.hullFraction = 0;
        metrics.fringeSpacing_mas = Inf;
        metrics.beamAxisRatioProxy = Inf;
        return

    end


    metrics.rhoMax_Glambda = ...
        max(rho);


    metrics.rhoRMS_Glambda = ...
        sqrt(mean(rho.^2));


    %% --------------------------------------------------------
    %  UV covariance isotropy
    %  --------------------------------------------------------
    %
    % For eigenvalues lambda_max >= lambda_min:
    %
    %   isotropy = 2 sqrt(lambda_max lambda_min)/(lambda_max+lambda_min)
    %
    % equals 1 for an isotropic covariance and tends to zero for a line.

    U = [ ...
        u; ...
       -u];


    V = [ ...
        v; ...
       -v];


    C = ...
        cov([U,V],1);


    eigC = ...
        sort(real(eig(C)),'descend');


    eigC = ...
        max(eigC,0);


    if numel(eigC) < 2 || sum(eigC) <= 0

        metrics.isotropy = 0;
        metrics.beamAxisRatioProxy = Inf;

    else

        metrics.isotropy = ...
            2*sqrt(eigC(1)*eigC(2))/ ...
            (eigC(1)+eigC(2));


        if eigC(2) > 0

            % Fourier elongation maps inversely into image-plane resolution.
            metrics.beamAxisRatioProxy = ...
                sqrt(eigC(1)/eigC(2));

        else

            metrics.beamAxisRatioProxy = ...
                Inf;

        end

    end


    %% --------------------------------------------------------
    %  Angular fill modulo pi
    %  --------------------------------------------------------

    beta = ...
        mod( ...
            atan2(v,u), ...
            pi);


    betaEdges = ...
        linspace( ...
            0, ...
            pi, ...
            nAngularBins+1);


    betaCounts = ...
        histcounts( ...
            beta, ...
            betaEdges);


    metrics.angularFill = ...
        nnz(betaCounts > 0)/ ...
        nAngularBins;


    %% --------------------------------------------------------
    %  Radial fill
    %  --------------------------------------------------------

    rhoEdges = ...
        linspace( ...
            0, ...
            metrics.rhoMax_Glambda, ...
            nRadialBins+1);


    rhoCounts = ...
        histcounts( ...
            rho, ...
            rhoEdges);


    metrics.radialFill = ...
        nnz(rhoCounts > 0)/ ...
        nRadialBins;


    %% --------------------------------------------------------
    %  Normalized UV-cell occupancy
    %  --------------------------------------------------------
    %
    % Divide u and v by this case's rho_max.  This removes pure scale and
    % measures normalized Fourier morphology.

    qU = ...
        U/metrics.rhoMax_Glambda;


    qV = ...
        V/metrics.rhoMax_Glambda;


    cellEdges = ...
        linspace( ...
            -1, ...
            1, ...
            nUVCellBins+1);


    cellCounts = ...
        histcounts2( ...
            qU, ...
            qV, ...
            cellEdges, ...
            cellEdges);


    metrics.cellOccupancy = ...
        nnz(cellCounts > 0)/ ...
        (nUVCellBins^2);


    %% --------------------------------------------------------
    %  Normalized convex-hull area
    %  --------------------------------------------------------
    %
    % Normalize by the area pi of the unit disk after rho_max scaling.

    pointsUnique = ...
        unique( ...
            [qU,qV], ...
            'rows');


    if size(pointsUnique,1) >= 3

        try

            iHull = ...
                convhull( ...
                    pointsUnique(:,1), ...
                    pointsUnique(:,2));


            areaHull = ...
                polyarea( ...
                    pointsUnique(iHull,1), ...
                    pointsUnique(iHull,2));


            metrics.hullFraction = ...
                min( ...
                    1, ...
                    areaHull/pi);

        catch

            metrics.hullFraction = ...
                0;

        end

    else

        metrics.hullFraction = ...
            0;

    end


    %% --------------------------------------------------------
    %  Angular scale associated with rho_max
    %  --------------------------------------------------------
    %
    % This is fringe spacing 1/rho_max, not CLEAN-beam FWHM.

    masPerRad = ...
        (180/pi)*3600*1000;


    metrics.fringeSpacing_mas = ...
        masPerRad/ ...
        (metrics.rhoMax_Glambda*1e9);

end


function make_parent_normal_diagnostic_figure( ...
    ip, ...
    catalogRow, ...
    torusSeeds, ...
    torusRadius_km, ...
    phaseStartIndices, ...
    refSizeIndex, ...
    normalTrackStore, ...
    normalRMSDifference_deg, ...
    baselineToTorusRatio, ...
    caseSuccess)
%MAKE_PARENT_NORMAL_DIAGNOSTIC_FIGURE

    if ~isfinite(refSizeIndex)
        return
    end


    refSizeIndex = ...
        round(refSizeIndex);


    nSizes = ...
        numel(torusSeeds);


    nPhases = ...
        numel(phaseStartIndices);


    fig = ...
        figure( ...
            'Name',sprintf('Parent %d normal morphology',ip), ...
            'Color','w', ...
            'Position',[80 80 1450 900]);


    tl = ...
        tiledlayout( ...
            fig, ...
            2,2, ...
            'TileSpacing','compact', ...
            'Padding','compact');


    sizeColors = ...
        lines(nSizes);


    %% Size dependence at phase 1

    ax = ...
        nexttile(tl);


    hold(ax,'on');
    grid(ax,'on');
    axis(ax,'equal');


    draw_unit_sphere_reference(ax);


    for is = 1:nSizes

        if ~caseSuccess(ip,is,1)
            continue
        end


        nHat = ...
            normalTrackStore{ip,is,1};


        plot3( ...
            ax, ...
            nHat(:,1), ...
            nHat(:,2), ...
            nHat(:,3), ...
            'Color',sizeColors(is,:), ...
            'DisplayName',sprintf( ...
                '$R_{\\rm rms}=%.3g$ km', ...
                torusRadius_km(is)));

    end


    xlabel(ax,'$n_x$');
    ylabel(ax,'$n_y$');
    zlabel(ax,'$n_z$');

    view(ax,35,25);

    title(ax,'Normal tracks: varying torus size, fixed phase');

    legend(ax,'Location','best');


    %% Phase dependence at reference size

    ax = ...
        nexttile(tl);


    hold(ax,'on');
    grid(ax,'on');
    axis(ax,'equal');


    draw_unit_sphere_reference(ax);


    phaseColors = ...
        turbo(max(nPhases,2));


    for iph = 1:nPhases

        if ~caseSuccess(ip,refSizeIndex,iph)
            continue
        end


        nHat = ...
            normalTrackStore{ip,refSizeIndex,iph};


        plot3( ...
            ax, ...
            nHat(:,1), ...
            nHat(:,2), ...
            nHat(:,3), ...
            'Color',phaseColors(iph,:), ...
            'HandleVisibility','off');

    end


    xlabel(ax,'$n_x$');
    ylabel(ax,'$n_y$');
    zlabel(ax,'$n_z$');

    view(ax,35,25);

    title(ax,sprintf( ...
        'Normal tracks: all phases at $R_{\\rm rms}=%.3g$ km', ...
        torusRadius_km(refSizeIndex)));


    %% Normal-difference heat map

    ax = ...
        nexttile(tl);


    imagesc( ...
        ax, ...
        phaseStartIndices, ...
        1:nSizes, ...
        squeeze(normalRMSDifference_deg(ip,:,:)));


    set(ax,'YDir','normal');

    xlabel(ax,'Phase start index');
    ylabel(ax,'Torus-size case');

    yticks(ax,1:nSizes);

    yticklabels( ...
        ax, ...
        compose('%.3g km',torusRadius_km));


    title(ax,'RMS normal difference from reference [deg]');

    colorbar(ax);


    %% Baseline / torus radius heat map

    ax = ...
        nexttile(tl);


    imagesc( ...
        ax, ...
        phaseStartIndices, ...
        1:nSizes, ...
        squeeze(baselineToTorusRatio(ip,:,:)));


    set(ax,'YDir','normal');

    xlabel(ax,'Phase start index');
    ylabel(ax,'Torus-size case');

    yticks(ax,1:nSizes);

    yticklabels( ...
        ax, ...
        compose('%.3g km',torusRadius_km));


    title(ax,'$\mathrm{median}(B)/R_{\rm torus}$');

    colorbar(ax);


    title( ...
        tl, ...
        sprintf( ...
            'Parent %d: JPL row %d formation morphology', ...
            ip,catalogRow));

end


function make_parent_uv_phase_figure( ...
    ip, ...
    catalogRow, ...
    refSizeIndex, ...
    phaseStartIndices, ...
    uvPointStore, ...
    rhoMax_Glambda, ...
    nPhasesToPlot, ...
    caseSuccess)
%MAKE_PARENT_UV_PHASE_FIGURE
% Compare actual and rho-normalized UV morphology for representative phases
% at the reference torus size.

    if ~isfinite(refSizeIndex)
        return
    end


    refSizeIndex = ...
        round(refSizeIndex);


    nPhases = ...
        numel(phaseStartIndices);


    phaseCasesToPlot = ...
        unique( ...
            round( ...
                linspace( ...
                    1, ...
                    nPhases, ...
                    min(nPhasesToPlot,nPhases))), ...
            'stable');


    nPlot = ...
        numel(phaseCasesToPlot);


    if nPlot < 1
        return
    end


    rhoPlotMax = ...
        max_finite( ...
            squeeze( ...
                rhoMax_Glambda(ip,refSizeIndex,phaseCasesToPlot)));


    if ~isfinite(rhoPlotMax) || rhoPlotMax <= 0
        return
    end


    fig = ...
        figure( ...
            'Name',sprintf('Parent %d UV phase morphology',ip), ...
            'Color','w', ...
            'Position',[60 80 1700 760]);


    tl = ...
        tiledlayout( ...
            fig, ...
            2,nPlot, ...
            'TileSpacing','compact', ...
            'Padding','compact');


    baselineColors = ...
        lines(3);


    for j = 1:nPlot

        iph = ...
            phaseCasesToPlot(j);


        %% Physical UV scale

        ax = ...
            nexttile(tl,j);


        hold(ax,'on');
        grid(ax,'on');


        if caseSuccess(ip,refSizeIndex,iph)

            P = ...
                uvPointStore{ip,refSizeIndex,iph};


            for b = 1:3

                mask = ...
                    P.baselineID == b;


                plot( ...
                    ax, ...
                    P.u_Glambda(mask), ...
                    P.v_Glambda(mask), ...
                    '.', ...
                    'Color',baselineColors(b,:), ...
                    'MarkerSize',6);


                plot( ...
                    ax, ...
                   -P.u_Glambda(mask), ...
                   -P.v_Glambda(mask), ...
                    '.', ...
                    'Color',baselineColors(b,:), ...
                    'MarkerSize',6, ...
                    'HandleVisibility','off');

            end

        end


        axis(ax,'equal');

        xlim(ax,1.05*[-rhoPlotMax rhoPlotMax]);
        ylim(ax,1.05*[-rhoPlotMax rhoPlotMax]);

        xlabel(ax,'$u$ [G$\lambda$]');
        ylabel(ax,'$v$ [G$\lambda$]');

        title(ax,sprintf( ...
            'start %d: physical scale', ...
            phaseStartIndices(iph)));


        %% Normalized UV morphology

        ax = ...
            nexttile(tl,nPlot+j);


        hold(ax,'on');
        grid(ax,'on');


        if caseSuccess(ip,refSizeIndex,iph)

            P = ...
                uvPointStore{ip,refSizeIndex,iph};


            rhoMax = ...
                rhoMax_Glambda(ip,refSizeIndex,iph);


            if isfinite(rhoMax) && rhoMax > 0

                for b = 1:3

                    mask = ...
                        P.baselineID == b;


                    plot( ...
                        ax, ...
                        P.u_Glambda(mask)/rhoMax, ...
                        P.v_Glambda(mask)/rhoMax, ...
                        '.', ...
                        'Color',baselineColors(b,:), ...
                        'MarkerSize',6);


                    plot( ...
                        ax, ...
                       -P.u_Glambda(mask)/rhoMax, ...
                       -P.v_Glambda(mask)/rhoMax, ...
                        '.', ...
                        'Color',baselineColors(b,:), ...
                        'MarkerSize',6, ...
                        'HandleVisibility','off');

                end

            end

        end


        axis(ax,'equal');

        xlim(ax,[-1.05 1.05]);
        ylim(ax,[-1.05 1.05]);

        xlabel(ax,'$u/\rho_{\max}$');
        ylabel(ax,'$v/\rho_{\max}$');

        title(ax,sprintf( ...
            'start %d: normalized morphology', ...
            phaseStartIndices(iph)));

    end


    title( ...
        tl, ...
        sprintf( ...
            'Parent %d: JPL row %d --- effect of common starting phase on UV coverage', ...
            ip,catalogRow));

end


function make_parent_uv_metric_heatmaps( ...
    ip, ...
    catalogRow, ...
    torusRadius_km, ...
    phaseStartIndices, ...
    accessibleFraction, ...
    rhoMax_Glambda, ...
    uvIsotropy, ...
    uvAngularFill, ...
    uvRadialFill, ...
    uvCellOccupancy)
%MAKE_PARENT_UV_METRIC_HEATMAPS

    fig = ...
        figure( ...
            'Name',sprintf('Parent %d UV metric heatmaps',ip), ...
            'Color','w', ...
            'Position',[70 70 1550 850]);


    tl = ...
        tiledlayout( ...
            fig, ...
            2,3, ...
            'TileSpacing','compact', ...
            'Padding','compact');


    dataList = { ...
        squeeze(accessibleFraction(ip,:,:)), ...
        squeeze(rhoMax_Glambda(ip,:,:)), ...
        squeeze(uvIsotropy(ip,:,:)), ...
        squeeze(uvAngularFill(ip,:,:)), ...
        squeeze(uvRadialFill(ip,:,:)), ...
        squeeze(uvCellOccupancy(ip,:,:))};


    titleList = { ...
        'Accessible fraction', ...
        '$\rho_{\max}$ [G$\lambda$]', ...
        'UV covariance isotropy', ...
        'Angular fill', ...
        'Radial fill', ...
        'Normalized UV-cell occupancy'};


    for k = 1:6

        ax = ...
            nexttile(tl);


        imagesc( ...
            ax, ...
            phaseStartIndices, ...
            1:numel(torusRadius_km), ...
            dataList{k});


        set(ax,'YDir','normal');

        xlabel(ax,'Phase start index');
        ylabel(ax,'Torus-size case');

        yticks(ax,1:numel(torusRadius_km));

        yticklabels( ...
            ax, ...
            compose('%.3g km',torusRadius_km));


        title(ax,titleList{k});

        colorbar(ax);

    end


    title( ...
        tl, ...
        sprintf( ...
            'Parent %d: JPL row %d --- fixed-target interferometric morphology', ...
            ip,catalogRow));

end


function draw_unit_sphere_reference(ax)
%DRAW_UNIT_SPHERE_REFERENCE

    [xs,ys,zs] = ...
        sphere(24);


    surf( ...
        ax, ...
        xs,ys,zs, ...
        'FaceColor',[0.8 0.8 0.8], ...
        'FaceAlpha',0.06, ...
        'EdgeColor',[0.75 0.75 0.75], ...
        'EdgeAlpha',0.18, ...
        'HandleVisibility','off');

end


function [ok,value] = try_get_nested_field(S,pathString)
%TRY_GET_NESTED_FIELD

    parts = ...
        strsplit(pathString,'.');


    value = ...
        S;


    ok = ...
        true;


    for k = 1:numel(parts)

        if isstruct(value) && isfield(value,parts{k})

            value = ...
                value.(parts{k});

        else

            ok = ...
                false;

            value = ...
                [];

            return

        end

    end

end


function tf = is_numeric_xyz_history(X)
%IS_NUMERIC_XYZ_HISTORY

    tf = ...
        isnumeric(X) && ...
        ismatrix(X) && ...
        ~isempty(X) && ...
        (size(X,2) == 3 || size(X,1) == 3);

end


function X = orient_xyz_history(X)
%ORIENT_XYZ_HISTORY

    if size(X,2) == 3

        return

    elseif size(X,1) == 3

        X = X.';

    else

        error('Expected N-by-3 or 3-by-N array.');

    end

end


function X = normalize_rows_safe(X)
%NORMALIZE_ROWS_SAFE

    n = ...
        vecnorm(X,2,2);


    valid = ...
        n > 0 & ...
        isfinite(n);


    X(valid,:) = ...
        X(valid,:)./ ...
        n(valid);


    X(~valid,:) = ...
        NaN;

end


function X = enforce_row_sign_continuity(X)
%ENFORCE_ROW_SIGN_CONTINUITY

    for k = 2:size(X,1)

        if all(isfinite(X(k-1,:))) && ...
           all(isfinite(X(k,:))) && ...
           dot(X(k-1,:),X(k,:)) < 0

            X(k,:) = ...
                -X(k,:);

        end

    end

end


function candidates = collect_numeric_xyz_fields(S,rootName)
%COLLECT_NUMERIC_XYZ_FIELDS
% Recursively collect numeric N-by-3 / 3-by-N arrays.

    candidates = struct( ...
        'path',{}, ...
        'value',{});


    if ~isstruct(S)
        return
    end


    f = ...
        fieldnames(S);


    for k = 1:numel(f)

        name = ...
            f{k};


        value = ...
            S.(name);


        path = ...
            [rootName,'.',name];


        if is_numeric_xyz_history(value)

            candidates(end+1).path = path; %#ok<AGROW>
            candidates(end).value = value;

        elseif isstruct(value) && isscalar(value)

            child = ...
                collect_numeric_xyz_fields( ...
                    value, ...
                    path);


            if ~isempty(child)

                candidates = [ ...
                    candidates, ...
                    child]; %#ok<AGROW>

            end

        end

    end

end


function x = max_finite(xIn)
%MAX_FINITE

    xIn = ...
        xIn(isfinite(xIn));


    if isempty(xIn)

        x = NaN;

    else

        x = max(xIn);

    end

end


function x = median_finite(xIn)
%MEDIAN_FINITE

    xIn = ...
        xIn(isfinite(xIn));


    if isempty(xIn)

        x = NaN;

    else

        x = median(xIn);

    end

end


function cv = coefficient_of_variation_finite(x)
%COEFFICIENT_OF_VARIATION_FINITE

    x = ...
        x(isfinite(x));


    if numel(x) < 2

        cv = NaN;
        return

    end


    mu = ...
        mean(x);


    if abs(mu) <= eps

        cv = NaN;

    else

        cv = ...
            std(x,0)/abs(mu);

    end

end
