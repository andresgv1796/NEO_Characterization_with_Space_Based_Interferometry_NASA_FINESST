% ========================================================================
% RUN_FINESST_FloquetSky_MissionAverage.m
%
% Family-wide FINESST Floquet sky-suitability study
%
%   every row of every configured JPL periodic-orbit family
%       -> complete monodromy spectrum
%       -> all nontrivial center conjugate pairs
%       -> center-branch and stability-structure analysis
%       -> structure-informed representative parent selection
%       -> one case per available center mode at every selected parent
%       -> periodic Floquet position-plane normal
%       -> one-period Floquet-normal template
%       -> repeated model-inertial mission-average accessibility
%
% NO fixed family fractions are used.
% NO planar/vertical center label is requested.
% NO competing center mode or fallback normal is substituted.
% ========================================================================

clear;
clc;
close all;

%% Plot defaults
set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
set(groot,'DefaultAxesFontSize',15);
set(groot,'DefaultTextFontSize',17);
set(groot,'DefaultLineLineWidth',1.4);
set(groot,'DefaultFigureColor','w');
set(groot,'DefaultAxesColor','w');

% ========================================================================
%% USER CONFIGURATION
% ========================================================================

systemName = 'EarthMoon';
sys = cr3bp_system_parameters(systemName);
mu = sys.mu;
studyVersion = '4.1-refined-equal-parent-visualization';

dataRoot = pwd;
familyFilePattern = '*_JPL_IC_Earth_Moon.csv';

% The family analysis is deliberately exhaustive: every JPL row is used.
familyCfg = struct();
familyCfg.analysisVersion = studyVersion;
familyCfg.stateClosureTol = 1e-8;
familyCfg.branchTrackMaxCost = 0.35;
familyCfg.branchTrackUniquenessTol = 0.02;
familyCfg.signatureExcursionMaxRows = 2;
familyCfg.topologyStateWeight = 1.0;
familyCfg.topologySpectrumWeight = 0.5;
familyCfg.topologyEdgeRatioMax = 8.0;
familyCfg.topologyNeighborCount = 12;
familyCfg.topologyMinSegmentRows = 5;
familyCfg.minParentsPerSegment = 3;
familyCfg.maxParents = 24;
familyCfg.parentCoverageTol = 0.18;
familyCfg.continuationWeight = 1.0;
familyCfg.spectrumWeight = 1.0;
familyCfg.printProgress = true;

% Center-mode enumeration. The +1 trivial pair is explicitly excluded.
floquetOpts = struct();
floquetOpts.nTrivialMultipliers = 2;
floquetOpts.trivialPairMaxDistance = 1e-3;
floquetOpts.trivialFlowAlignmentMin = 0.50;
floquetOpts.trivialMultiplierTol = 1e-5;
floquetOpts.imagTol = 1e-10;
floquetOpts.centerRadialLogTol = 5e-3;
floquetOpts.centerConjugacyTol = 5e-3;
floquetOpts.stabilityRadialLogTol = 5e-3;
floquetOpts.complexQuartetAngleTol = 5e-3;
floquetOpts.complexQuartetRadialTol = 1e-4;
floquetOpts.eigenpairResidualRelTol = 1e-10;

% Dense transport and geometry validation.
transportOpts = floquetOpts;
transportOpts.planeRankRelTol = 1e-10;
transportOpts.continuityAbsDotTol = 1e-6;
transportOpts.floquetClosureRelTol = 1e-7;
transportOpts.modeMatchMaxCost = 5e-3;
transportOpts.modeMatchUniquenessTol = 1e-6;
transportOpts.normalNormTol = 1e-10;

% One periodic template is propagated once. Mission averages are constructed
% by deterministic low-discrepancy sampling of that repeated template.
nPrimaryPeriods = 1;
nSavePerPrimaryPeriod = 500;
synodicPhase0_rad = 0.0;
missionDurationYears = 10.0;
missionAverageSamples = 6000;
convergenceDurationYears = 5.0;
convergenceAverageSamples = 3000;
earthMoonSynodicPeriodDays = 27.321661;
daysPerJulianYear = 365.25;
yearsPerTU = earthMoonSynodicPeriodDays/(2*pi*daysPerJulianYear);
missionDurationTU = missionDurationYears/yearsPerTU;
convergenceDurationTU = convergenceDurationYears/yearsPerTU;

% Sky and normal-accessibility metric.
nTargets = 10000;
alphaThresholdsDeg = [10 20 30];
alphaUsefulDeg = 30.0;

if nPrimaryPeriods ~= 1
    error('FINESST:AllModesStudy:TemplateMustBe1Tp', ...
        'The Floquet-normal template must contain exactly one primary period.');
end
if convergenceDurationYears >= missionDurationYears
    error('FINESST:AllModesStudy:InvalidConvergenceHorizon', ...
        'The convergence-check duration must be shorter than the mission duration.');
end
if ~isequal(alphaThresholdsDeg,[10 20 30])
    error('FINESST:AllModesStudy:RequiredAngularDiagnostics', ...
        'This driver''s fixed summary columns require [10 20 30] deg.');
end

% Explicit integrator; no solver fallback is used.
odeOpts = odeset('RelTol',1e-12,'AbsTol',1e-12);

% Output.
resultsDir = fullfile(pwd,'FINESST_FloquetSky_RefinedStudy_results');
saveResults = true;
saveFamilyFigures = true;
savePerModeSkyFigures = false;
saveFamilySkyFigures = true;
saveEditableFigures = true;
reuseSavedFamilyAnalysis = true;
keepFiguresOpen = true;
missionTag = local_duration_tag(missionDurationYears);
generatedFigureHandles = gobjects(0,1);

if ~isfolder(resultsDir)
    mkdir(resultsDir);
end
familyResultsDir = fullfile(resultsDir,'family_floquet_structure');
if ~isfolder(familyResultsDir)
    mkdir(familyResultsDir);
end
targetMapDir = fullfile(resultsDir,'target_accessibility_maps');
perModeMapDir = fullfile(targetMapDir,'per_mode');
familyMapDir = fullfile(targetMapDir,'per_family');
familyPathDir = fullfile(targetMapDir,'family_normal_paths');
if (savePerModeSkyFigures || saveFamilySkyFigures || ...
        saveResults) && ...
        ~isfolder(targetMapDir)
    mkdir(targetMapDir);
end
if savePerModeSkyFigures && ~isfolder(perModeMapDir)
    mkdir(perModeMapDir);
end
if (saveFamilySkyFigures || saveResults) && ~isfolder(familyMapDir)
    mkdir(familyMapDir);
end
if (saveFamilySkyFigures || saveResults) && ~isfolder(familyPathDir)
    mkdir(familyPathDir);
end

% ========================================================================
%% DISCOVER FAMILIES AND BUILD SKY
% ========================================================================

familyFiles = dir(fullfile(dataRoot,familyFilePattern));
if isempty(familyFiles)
    error('FINESST:AllModesStudy:NoFamilyFiles', ...
        'No JPL family files match "%s" under %s.',familyFilePattern,dataRoot);
end
[~,iSort] = sort(string({familyFiles.name}));
familyFiles = familyFiles(iSort);

targets = fibonacci_sphere_targets(nTargets);

fprintf('\n============================================================\n');
fprintf('FINESST ALL-CENTER-MODE FLOQUET SKY STUDY\n');
fprintf('============================================================\n');
fprintf('System                       : %s\n',systemName);
fprintf('JPL family files             : %d\n',numel(familyFiles));
fprintf('Family analysis              : every JPL row\n');
fprintf('Parent selection             : Floquet-structure coverage\n');
fprintf('Family topology              : graph recovery; source-row order ignored\n');
fprintf('Center-mode policy           : enumerate every nontrivial pair\n');
fprintf('Sky frame                    : model inertial; J2000-aligned at phase zero\n');
fprintf('Floquet template             : %d primary period\n',nPrimaryPeriods);
fprintf('Mission-average horizon      : %.3f years = %.6f TU\n', ...
    missionDurationYears,missionDurationTU);
fprintf('Convergence-check horizon    : %.3f years = %.6f TU\n', ...
    convergenceDurationYears,convergenceDurationTU);
fprintf('Mission / check samples      : %d / %d\n', ...
    missionAverageSamples,convergenceAverageSamples);
fprintf('Midpoint samples / period    : %d\n',nSavePerPrimaryPeriod);
fprintf('Fibonacci targets            : %d\n',nTargets);
fprintf('============================================================\n');

% ========================================================================
%% STAGE 1: COMPLETE FAMILY FLOQUET ANALYSIS AND PARENT SELECTION
% ========================================================================

familyAnalyses = cell(numel(familyFiles),1);

emptyCase = struct( ...
    'familyIndex',{}, ...
    'familyFile',{}, ...
    'familyName',{}, ...
    'parentRow',{}, ...
    'continuationSegmentId',{}, ...
    'segmentCoordinate',{}, ...
    'parentSelectionReason',{}, ...
    'TPrimary',{}, ...
    'xPO0',{}, ...
    'branchId',{}, ...
    'segmentBranchId',{}, ...
    'localModeIndex',{}, ...
    'modeReference',{});
cases = emptyCase;

for iFamily = 1:numel(familyFiles)
    familyPath = fullfile(familyFiles(iFamily).folder,familyFiles(iFamily).name);
    [ICfamily,Tfamily] = load_JPL_orbit_data(familyPath);
    [~,familyName] = fileparts(familyFiles(iFamily).name);

    fprintf('\n============================================================\n');
    fprintf('FAMILY %d / %d: %s\n',iFamily,numel(familyFiles),familyName);
    fprintf('============================================================\n');

    familyCfgCurrent = familyCfg;
    familyCfgCurrent.familyName = familyName;

    analysisPath = fullfile(familyResultsDir, ...
        sprintf('%s__floquet_analysis.mat',familyName));
    reusedAnalysis = false;
    if reuseSavedFamilyAnalysis && isfile(analysisPath)
        cached = load(analysisPath,'analysis');
        if isfield(cached,'analysis') && local_cache_matches( ...
                cached.analysis,ICfamily,Tfamily,floquetOpts,familyCfgCurrent)
            analysis = cached.analysis;
            reusedAnalysis = true;
            fprintf('  REUSED validated saved all-row Floquet analysis.\n');
        else
            fprintf('  Saved analysis is incompatible; recomputing family.\n');
        end
    end
    if ~reusedAnalysis
        analysis = analyze_cr3bp_family_floquet_structure( ...
            ICfamily,Tfamily,mu,odeOpts,floquetOpts,familyCfgCurrent);
    end
    familyAnalyses{iFamily} = analysis;

    fprintf('  analyzed family rows       : %d\n',height(analysis.scanTable));
    fprintf('  selected parent rows       : %d\n',numel(analysis.candidateParents));
    fprintf('  recovered continuation segments: %d\n', ...
        analysis.topology.nContinuationSegments);
    fprintf('  center modes across family : %d\n',height(analysis.modeTable));
    fprintf('  tracked center branches    : %d\n', ...
        analysis.branchDiagnostics.nBranches);
    fprintf('  trivial-pair warnings      : %d\n', ...
        nnz(~analysis.scanTable.TrivialPairResolved));

    writetable(analysis.scanTable,fullfile(familyResultsDir, ...
        sprintf('%s__family_scan.csv',familyName)));
    if ~isempty(analysis.modeTable)
        writetable(analysis.modeTable,fullfile(familyResultsDir, ...
            sprintf('%s__center_modes.csv',familyName)));
    end
    save(analysisPath,'analysis','-v7.3');

    if saveFamilyFigures
        fig = local_plot_family_floquet_structure(analysis);
        generatedFigureHandles(end+1,1) = fig;
        local_export_figure(fig,fullfile(familyResultsDir, ...
            sprintf('%s__floquet_structure.png',familyName)),300, ...
            saveEditableFigures);
        if ~keepFiguresOpen
            close(fig);
        end
    end

    for iParent = 1:numel(analysis.candidateParents)
        parent = analysis.candidateParents(iParent);
        for iMode = 1:numel(parent.centerModes)
            j = numel(cases)+1;
            cases(j).familyIndex = iFamily;
            cases(j).familyFile = familyPath;
            cases(j).familyName = familyName;
            cases(j).parentRow = parent.row;
            cases(j).continuationSegmentId = parent.continuationSegmentId;
            cases(j).segmentCoordinate = parent.segmentCoordinate;
            cases(j).parentSelectionReason = parent.selectionReason;
            cases(j).TPrimary = parent.TPrimary;
            cases(j).xPO0 = parent.state;
            cases(j).branchId = parent.centerModes(iMode).branchId;
            cases(j).segmentBranchId = ...
                parent.centerModes(iMode).segmentBranchId;
            cases(j).localModeIndex = parent.centerModes(iMode).localIndex;
            cases(j).modeReference = parent.centerModes(iMode);
        end
    end
end

nCases = numel(cases);
fprintf('\n============================================================\n');
fprintf('FAMILY ANALYSIS COMPLETE\n');
fprintf('Expanded parent-mode cases: %d\n',nCases);
fprintf('Target-map directory      : %s\n',targetMapDir);
fprintf('============================================================\n');

if nCases == 0
    error('FINESST:AllModesStudy:NoCenterModeCases', ...
        'No selected family parent possesses a nontrivial center mode.');
end

% ========================================================================
%% STAGE 2: TRANSPORT EVERY SELECTED CENTER MODE AND SCAN THE SKY
% ========================================================================

nThresholds = numel(alphaThresholdsDeg);
continuousScoreAll = nan(nTargets,nCases);
fractionWithinAll = nan(nTargets,nThresholds,nCases);
convergenceScoreAll = nan(nTargets,nCases);
convergenceF30All = nan(nTargets,nCases);
alphaMinAll = nan(nTargets,nCases);
alphaMeanAll = nan(nTargets,nCases);
alphaRMSAll = nan(nTargets,nCases);
normalPathAll = cell(nCases,1);
normalSynodicPathAll = cell(nCases,1);
missionNormalSamplesAll = cell(nCases,1);
normalTimeAll = cell(nCases,1);
chiPathAll = cell(nCases,1);
parentOrbitSynodicAll = cell(nCases,1);
bestTargetIndexAll = nan(nCases,1);
status = repmat("not-run",nCases,1);
statusMessage = strings(nCases,1);
summaryRows = cell(nCases,1);

for iCase = 1:nCases
    c = cases(iCase);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('MODE CASE %d / %d\n',iCase,nCases);
    fprintf('Family / segment / parent / branch : %s / %d / %d / %d\n', ...
        c.familyName,c.continuationSegmentId,c.parentRow,c.segmentBranchId);
    fprintf('Segment coordinate     : %.8f\n',c.segmentCoordinate);
    fprintf('TPrimary                 : %.12f TU\n',c.TPrimary);
    fprintf('Floquet template         : %.12f TU = 1 TPrimary\n',c.TPrimary);
    fprintf('Mission average          : %.6f years = %.6f TU\n', ...
        missionDurationYears,missionDurationTU);

    try
        floquet = transport_cr3bp_floquet_center_mode( ...
            c.xPO0,c.TPrimary,mu,c.modeReference,nPrimaryPeriods, ...
            nSavePerPrimaryPeriod,synodicPhase0_rad,odeOpts,transportOpts);
    catch ME
        if strcmp(ME.identifier,'FINESST:ModeTransport:PositionPlaneRankDeficient')
            status(iCase) = "position-plane-inadmissible";
            statusMessage(iCase) = string(ME.message);
            fprintf('POSITION-PLANE-INADMISSIBLE: %s\n',ME.message);
            continue
        elseif startsWith(string(ME.identifier),"FINESST:ModeTransport:")
            status(iCase) = "floquet-transport-inadmissible";
            statusMessage(iCase) = string(ME.message);
            fprintf('FLOQUET-TRANSPORT-INADMISSIBLE: %s\n',ME.message);
            continue
        end
        rethrow(ME)
    end

    mission = sample_repeated_floquet_normal_mission( ...
        floquet.nSynodicOnePeriod,c.TPrimary,missionDurationTU, ...
        missionAverageSamples,synodicPhase0_rad);
    convergenceMission = sample_repeated_floquet_normal_mission( ...
        floquet.nSynodicOnePeriod,c.TPrimary,convergenceDurationTU, ...
        convergenceAverageSamples,synodicPhase0_rad);

    scan = evaluate_floquet_normal_sky_suitability( ...
        mission.nInertial,targets.sInertial, ...
        alphaThresholdsDeg,alphaUsefulDeg);
    convergenceScan = evaluate_floquet_normal_sky_suitability( ...
        convergenceMission.nInertial,targets.sInertial, ...
        alphaThresholdsDeg,alphaUsefulDeg);

    continuousScoreAll(:,iCase) = scan.continuousScore;
    fractionWithinAll(:,:,iCase) = scan.fractionWithin;
    convergenceScoreAll(:,iCase) = convergenceScan.continuousScore;
    convergenceF30All(:,iCase) = convergenceScan.fractionWithin(:,3);
    alphaMinAll(:,iCase) = scan.alphaMinDeg;
    alphaMeanAll(:,iCase) = scan.alphaMeanDeg;
    alphaRMSAll(:,iCase) = scan.alphaRMSDeg;
    normalPathAll{iCase} = floquet.nInertial;
    normalSynodicPathAll{iCase} = floquet.nSynodicOnePeriod;
    missionNormalSamplesAll{iCase} = mission.nInertial;
    normalTimeAll{iCase} = floquet.t;
    chiPathAll{iCase} = floquet.chi;
    parentOrbitSynodicAll{iCase} = floquet.parentOrbitSynodic;
    status(iCase) = "valid";

    [bestScore,bestTargetIdx] = max(scan.continuousScore);
    bestTargetIndexAll(iCase) = bestTargetIdx;

    fprintf('  lambda                   = %.12g %+.12gi\n', ...
        real(floquet.lambda),imag(floquet.lambda));
    fprintf('  rho/(2pi)                = %.10f cycles\n',floquet.rho/(2*pi));
    fprintf('  planar / vertical proj.  = %.8f / %.8f\n', ...
        floquet.selectedMode.planarProjection, ...
        floquet.selectedMode.verticalProjection);
    fprintf('  unit-circle defect       = %.6e\n',floquet.unitCircleDefect);
    fprintf('  eigenpair residual       = %.6e\n',floquet.eigenpairResidualRel);
    fprintf('  Floquet closure          = %.6e\n',floquet.floquetClosureRel);
    fprintf('  chi min / median         = %.6e / %.6e\n', ...
        min(floquet.chiOnePeriod),median(floquet.chiOnePeriod));
    fprintf('  mission-best lon / lat   = %.6f / %.6f deg\n', ...
        targets.longitudeDeg(bestTargetIdx), ...
        targets.latitudeDeg(bestTargetIdx));
    fprintf('  best normal score        = %.8f\n',bestScore);

    summaryRows{iCase} = table( ...
        string(c.familyName),c.continuationSegmentId,c.parentRow, ...
        c.segmentCoordinate,string(c.parentSelectionReason), ...
        c.segmentBranchId,c.branchId,c.localModeIndex, ...
        c.TPrimary,missionDurationYears,real(floquet.lambda),imag(floquet.lambda), ...
        floquet.rho/(2*pi),floquet.selectedMode.planarProjection, ...
        floquet.selectedMode.verticalProjection,floquet.unitCircleDefect, ...
        floquet.eigenpairResidualRel,floquet.floquetClosureRel, ...
        floquet.stateClosureNorm,floquet.modeMatchCost, ...
        floquet.periodicSeamStepDeg,min(floquet.chiOnePeriod), ...
        median(floquet.chiOnePeriod),targets.longitudeDeg(bestTargetIdx), ...
        targets.latitudeDeg(bestTargetIdx),bestScore, ...
        scan.fractionWithin(bestTargetIdx,1), ...
        scan.fractionWithin(bestTargetIdx,2), ...
        scan.fractionWithin(bestTargetIdx,3), ...
        scan.alphaMinDeg(bestTargetIdx), ...
        scan.alphaMeanDeg(bestTargetIdx), ...
        scan.alphaRMSDeg(bestTargetIdx), ...
        'VariableNames',{ ...
            'Family','ContinuationSegmentId','ParentRow','SegmentCoordinate', ...
            'ParentSelectionReason','SegmentBranchId','CenterBranchId', ...
            'LocalModeIndex','TPrimary_TU','MissionDuration_years', ...
            'LambdaReal','LambdaImag','RhoOver2Pi', ...
            'PlanarProjection','VerticalProjection','UnitCircleDefect', ...
            'EigenpairResidualRel','FloquetClosureRel','POStateClosureNorm', ...
            'ModeMatchCost','PeriodicSeamStep_deg','ChiMin','ChiMedian', ...
            'MissionBestModelLongitude_deg','MissionBestModelLatitude_deg', ...
            'BestContinuousScore', ...
            'BestF10','BestF20','BestF30', ...
            'BestAlphaMin_deg','BestAlphaMean_deg','BestAlphaRMS_deg'});

    if savePerModeSkyFigures
        i30 = find(abs(alphaThresholdsDeg-30) <= 10*eps(30),1);
        fig = figure('Color','w','Position',[100 60 1250 950]);
        generatedFigureHandles(end+1,1) = fig;
        tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

        nexttile;
        scatter(targets.longitudeDeg,targets.latitudeDeg,18, ...
            scan.continuousScore,'filled');
        hold on;
        local_plot_axis_marker(gca,targets.sInertial(bestTargetIdx,:), ...
            [0.85 0.05 0.05],10,'p');
        local_format_sky_axes();
        title(sprintf( ...
            ['%s, segment %d, row %d, branch %d: %.3g-year mission average ' ...
             '(red: mission optimum)'], ...
            strrep(c.familyName,'_','\_'),c.continuationSegmentId,c.parentRow, ...
            c.segmentBranchId,missionDurationYears));
        cb = colorbar;
        ylabel(cb,'Continuous normal suitability');
        clim([0 1]);

        nexttile;
        scatter(targets.longitudeDeg,targets.latitudeDeg,18, ...
            scan.fractionWithin(:,i30),'filled');
        hold on;
        local_plot_axis_marker(gca,targets.sInertial(bestTargetIdx,:), ...
            [0.85 0.05 0.05],10,'p');
        local_format_sky_axes();
        title(sprintf('$f_{30}$ over %.3g years',missionDurationYears));
        cb = colorbar;
        ylabel(cb,'$f_{30}$');
        clim([0 1]);

        local_export_figure(fig,fullfile(perModeMapDir,sprintf( ...
            '%s__segment_%02d__row_%04d__branch_%02d__target_accessibility_%s.png', ...
            c.familyName,c.continuationSegmentId,c.parentRow, ...
            c.segmentBranchId,missionTag)),250,saveEditableFigures);
        if ~keepFiguresOpen
            close(fig);
        end
    end
end

% ========================================================================
%% STATUS, FAMILY SUMMARIES, AND SAVE
% ========================================================================

validSummary = ~cellfun(@isempty,summaryRows);
if any(validSummary)
    summaryTable = vertcat(summaryRows{validSummary});
else
    summaryTable = table();
end

caseStatusTable = table( ...
    string({cases.familyName}).',[cases.parentRow].', ...
    [cases.continuationSegmentId].',[cases.segmentCoordinate].', ...
    [cases.segmentBranchId].',[cases.branchId].', ...
    [cases.localModeIndex].', ...
    status,statusMessage, ...
    'VariableNames',{ ...
        'Family','ParentRow','ContinuationSegmentId','SegmentCoordinate', ...
        'SegmentBranchId','CenterBranchId','LocalModeIndex', ...
        'Status','StatusMessage'});

validCaseMask = status == "valid";

% Make family inclusion and exclusion explicit. In particular, a family with
% center modes but no admissible fixed-parent-phase position plane remains in
% the study-status table rather than disappearing from the figures silently.
nFamilies = numel(familyFiles);
familyStudyStatusRows = cell(nFamilies,1);
for iFamily = 1:nFamilies
    familyAllMask = [cases.familyIndex].' == iFamily;
    familyRejectedStatus = status(familyAllMask & ~validCaseMask);
    familyStudyStatusRows{iFamily} = table( ...
        string(familyAnalyses{iFamily}.familyName), ...
        height(familyAnalyses{iFamily}.scanTable), ...
        numel(familyAnalyses{iFamily}.candidateParents), ...
        nnz(familyAllMask),nnz(familyAllMask & validCaseMask), ...
        nnz(familyAllMask & ~validCaseMask), ...
        local_main_rejection_reason(familyRejectedStatus), ...
        nnz(familyAllMask & validCaseMask) > 0, ...
        'VariableNames',{ ...
            'Family','ScannedRows','SelectedParents','CenterModeCases', ...
            'ValidPlaneModels','RejectedModels','MainRejectionReason', ...
            'SkyMapsAvailable'});
end
familyStudyStatusTable = vertcat(familyStudyStatusRows{:});
fprintf('\nFAMILY STUDY STATUS (families without maps remain explicit):\n');
disp(familyStudyStatusTable);

idx30 = find(abs(alphaThresholdsDeg-30) <= 10*eps(30),1);
if isempty(idx30)
    error('FINESST:AllModesStudy:Missing30DegThreshold', ...
        'alphaThresholdsDeg must explicitly contain 30 deg.');
end

% Family-level accessibility products. Modes are first maximized within each
% parent. The opportunity envelope then maximizes across parents, whereas the
% equal-parent maps average the parent maps so every representative parent has
% one vote. The two summaries answer deliberately different questions.
familyScore = nan(nTargets,nFamilies);
familyWinnerCase = nan(nTargets,nFamilies);
familyF30 = nan(nTargets,nFamilies);
familyWinnerF30Case = nan(nTargets,nFamilies);
familyMeanScore = nan(nTargets,nFamilies);
familyScoreStd = nan(nTargets,nFamilies);
familyMeanF30 = nan(nTargets,nFamilies);
familyPathSummaryRows = cell(nFamilies,1);
familyParentSummaryRows = cell(nFamilies,1);
familyPreferenceRows = cell(nFamilies,1);
familyDiagnosticsRows = cell(nFamilies,1);
familySegmentCoverageRows = cell(nFamilies,1);

for iFamily = 1:nFamilies
    familyCaseMask = validCaseMask & ([cases.familyIndex].' == iFamily);
    if ~any(familyCaseMask)
        fprintf('FAMILY MAP UNAVAILABLE: %s has no valid transported modes.\n', ...
            familyAnalyses{iFamily}.familyName);
        continue
    end

    familyCaseIndex = find(familyCaseMask);
    familyName = familyAnalyses{iFamily}.familyName;
    parentRows = unique([cases(familyCaseIndex).parentRow].','sorted');
    nParents = numel(parentRows);
    parentScore = nan(nTargets,nParents);
    parentF30 = nan(nTargets,nParents);
    parentConvergenceScore = nan(nTargets,nParents);
    parentConvergenceF30 = nan(nTargets,nParents);
    parentWinnerCase = nan(nTargets,nParents);
    parentWinnerF30Case = nan(nTargets,nParents);
    parentCaseIndex = cell(nParents,1);
    parentBestTargetIndex = nan(nParents,1);
    parentBestScore = nan(nParents,1);
    parentF30AtContinuousBest = nan(nParents,1);
    parentBestF30TargetIndex = nan(nParents,1);
    parentBestF30Score = nan(nParents,1);
    parentBestCase = nan(nParents,1);

    for jParent = 1:nParents
        inParent = [cases(familyCaseIndex).parentRow].' == parentRows(jParent);
        parentCaseIndex{jParent} = familyCaseIndex(inParent);
        [parentScore(:,jParent),winnerMode] = max( ...
            continuousScoreAll(:,parentCaseIndex{jParent}),[],2);
        parentWinnerCase(:,jParent) = parentCaseIndex{jParent}(winnerMode);
        f30Parent = reshape( ...
            fractionWithinAll(:,idx30,parentCaseIndex{jParent}),nTargets,[]);
        [parentF30(:,jParent),winnerF30Mode] = max(f30Parent,[],2);
        parentWinnerF30Case(:,jParent) = ...
            parentCaseIndex{jParent}(winnerF30Mode);
        parentConvergenceScore(:,jParent) = max( ...
            convergenceScoreAll(:,parentCaseIndex{jParent}),[],2);
        parentConvergenceF30(:,jParent) = max( ...
            convergenceF30All(:,parentCaseIndex{jParent}),[],2);
        [parentBestScore(jParent),parentBestTargetIndex(jParent)] = ...
            max(parentScore(:,jParent));
        parentBestCase(jParent) = parentWinnerCase( ...
            parentBestTargetIndex(jParent),jParent);
        parentF30AtContinuousBest(jParent) = parentF30( ...
            parentBestTargetIndex(jParent),jParent);
        [parentBestF30Score(jParent),parentBestF30TargetIndex(jParent)] = ...
            max(parentF30(:,jParent));
    end

    [familyScore(:,iFamily),winnerParent] = max(parentScore,[],2);
    linearWinner = sub2ind([nTargets nParents],(1:nTargets).',winnerParent);
    familyWinnerCase(:,iFamily) = parentWinnerCase(linearWinner);
    familyTieCount = sum(abs(parentScore-familyScore(:,iFamily)) <= 1e-12,2);
    [familyF30(:,iFamily),winnerF30Parent] = max(parentF30,[],2);
    linearWinnerF30 = sub2ind( ...
        [nTargets nParents],(1:nTargets).',winnerF30Parent);
    familyWinnerF30Case(:,iFamily) = parentWinnerF30Case(linearWinnerF30);
    familyF30TieCount = sum(abs(parentF30-familyF30(:,iFamily)) <= 1e-12,2);
    familyMeanScore(:,iFamily) = mean(parentScore,2);
    familyScoreStd(:,iFamily) = std(parentScore,0,2);
    familyMeanF30(:,iFamily) = mean(parentF30,2);

    convergenceFamilyScore = max(parentConvergenceScore,[],2);
    convergenceFamilyF30 = max(parentConvergenceF30,[],2);
    convergenceFamilyMeanScore = mean(parentConvergenceScore,2);
    convergenceFamilyMeanF30 = mean(parentConvergenceF30,2);
    scoreConvergenceRMS = sqrt(mean( ...
        (familyScore(:,iFamily)-convergenceFamilyScore).^2));
    scoreConvergenceMax = max(abs( ...
        familyScore(:,iFamily)-convergenceFamilyScore));
    f30ConvergenceRMS = sqrt(mean( ...
        (familyF30(:,iFamily)-convergenceFamilyF30).^2));
    f30ConvergenceMax = max(abs( ...
        familyF30(:,iFamily)-convergenceFamilyF30));
    continuousVsF30RMS = sqrt(mean( ...
        (familyScore(:,iFamily)-familyF30(:,iFamily)).^2));
    sameWinningParentFraction = mean(winnerParent == winnerF30Parent);
    meanScoreConvergenceRMS = sqrt(mean( ...
        (familyMeanScore(:,iFamily)-convergenceFamilyMeanScore).^2));
    meanF30ConvergenceRMS = sqrt(mean( ...
        (familyMeanF30(:,iFamily)-convergenceFamilyMeanF30).^2));
    familyDiagnosticsRows{iFamily} = table( ...
        string(familyName),convergenceDurationYears,missionDurationYears, ...
        scoreConvergenceRMS,scoreConvergenceMax, ...
        f30ConvergenceRMS,f30ConvergenceMax,continuousVsF30RMS, ...
        sameWinningParentFraction,meanScoreConvergenceRMS, ...
        meanF30ConvergenceRMS, ...
        'VariableNames',{ ...
            'Family','ShortHorizon_years','LongHorizon_years', ...
            'ContinuousEnvelopeRMSDifference', ...
            'ContinuousEnvelopeMaxDifference', ...
            'F30EnvelopeRMSDifference','F30EnvelopeMaxDifference', ...
            'ContinuousVsF30EnvelopeRMSDifference', ...
            'SameWinningParentFraction', ...
            'EqualParentMeanContinuousRMSDifference', ...
            'EqualParentMeanF30RMSDifference'});

    winnerCase = familyWinnerCase(:,iFamily);
    winnerF30Case = familyWinnerF30Case(:,iFamily);
    familySkyTable = table( ...
        (1:nTargets).',targets.longitudeDeg,targets.latitudeDeg, ...
        familyMeanScore(:,iFamily),familyScoreStd(:,iFamily), ...
        familyMeanF30(:,iFamily), ...
        familyScore(:,iFamily),familyTieCount,winnerCase, ...
        [cases(winnerCase).continuationSegmentId].', ...
        [cases(winnerCase).parentRow].', ...
        [cases(winnerCase).segmentBranchId].', ...
        familyF30(:,iFamily),familyF30TieCount,winnerF30Case, ...
        [cases(winnerF30Case).continuationSegmentId].', ...
        [cases(winnerF30Case).parentRow].', ...
        [cases(winnerF30Case).segmentBranchId].', ...
        'VariableNames',{ ...
            'TargetIndex','ModelLongitude_deg','ModelLatitude_deg', ...
            'EqualParentMeanContinuousScore','ParentContinuousStd', ...
            'EqualParentMeanF30', ...
            'OpportunityEnvelopeContinuousScore','ContinuousTieCount', ...
            'ExampleContinuousWinnerCase','ExampleContinuousWinnerSegment', ...
            'ExampleContinuousWinnerParentRow', ...
            'ExampleContinuousWinnerSegmentBranch', ...
            'OpportunityEnvelopeF30','F30TieCount','ExampleF30WinnerCase', ...
            'ExampleF30WinnerSegment','ExampleF30WinnerParentRow', ...
            'ExampleF30WinnerSegmentBranch'});

    modeBestIdx = bestTargetIndexAll(familyCaseIndex);
    parentBestDirections = targets.sInertial(parentBestTargetIndex,:);
    [preferredAxis,axisConcentration,axisMedianDeg,axisRMSDeg] = ...
        local_axial_preferred_direction(parentBestDirections);
    preferredLongitude = mod(atan2d(preferredAxis(2),preferredAxis(1)),360);
    preferredLatitude = asind(preferredAxis(3));

    chiMinimum = cellfun(@min,chiPathAll(familyCaseIndex));
    chiMedian = cellfun(@median,chiPathAll(familyCaseIndex));
    familyPathSummary = table( ...
        familyCaseIndex,[cases(familyCaseIndex).continuationSegmentId].', ...
        [cases(familyCaseIndex).parentRow].', ...
        [cases(familyCaseIndex).segmentBranchId].', ...
        [cases(familyCaseIndex).branchId].', ...
        [cases(familyCaseIndex).localModeIndex].', ...
        [cases(familyCaseIndex).segmentCoordinate].', ...
        targets.longitudeDeg(modeBestIdx),targets.latitudeDeg(modeBestIdx), ...
        arrayfun(@(i) continuousScoreAll(modeBestIdx(i),familyCaseIndex(i)), ...
            (1:numel(familyCaseIndex)).'), ...
        chiMinimum(:),chiMedian(:), ...
        'VariableNames',{ ...
            'CaseIndex','ContinuationSegmentId','ParentRow','SegmentBranchId', ...
            'CenterBranchId','LocalModeIndex','SegmentCoordinate', ...
            'MissionBestModelLongitude_deg', ...
            'MissionBestModelLatitude_deg','MissionBestContinuousScore', ...
            'ChiMin','ChiMedian'});
    familyPathSummaryRows{iFamily} = familyPathSummary;

    parentFirstCase = cellfun(@(v) v(1),parentCaseIndex);
    parentPeriods = [cases(parentFirstCase).TPrimary].';
    parentSummary = table( ...
        repmat(string(familyName),nParents,1),parentRows, ...
        [cases(parentFirstCase).continuationSegmentId].', ...
        cellfun(@numel,parentCaseIndex), ...
        parentPeriods, ...
        targets.longitudeDeg(parentBestTargetIndex), ...
        targets.latitudeDeg(parentBestTargetIndex), ...
        parentBestScore,parentF30AtContinuousBest, ...
        targets.longitudeDeg(parentBestF30TargetIndex), ...
        targets.latitudeDeg(parentBestF30TargetIndex),parentBestF30Score, ...
        parentBestCase, ...
        [cases(parentBestCase).segmentBranchId].', ...
        [cases(parentBestCase).branchId].', ...
        string({cases(parentFirstCase).parentSelectionReason}).', ...
        'VariableNames',{ ...
            'Family','ParentRow','ContinuationSegmentId','ValidModes', ...
            'TPrimary_TU','MissionBestModelLongitude_deg', ...
            'MissionBestModelLatitude_deg','MissionBestContinuousScore', ...
            'MissionF30AtContinuousBest','F30BestModelLongitude_deg', ...
            'F30BestModelLatitude_deg','MissionBestF30','BestCaseIndex', ...
            'BestSegmentBranchId','BestCenterBranchId', ...
            'ParentSelectionReason'});
    familyParentSummaryRows{iFamily} = parentSummary;
    familyPreferenceRows{iFamily} = table( ...
        string(familyName),nParents,numel(familyCaseIndex),preferredLongitude, ...
        preferredLatitude,axisConcentration,axisMedianDeg,axisRMSDeg, ...
        'VariableNames',{ ...
            'Family','ValidRepresentativeParents','ValidRepresentativeModes', ...
            'PreferredAxisLongitude_deg', ...
            'PreferredAxisLatitude_deg','AxialConcentration', ...
            'MedianBestDirectionOffset_deg','RMSBestDirectionOffset_deg'});

    scanTable = familyAnalyses{iFamily}.scanTable;
    continuationSegmentIds = (1:familyAnalyses{iFamily}.topology.nContinuationSegments).';
    segmentCoverageRows = cell(numel(continuationSegmentIds),1);
    for jSegment = 1:numel(continuationSegmentIds)
        continuationSegmentId = continuationSegmentIds(jSegment);
        sourceRows = scanTable.ContinuationSegmentId == continuationSegmentId;
        selectedRows = sourceRows & scanTable.SelectedParent;
        validParentMask = [cases(familyCaseIndex).continuationSegmentId].' == ...
            continuationSegmentId;
        Xsegment = familyAnalyses{iFamily}.ICfamily(sourceRows,:);
        Tsegment = familyAnalyses{iFamily}.Tfamily(sourceRows);
        segmentCoverageRows{jSegment} = table( ...
            string(familyName),continuationSegmentId,nnz(sourceRows), ...
            nnz(scanTable.CenterModeCount(sourceRows) > 0), ...
            nnz(selectedRows),numel(unique( ...
                [cases(familyCaseIndex(validParentMask)).parentRow])), ...
            min(Tsegment),max(Tsegment),min(Xsegment(:,3)), ...
            median(Xsegment(:,3)),max(Xsegment(:,3)), ...
            'VariableNames',{ ...
                'Family','ContinuationSegmentId','FamilyRows', ...
                'CenterBearingRows','SelectedParentRows', ...
                'ValidRepresentativeParents','TPrimaryMin_TU', ...
                'TPrimaryMax_TU','InitialZMin_ND','InitialZMedian_ND', ...
                'InitialZMax_ND'});
    end
    familySegmentCoverage = vertcat(segmentCoverageRows{:});
    familySegmentCoverageRows{iFamily} = familySegmentCoverage;

    if saveResults
        writetable(familySkyTable,fullfile(familyMapDir,sprintf( ...
            '%s__family_target_accessibility_%s.csv',familyName,missionTag)));
        writetable(familyPathSummary,fullfile(familyPathDir,sprintf( ...
            '%s__mode_best_directions_%s.csv',familyName,missionTag)));
        writetable(parentSummary,fullfile(familyPathDir,sprintf( ...
            '%s__parent_best_directions_%s.csv',familyName,missionTag)));
        writetable(familySegmentCoverage,fullfile(familyPathDir,sprintf( ...
            '%s__continuation_segment_coverage.csv',familyName)));
    end

    if saveFamilySkyFigures
        fig = local_plot_family_summary( ...
            parentRows,parentCaseIndex,parentOrbitSynodicAll, ...
            parentPeriods,targets,familyMeanScore(:,iFamily), ...
            familyScore(:,iFamily),familyMeanF30(:,iFamily), ...
            parentBestDirections,preferredAxis,familyName, ...
            missionDurationYears);
        generatedFigureHandles(end+1,1) = fig;
        local_export_figure(fig,fullfile(familyMapDir,sprintf( ...
            '%s__family_target_accessibility_%s.png',familyName,missionTag)), ...
            300,saveEditableFigures);
        if ~keepFiguresOpen
            close(fig);
        end

        fig = local_plot_family_geometry( ...
            parentCaseIndex,normalSynodicPathAll, ...
            missionNormalSamplesAll,parentOrbitSynodicAll,parentPeriods, ...
            parentBestDirections,preferredAxis, ...
            axisConcentration,axisMedianDeg,familyName,missionDurationYears);
        generatedFigureHandles(end+1,1) = fig;
        local_export_figure(fig,fullfile(familyPathDir,sprintf( ...
            '%s__parent_orbits_normal_paths_and_best_directions_%s.png', ...
            familyName,missionTag)),300,saveEditableFigures);
        if ~keepFiguresOpen
            close(fig);
        end
    end
end

fprintf('\nNo cross-family maximum envelope is generated.\n');
fprintf(['Families use a common %.3g-year physical horizon and remain separate ' ...
    'until a Fourier-plane metric is introduced.\n'],missionDurationYears);

validPreferenceRows = ~cellfun(@isempty,familyPreferenceRows);
if any(validPreferenceRows)
    familyPreferenceTable = vertcat(familyPreferenceRows{validPreferenceRows});
else
    familyPreferenceTable = table();
end

validParentSummaryRows = ~cellfun(@isempty,familyParentSummaryRows);
if any(validParentSummaryRows)
    parentSummaryTable = vertcat(familyParentSummaryRows{validParentSummaryRows});
else
    parentSummaryTable = table();
end
validDiagnosticsRows = ~cellfun(@isempty,familyDiagnosticsRows);
if any(validDiagnosticsRows)
    familyDiagnosticsTable = vertcat(familyDiagnosticsRows{validDiagnosticsRows});
else
    familyDiagnosticsTable = table();
end
validSegmentCoverageRows = ~cellfun(@isempty,familySegmentCoverageRows);
if any(validSegmentCoverageRows)
    segmentCoverageTable = vertcat( ...
        familySegmentCoverageRows{validSegmentCoverageRows});
else
    segmentCoverageTable = table();
end

if saveResults
    writetable(caseStatusTable,fullfile(resultsDir, ...
        'FINESST_MissionAverage_case_status.csv'));
    writetable(familyStudyStatusTable,fullfile(resultsDir, ...
        'FINESST_family_study_status.csv'));
    if ~isempty(summaryTable)
        writetable(summaryTable,fullfile(resultsDir, ...
            'FINESST_MissionAverage_mode_summary.csv'));
    end
    if ~isempty(familyPreferenceTable)
        writetable(familyPreferenceTable,fullfile(familyPathDir, ...
            sprintf('FAMILY_equal_parent_preferred_axes_%s.csv',missionTag)));
    end
    if ~isempty(parentSummaryTable)
        writetable(parentSummaryTable,fullfile(familyPathDir, ...
            sprintf('ALL_parent_best_directions_%s.csv',missionTag)));
    end
    if ~isempty(familyDiagnosticsTable)
        writetable(familyDiagnosticsTable,fullfile(familyPathDir, ...
            sprintf('FAMILY_%gyr_to_%gyr_and_metric_diagnostics.csv', ...
            convergenceDurationYears,missionDurationYears)));
    end
    if ~isempty(segmentCoverageTable)
        writetable(segmentCoverageTable,fullfile(familyPathDir, ...
            'CONTINUATION_SEGMENT_coverage.csv'));
    end

    save(fullfile(resultsDir,'FINESST_MissionAverage.mat'), ...
        'sys','familyFiles','familyAnalyses','cases','caseStatusTable', ...
        'familyStudyStatusTable', ...
        'summaryTable','familyPreferenceTable','parentSummaryTable', ...
        'familyDiagnosticsTable','segmentCoverageTable', ...
        'familyPathSummaryRows','familyParentSummaryRows', ...
        'familySegmentCoverageRows', ...
        'targets','familyCfg','floquetOpts','transportOpts', ...
        'nPrimaryPeriods','nSavePerPrimaryPeriod','synodicPhase0_rad', ...
        'missionDurationYears','missionDurationTU','missionAverageSamples', ...
        'earthMoonSynodicPeriodDays','daysPerJulianYear','yearsPerTU', ...
        'convergenceDurationYears','convergenceDurationTU', ...
        'convergenceAverageSamples', ...
        'alphaThresholdsDeg','alphaUsefulDeg','continuousScoreAll', ...
        'fractionWithinAll','convergenceScoreAll','convergenceF30All', ...
        'alphaMinAll','alphaMeanAll','alphaRMSAll', ...
        'normalPathAll','normalSynodicPathAll','missionNormalSamplesAll', ...
        'normalTimeAll','chiPathAll','parentOrbitSynodicAll', ...
        'bestTargetIndexAll', ...
        'familyScore','familyWinnerCase','familyF30','familyWinnerF30Case', ...
        'familyMeanScore','familyScoreStd','familyMeanF30', ...
        '-v7.3');
end

if ~any(validCaseMask)
    error('FINESST:AllModesStudy:NoValidPositionPlanes', ...
        ['Center modes were found, but none defined an admissible transported ' ...
         '2-D position plane. Inspect the saved case-status table.']);
end

fprintf('\nStudy complete.\n');
fprintf('Results directory: %s\n',resultsDir);
if keepFiguresOpen
    fprintf(['Live figures retained: %d. Handles are available in ' ...
        'generatedFigureHandles.\n'],numel(generatedFigureHandles));
end


function fig = local_plot_family_floquet_structure(analysis)
    scan = analysis.scanTable;
    selected = scan.SelectedParent;
    nContinuationSegments = analysis.topology.nContinuationSegments;
    colors = lines(nContinuationSegments);

    fig = figure('Color','w','Position',[80 80 1250 900]);
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    ax1 = nexttile;
    hold on;
    hSegment = gobjects(nContinuationSegments,1);
    for iSegment = 1:nContinuationSegments
        rows = find(scan.ContinuationSegmentId == iSegment);
        [~,order] = sort(scan.SegmentCoordinate(rows));
        rows = rows(order);
        hSegment(iSegment) = semilogy( ...
            scan.SegmentCoordinate(rows),scan.SpectralRadius(rows),'-', ...
            'Color',colors(iSegment,:), ...
            'DisplayName',sprintf('segment %d',iSegment));
    end
    semilogy(scan.SegmentCoordinate(selected),scan.SpectralRadius(selected), ...
        'ko','MarkerFaceColor','y','HandleVisibility','off');
    grid on; box on;
    xlim([0 1]);
    ylabel('Spectral radius');
    title(sprintf([ ...
        '%s: %d numerical continuation segments (not physical subfamilies), ' ...
        '%d center branches'], ...
        strrep(analysis.familyName,'_','\_'),nContinuationSegments, ...
        analysis.branchDiagnostics.nBranches));
    legend(hSegment,'Location','bestoutside');

    ax2 = nexttile;
    hold on;
    if ~isempty(analysis.modeTable)
        branches = unique(analysis.modeTable.BranchId);
        for k = 1:numel(branches)
            mask = analysis.modeTable.BranchId == branches(k);
            iSegment = analysis.modeTable.ContinuationSegmentId(find(mask,1));
            [~,order] = sort(analysis.modeTable.SegmentCoordinate(mask));
            x = analysis.modeTable.SegmentCoordinate(mask);
            y = analysis.modeTable.Rho_rad(mask)/pi;
            plot(x(order),y(order),'.-', ...
                'Color',colors(iSegment,:),'HandleVisibility','off');
        end
        selectedMode = analysis.modeTable.SelectedParent;
        hSelected = plot(analysis.modeTable.SegmentCoordinate(selectedMode), ...
            analysis.modeTable.Rho_rad(selectedMode)/pi, ...
            'ko','MarkerFaceColor','y','DisplayName','selected parents');
        legend(hSelected,'Location','bestoutside');
    end
    grid on; box on;
    xlim([0 1]);
    ylabel('$\rho/\pi$');

    ax3 = nexttile;
    hold on;
    for iSegment = 1:nContinuationSegments
        rows = find(scan.ContinuationSegmentId == iSegment);
        [~,order] = sort(scan.SegmentCoordinate(rows));
        rows = rows(order);
        stairs(scan.SegmentCoordinate(rows),scan.CenterModeCount(rows), ...
            '-','Color',colors(iSegment,:),'HandleVisibility','off');
    end
    plot(scan.SegmentCoordinate(selected),scan.CenterModeCount(selected), ...
        'ro','MarkerFaceColor','r');
    grid on; box on;
    xlim([0 1]);
    xlabel('Numerical continuation-segment graph coordinate');
    ylabel('Center-mode count');
    linkaxes([ax1 ax2 ax3],'x');
end


function local_format_sky_axes()
    grid on; box on;
    xlim([0 360]); ylim([-90 90]);
    xticks(0:60:360); yticks(-90:30:90);
    xlabel('Model-inertial longitude [deg]');
    ylabel('Model-inertial latitude [deg]');
end


function fig = local_plot_family_summary( ...
        parentRows,parentCaseIndex,parentOrbitSynodicAll,parentPeriods, ...
        targets,meanScore,opportunityScore,meanF30,parentBestDirections, ...
        preferredAxis,familyName,missionDurationYears)
    nParents = numel(parentRows);
    [parentColors,parentColorMap,periodLimits] = ...
        local_parent_period_colors(parentPeriods);
    familyLabel = strrep(familyName,'_','\_');

    fig = figure('Color','w','Position',[40 40 1750 1050]);
    layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    axOrbit = nexttile(layout);
    hold(axOrbit,'on');
    for jParent = 1:nParents
        iCase = parentCaseIndex{jParent}(1);
        r = parentOrbitSynodicAll{iCase};
        plot3(axOrbit,r(1,:),r(2,:),r(3,:),'-', ...
            'Color',parentColors(jParent,:),'LineWidth',1.5, ...
            'HandleVisibility','off');
    end
    axis(axOrbit,'equal'); grid(axOrbit,'on'); box(axOrbit,'on');
    xlabel(axOrbit,'synodic $x$ [ND]');
    ylabel(axOrbit,'synodic $y$ [ND]');
    zlabel(axOrbit,'synodic $z$ [ND]');
    title(axOrbit,sprintf('Representative parents ($N_p=%d$)',nParents));
    view(axOrbit,38,24);
    colormap(axOrbit,parentColorMap);
    clim(axOrbit,periodLimits);
    cbPeriod = colorbar(axOrbit);
    ylabel(cbPeriod,'Parent period $T_p$ [TU]');

    axMean = nexttile(layout);
    local_scatter_hammer_map(axMean,targets,meanScore, ...
        'Equal-parent mean continuous suitability', ...
        'Mean continuous suitability');
    hold(axMean,'on');
    for jParent = 1:nParents
        s = local_orient_direction(parentBestDirections(jParent,:),preferredAxis);
        local_plot_hammer_marker(axMean,s,parentColors(jParent,:),7,'o');
    end
    local_plot_hammer_marker(axMean,preferredAxis,[0 0 0],12,'p');
    text(axMean,0,-1.64, ...
        'colored circles: parent optima; black star: equal-parent preferred axis', ...
        'HorizontalAlignment','center','FontSize',9,'Interpreter','none');

    axOpportunity = nexttile(layout);
    local_scatter_hammer_map(axOpportunity,targets,opportunityScore, ...
        'Best-parent opportunity envelope', ...
        'Maximum continuous suitability');

    axMeanF30 = nexttile(layout);
    local_scatter_hammer_map(axMeanF30,targets,meanF30, ...
        'Equal-parent mean $f_{30}$', ...
        'Mean $f_{30}$');

    title(layout,sprintf([ ...
        '%s: %.3g-year model-inertial observability; ' ...
        'axes J2000-aligned only at phase zero'], ...
        familyLabel,missionDurationYears));
end


function local_scatter_hammer_map(ax,targets,values,titleText,colorbarText)
    [x,y] = local_hammer_project( ...
        targets.longitudeDeg,targets.latitudeDeg);
    scatter(ax,x,y,16,values,'filled');
    local_format_hammer_axes(ax);
    colormap(ax,parula(256));
    clim(ax,[0 1]);
    cb = colorbar(ax);
    ylabel(cb,colorbarText);
    title(ax,titleText);
end


function local_format_hammer_axes(ax)
    holdState = ishold(ax);
    hold(ax,'on');
    gridColor = [0.88 0.88 0.88];
    latitudeGrid = linspace(-89.9,89.9,241);
    for longitude = -120:60:120
        [x,y] = local_hammer_project( ...
            longitude*ones(size(latitudeGrid)),latitudeGrid);
        plot(ax,x,y,'-','Color',gridColor, ...
            'LineWidth',0.55,'HandleVisibility','off');
    end
    longitudeGrid = linspace(-180,180,481);
    for latitude = -60:30:60
        [x,y] = local_hammer_project( ...
            longitudeGrid,latitude*ones(size(longitudeGrid)));
        plot(ax,x,y,'-','Color',gridColor, ...
            'LineWidth',0.55,'HandleVisibility','off');
    end
    boundaryLatitude = linspace(-90,90,361);
    [xEast,yEast] = local_hammer_project( ...
        180*ones(size(boundaryLatitude)),boundaryLatitude);
    [xWest,yWest] = local_hammer_project( ...
        -180*ones(size(boundaryLatitude)),boundaryLatitude);
    plot(ax,[xEast fliplr(xWest) xEast(1)], ...
        [yEast fliplr(yWest) yEast(1)],'k-', ...
        'LineWidth',0.85,'HandleVisibility','off');
    axis(ax,'equal');
    xlim(ax,[-2.95 2.95]);
    ylim(ax,[-1.72 1.58]);
    set(ax,'XTick',[],'YTick',[],'Box','off');
    xlabel(ax,'Model-inertial longitude (center $0^\circ$, edge $\pm180^\circ$)');
    ylabel(ax,'Model-inertial latitude');
    if ~holdState
        hold(ax,'off');
    end
end


function [x,y] = local_hammer_project(longitudeDeg,latitudeDeg)
    longitudeDeg = mod(longitudeDeg+180,360)-180;
    lambda = deg2rad(longitudeDeg);
    phi = deg2rad(latitudeDeg);
    denominator = sqrt(1+cos(phi).*cos(lambda/2));
    x = 2*sqrt(2)*cos(phi).*sin(lambda/2)./denominator;
    y = sqrt(2)*sin(phi)./denominator;
end


function local_plot_hammer_marker(ax,direction,color,markerSize,marker)
    direction = direction(:)/norm(direction);
    longitude = atan2d(direction(2),direction(1));
    latitude = asind(min(1,max(-1,direction(3))));
    [x,y] = local_hammer_project(longitude,latitude);
    plot(ax,x,y,marker,'MarkerSize',markerSize, ...
        'MarkerFaceColor',color,'MarkerEdgeColor','k', ...
        'HandleVisibility','off');
end


function [colors,colorMap,periodLimits] = ...
        local_parent_period_colors(parentPeriods)
    parentPeriods = parentPeriods(:);
    colorMap = parula(256);
    periodMin = min(parentPeriods);
    periodMax = max(parentPeriods);
    if periodMax-periodMin <= 100*eps(max(1,max(abs(parentPeriods))))
        halfWidth = max(1,abs(periodMin))*1e-6;
        periodLimits = [periodMin-halfWidth periodMax+halfWidth];
        normalizedPeriod = 0.5*ones(size(parentPeriods));
    else
        periodLimits = [periodMin periodMax];
        normalizedPeriod = (parentPeriods-periodMin)/(periodMax-periodMin);
    end
    colorIndex = 1+round(normalizedPeriod*(size(colorMap,1)-1));
    colors = colorMap(colorIndex,:);
end


function local_plot_axis_marker(ax,direction,color,markerSize,marker)
    direction = direction(:)/norm(direction);
    longitude = mod(atan2d(direction(2),direction(1)),360);
    latitude = asind(min(1,max(-1,direction(3))));
    plot(ax,longitude,latitude,marker,'MarkerSize',markerSize, ...
        'MarkerFaceColor',color,'MarkerEdgeColor','k', ...
        'HandleVisibility','off');
end


function [axisDirection,concentration,medianOffsetDeg,rmsOffsetDeg] = ...
        local_axial_preferred_direction(directions)
    if isempty(directions) || size(directions,2) ~= 3 || ...
            any(~isfinite(directions(:)))
        error('FINESST:AllModesStudy:InvalidBestDirections', ...
            'Best-direction array must be finite and N-by-3.');
    end
    directions = directions./vecnorm(directions,2,2);
    orientationTensor = directions.'*directions/size(directions,1);
    [V,D] = eig((orientationTensor+orientationTensor.')/2,'vector');
    [concentration,iMax] = max(real(D));
    axisDirection = real(V(:,iMax));
    axisDirection = axisDirection/norm(axisDirection);
    if axisDirection(3) < 0 || ...
            (abs(axisDirection(3)) <= 10*eps && axisDirection(1) < 0)
        axisDirection = -axisDirection;
    end
    offsetDeg = acosd(min(1,max(0,abs(directions*axisDirection))));
    medianOffsetDeg = median(offsetDeg);
    rmsOffsetDeg = sqrt(mean(offsetDeg.^2));
end


function fig = local_plot_family_geometry( ...
        parentCaseIndex,normalSynodicPathAll, ...
        missionNormalSamplesAll,parentOrbitSynodicAll,parentPeriods, ...
        parentBestDirections,preferredAxis, ...
        concentration,medianOffsetDeg,familyName,missionDurationYears)
    nParents = numel(parentCaseIndex);
    [parentColors,parentColorMap,periodLimits] = ...
        local_parent_period_colors(parentPeriods);
    lineStyles = {'-','--',':','-.'};
    orientedBest = zeros(nParents,3);
    for jParent = 1:nParents
        orientedBest(jParent,:) = local_orient_direction( ...
            parentBestDirections(jParent,:),preferredAxis);
    end

    fig = figure('Color','w','Position',[30 100 1900 650]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    axOrbit = nexttile;
    hold(axOrbit,'on');
    for jParent = 1:nParents
        iCase = parentCaseIndex{jParent}(1);
        r = parentOrbitSynodicAll{iCase};
        plot3(axOrbit,r(1,:),r(2,:),r(3,:),'-', ...
            'Color',parentColors(jParent,:),'LineWidth',1.45, ...
            'HandleVisibility','off');
    end
    axis(axOrbit,'equal'); grid(axOrbit,'on'); box(axOrbit,'on');
    xlabel(axOrbit,'synodic $x$ [ND]');
    ylabel(axOrbit,'synodic $y$ [ND]');
    zlabel(axOrbit,'synodic $z$ [ND]');
    title(axOrbit,'Representative parent periodic orbits');
    view(axOrbit,38,24);
    colormap(axOrbit,parentColorMap);
    clim(axOrbit,periodLimits);
    cbPeriod = colorbar(axOrbit);
    ylabel(cbPeriod,'Parent period $T_p$ [TU]');

    axSynodic = nexttile;
    hold(axSynodic,'on');
    local_add_reference_sphere(axSynodic);
    seamAngles = zeros(0,1);
    for jParent = 1:nParents
        caseList = parentCaseIndex{jParent};
        for jMode = 1:numel(caseList)
            n = local_orient_path( ...
                normalSynodicPathAll{caseList(jMode)},preferredAxis);
            seamAngles(end+1,1) = acosd(min(1,max(-1, ...
                abs(dot(n(:,end),n(:,1)))))); %#ok<AGROW>
            nClosed = [n n(:,1)];
            plot3(axSynodic,nClosed(1,:),nClosed(2,:),nClosed(3,:), ...
                lineStyles{1+mod(jMode-1,numel(lineStyles))}, ...
                'Color',parentColors(jParent,:),'LineWidth',1.2, ...
                'HandleVisibility','off');
            plot3(axSynodic,n(1,1),n(2,1),n(3,1),'o', ...
                'MarkerSize',3,'MarkerFaceColor',parentColors(jParent,:), ...
                'MarkerEdgeColor','none','HandleVisibility','off');
        end
    end
    local_format_unit_sphere(axSynodic,'synodic');
    title(axSynodic,sprintf([ ...
        'Synodic one-$T_p$ normal lines, explicitly closed\n' ...
        'one axial sign; maximum sampled seam %.3g$^\circ$'], ...
        max(seamAngles)));

    axCoverage = nexttile;
    hold(axCoverage,'on');
    local_add_reference_sphere(axCoverage);
    for jParent = 1:nParents
        caseList = parentCaseIndex{jParent};
        for jMode = 1:numel(caseList)
            n = missionNormalSamplesAll{caseList(jMode)};
            sampleIndex = unique(round(linspace(1,size(n,2), ...
                min(400,size(n,2)))));
            n = local_orient_samples(n(:,sampleIndex),preferredAxis);
            scatter3(axCoverage,n(1,:),n(2,:),n(3,:),6, ...
                parentColors(jParent,:),'filled', ...
                'MarkerFaceAlpha',0.13,'MarkerEdgeColor','none', ...
                'HandleVisibility','off');
        end
        s = orientedBest(jParent,:);
        plot3(axCoverage,s(1),s(2),s(3),'o','MarkerSize',7, ...
            'MarkerFaceColor',parentColors(jParent,:), ...
            'MarkerEdgeColor','k','HandleVisibility','off');
    end
    plot3(axCoverage,preferredAxis(1),preferredAxis(2),preferredAxis(3), ...
        'p','MarkerSize',14,'MarkerFaceColor','k','MarkerEdgeColor','w', ...
        'HandleVisibility','off');
    local_format_unit_sphere(axCoverage,'model-inertial');
    title(axCoverage,sprintf([ ...
        '%.3g-year model-inertial coverage (unconnected samples)\n' ...
        'concentration %.3f; median parent offset %.1f$^\circ$'], ...
        missionDurationYears,concentration,medianOffsetDeg));

    sgtitle(sprintf([ ...
        '%s: parent geometry, closed synodic template, and mission coverage; ' ...
        'color = $T_p$'],strrep(familyName,'_','\_')));
end


function local_add_reference_sphere(ax)
    [xs,ys,zs] = sphere(48);
    surf(ax,xs,ys,zs,'FaceColor',[0.75 0.82 0.90], ...
        'FaceAlpha',0.07,'EdgeColor',[0.72 0.72 0.72], ...
        'EdgeAlpha',0.18,'HandleVisibility','off');
end


function local_format_unit_sphere(ax,frameName)
    axis(ax,'equal');
    xlim(ax,[-1 1]); ylim(ax,[-1 1]); zlim(ax,[-1 1]);
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,sprintf('%s $x$',frameName));
    ylabel(ax,sprintf('%s $y$',frameName));
    zlabel(ax,sprintf('%s $z$',frameName));
    view(ax,38,24);
end


function n = local_orient_samples(n,referenceAxis)
    referenceAxis = referenceAxis(:)/norm(referenceAxis);
    alignment = referenceAxis.'*n;
    signs = ones(1,size(n,2));
    signs(alignment < 0) = -1;
    n = n.*signs;
end


function n = local_orient_path(n,referenceAxis)
    referenceAxis = referenceAxis(:)/norm(referenceAxis);
    alignment = dot(mean(n,2),referenceAxis);
    if abs(alignment) <= 100*eps
        alignment = dot(n(:,1),referenceAxis);
    end
    if alignment < 0
        n = -n;
    end
end


function s = local_orient_direction(s,referenceAxis)
    s = s(:)/norm(s);
    referenceAxis = referenceAxis(:)/norm(referenceAxis);
    if dot(s,referenceAxis) < 0
        s = -s;
    end
    s = s.';
end


function local_export_figure(fig,pngPath,resolution,saveEditableFigures)
    exportgraphics(fig,pngPath,'Resolution',resolution);
    if saveEditableFigures
        [folder,name] = fileparts(pngPath);
        savefig(fig,fullfile(folder,[name '.fig']));
    end
end


function reason = local_main_rejection_reason(rejectedStatus)
    if isempty(rejectedStatus)
        reason = "none";
        return
    end
    [uniqueStatus,~,group] = unique(rejectedStatus);
    counts = accumarray(group,1);
    [~,iMax] = max(counts);
    reason = uniqueStatus(iMax);
end


function tag = local_duration_tag(durationYears)
    tag = sprintf('%.6gyr',durationYears);
    tag = strrep(tag,'.','p');
    tag = strrep(tag,'-','m');
end


function tf = local_cache_matches( ...
        analysis,ICfamily,Tfamily,floquetOpts,familyCfgCurrent)
    required = {'ICfamily','Tfamily','floquetOptions','configuration', ...
        'scanTable','modeTable','candidateParents','spectra', ...
        'analysisVersion','topology','branchDiagnostics'};
    tf = isstruct(analysis) && all(isfield(analysis,required));
    if ~tf
        return
    end

    tf = isequaln(analysis.ICfamily,ICfamily) && ...
        isequaln(analysis.Tfamily(:),Tfamily(:)) && ...
        isequaln(analysis.floquetOptions,floquetOpts) && ...
        strcmp(string(analysis.analysisVersion), ...
            string(familyCfgCurrent.analysisVersion)) && ...
        isequaln(analysis.configuration,familyCfgCurrent);
end
