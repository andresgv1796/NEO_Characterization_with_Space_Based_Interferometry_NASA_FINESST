%% RUN_JPL_Halo_Floquet_Spectrum_Survey
%
% FINESST / GMOS
% Whole-catalog Floquet survey of the JPL L1 northern-halo family.
%
% This script reads the CSV exactly in the form downloaded from the
% JPL / Three-Body Periodic Orbit database:
%
%   col  1 : Id
%   cols 2:7 : [x0 y0 z0 vx0 vy0 vz0]
%   col  8 : Jacobi constant
%   col  9 : Period [TU]
%   col 10 : Period [days]
%   col 11 : JPL stability index
%   col 12 : mass ratio mu
%
% For every periodic orbit in the catalog it:
%
%   1. computes the monodromy eigensystem using the existing
%      monodromy_CR3BP routine;
%   2. records all six Floquet multipliers;
%   3. uses analyze_floquet_spectrum for generic Floquet bookkeeping;
%   4. uses gmos_center_eigpair_from_eigs for the SAME canonical center-mode
%      selection used by the production GMOS pipeline;
%   5. classifies the non-trivial four-multiplier spectrum;
%   6. compares the computed spectral-radius stability index with the
%      stability index already supplied by the JPL catalog;
%   7. identifies the contiguous parts of the halo family for which the
%      current GMOS center-mode selector can find a complex unit-circle pair;
%   8. saves both orbit-level and multiplier-level tables.
%
% The catalog ROW is treated as the primary continuation coordinate.
% z0 is useful physically but is not guaranteed to be monotonic over the
% entire nonlinear halo family.
%
% Required existing / refactored functions:
%
%       monodromy_CR3BP
%       analyze_floquet_spectrum
%       gmos_center_eigpair_from_eigs
%
% monodromy_CR3BP performs the expensive 42-state integration exactly once
% per catalog orbit.  The two analysis helpers then reuse the returned V,D
% eigensystem without repeating that integration.
%
% -------------------------------------------------------------------------

clear;
clc;
close all;


%% ============================================================
%  0. PLOT DEFAULTS
%  ============================================================

set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

set(groot,'DefaultAxesFontSize',14);
set(groot,'DefaultTextFontSize',16);

set(groot,'DefaultFigureColor','w');
set(groot,'DefaultAxesColor','w');
set(groot,'DefaultAxesXColor','k');
set(groot,'DefaultAxesYColor','k');
set(groot,'DefaultTextColor','k');


%% ============================================================
%  1. CONFIGURATION
%  ============================================================

catalogFile = ...
    'L1_northern_halo_orbits_JPL_IC_Earth_Moon.csv';


% Match the CURRENT GMOS center-mode selector.
tolUnit = 1e-6;
tolImag = 1e-10;


% Mode character is now quantified by the vertical energy fraction
%
%   f_v = S_v^2/(S_p^2 + S_v^2).
%
% A center eigenvector is called vertical-dominant when f_v > 0.5.
% No arbitrary score threshold is used.
verticalDominanceThreshold = 0.50;


% Use the whole catalog.
rowRange = [];


% Parallel monodromy integrations can substantially reduce runtime.
useParallel = true;


% Rows from the morphology study that are useful to highlight.
highlightRows = [40 82 125 197 278 299];


% Save all numerical results.
saveResults = true;

outputDir = ...
    fullfile( ...
        pwd, ...
        'JPL_Halo_Floquet_Survey');

if ~exist(outputDir,'dir')
    mkdir(outputDir);
end


%% ============================================================
%  2. READ THE JPL CATALOG
%  ============================================================

opts = ...
    detectImportOptions( ...
        catalogFile, ...
        'VariableNamingRule','preserve');


catalog = ...
    readtable( ...
        catalogFile, ...
        opts);


if width(catalog) < 12

    error( ...
        ['Expected at least 12 columns in the JPL halo catalog. ', ...
         'Found %d columns.'], ...
        width(catalog));

end


% Read by COLUMN POSITION so the script is robust to spaces and punctuation
% in the database headers.

jplId = ...
    catalog{:,1};


X0Catalog = ...
    catalog{:,2:7};


jacobiCatalog = ...
    catalog{:,8};


periodTU = ...
    catalog{:,9};


periodDays = ...
    catalog{:,10};


jplStabilityIndex = ...
    catalog{:,11};


muCatalog = ...
    catalog{:,12};


nCatalog = ...
    height(catalog);


catalogRow = ...
    (1:nCatalog).';


if isempty(rowRange)

    rowsToRun = ...
        (1:nCatalog).';

else

    rowsToRun = ...
        rowRange(:);

end


if any(rowsToRun < 1 | rowsToRun > nCatalog)

    error('rowRange contains indices outside the catalog.');

end


fprintf('\n');
fprintf('============================================================\n');
fprintf('JPL L1 NORTHERN-HALO FLOQUET SURVEY\n');
fprintf('============================================================\n');

fprintf('Catalog                    : %s\n',catalogFile);
fprintf('Number of catalog rows     : %d\n',nCatalog);
fprintf('Rows analyzed              : %d\n',numel(rowsToRun));
fprintf('Unit-circle tolerance      : %.3e\n',tolUnit);
fprintf('Imaginary-part tolerance   : %.3e\n',tolImag);
fprintf('Vertical-dominance boundary: f_v > %.3f\n',verticalDominanceThreshold);

fprintf('\nCatalog columns used:\n');
for k = 1:12
    fprintf('  %2d : %s\n',k,catalog.Properties.VariableNames{k});
end


muSpread = ...
    max(muCatalog) - ...
    min(muCatalog);


fprintf('\nMass ratio:\n');
fprintf('  median mu = %.16e\n',median(muCatalog));
fprintf('  spread    = %.3e\n',muSpread);


%% ============================================================
%  3. STORAGE
%  ============================================================

success = false(nCatalog,1);
failureMessage = cell(nCatalog,1);

lambdaAll = complex(nan(nCatalog,6),nan(nCatalog,6));

planarScoreAll = nan(nCatalog,6);
verticalScoreAll = nan(nCatalog,6);

spectralRadius = nan(nCatalog,1);
computedStabilityIndex = nan(nCatalog,1);

determinantError = nan(nCatalog,1);
reciprocalPairError = nan(nCatalog,1);
trivialPairDistance = nan(nCatalog,1);

nPositiveUnitComplexPairs = nan(nCatalog,1);

spectralClass = strings(nCatalog,1);
spectralClassCode = nan(nCatalog,1);

dominantUnstableLambda = complex(nan(nCatalog,1),nan(nCatalog,1));
dominantStableLambda = complex(nan(nCatalog,1),nan(nCatalog,1));

gmosWouldFindCenter = false(nCatalog,1);
verticalModeDominant = false(nCatalog,1);

selectedLambda = complex(nan(nCatalog,1),nan(nCatalog,1));
selectedRho_rad = nan(nCatalog,1);
selectedRho_cycles = nan(nCatalog,1);

selectedVerticalScore = nan(nCatalog,1);
selectedPlanarScore = nan(nCatalog,1);
selectedVerticalFraction = nan(nCatalog,1);
selectedPositionVerticalFraction = nan(nCatalog,1);
selectedVelocityVerticalFraction = nan(nCatalog,1);


% One set of tolerances shared by the generic Floquet analysis and the
% canonical GMOS center selector.

floquetOpts = struct();
floquetOpts.tolUnit = tolUnit;
floquetOpts.tolImag = tolImag;

centerOpts = floquetOpts;

% Do not issue hundreds of mode-character warnings inside the catalog sweep.
% The continuous vertical fraction and dominance flag are retained directly.
centerOpts.warnMismatch = false;


%% ============================================================
%  4. COMPUTE + ANALYZE ONE MONODROMY EIGENSYSTEM PER ORBIT
%  ============================================================

fprintf('\n');
fprintf('Computing monodromy eigensystems and Floquet diagnostics...\n');


runParallel = false;

if useParallel

    try

        if license('test','Distrib_Computing_Toolbox')

            pool = gcp('nocreate');

            if isempty(pool)
                parpool;
            end

            runParallel = true;

        end

    catch ME

        warning( ...
            'Parallel setup failed. Falling back to serial loop: %s', ...
            ME.message);

    end

end


nRun = numel(rowsToRun);


% PARFOR writes only into arrays directly sliced by kk.  Results are mapped
% back to the original catalog rows after the loop.

lambdaRun = complex(nan(nRun,6),nan(nRun,6));

planarScoreRun = nan(nRun,6);
verticalScoreRun = nan(nRun,6);

spectralRadiusRun = nan(nRun,1);
stabilityIndexRun = nan(nRun,1);

determinantErrorRun = nan(nRun,1);
reciprocalPairErrorRun = nan(nRun,1);
trivialPairDistanceRun = nan(nRun,1);

nCenterPairsRun = nan(nRun,1);

spectralClassRun = strings(nRun,1);
spectralClassCodeRun = nan(nRun,1);

unstableLambdaRun = complex(nan(nRun,1),nan(nRun,1));
stableLambdaRun = complex(nan(nRun,1),nan(nRun,1));

centerFoundRun = false(nRun,1);
verticalModeDominantRun = false(nRun,1);

selectedLambdaRun = complex(nan(nRun,1),nan(nRun,1));
selectedRhoRun = nan(nRun,1);
selectedRhoCyclesRun = nan(nRun,1);

selectedVerticalScoreRun = nan(nRun,1);
selectedPlanarScoreRun = nan(nRun,1);
selectedVerticalFractionRun = nan(nRun,1);
selectedPositionVerticalFractionRun = nan(nRun,1);
selectedVelocityVerticalFractionRun = nan(nRun,1);

successRun = false(nRun,1);
failureMessageRun = cell(nRun,1);


tic;


if runParallel

    parfor kk = 1:nRun

        i = rowsToRun(kk);

        try

            x0 = X0Catalog(i,:).';
            T = periodTU(i);
            mu = muCatalog(i);

            % EXPENSIVE operation: exactly one monodromy integration.
            M_eigs = monodromy_CR3BP([0,T],x0,mu);

            V = M_eigs(:,1:6);
            D = M_eigs(:,7:12);

            % Generic spectral analysis.
            spec = analyze_floquet_spectrum(V,D,floquetOpts);

            lambdaRun(kk,:) = spec.lambda.';
            planarScoreRun(kk,:) = spec.planarScore.';
            verticalScoreRun(kk,:) = spec.verticalScore.';

            spectralRadiusRun(kk) = spec.spectralRadius;
            stabilityIndexRun(kk) = spec.stabilityIndex;

            determinantErrorRun(kk) = spec.determinantError;
            reciprocalPairErrorRun(kk) = spec.reciprocalPairError;
            trivialPairDistanceRun(kk) = spec.trivialPairDistance;

            nCenterPairsRun(kk) = spec.nPositiveUnitComplexPairs;

            spectralClassRun(kk) = spec.spectralClass;
            spectralClassCodeRun(kk) = spec.spectralClassCode;

            unstableLambdaRun(kk) = spec.dominantUnstableLambda;
            stableLambdaRun(kk) = spec.dominantStableLambda;


            % Canonical production GMOS center-mode selection.
            %
            % This uses the SAME helper as
            % gmos_center_eigpair_from_monodromy, but reuses the already
            % computed V,D so the 42-state integration is NOT repeated.
            try

                cen = gmos_center_eigpair_from_eigs( ...
                    V, ...
                    D, ...
                    'vertical', ...
                    centerOpts);

                centerFoundRun(kk) = true;

                selectedLambdaRun(kk) = cen.lambda;
                selectedRhoRun(kk) = cen.rho0;
                selectedRhoCyclesRun(kk) = cen.rho0_cycles;

                selectedVerticalScoreRun(kk) = cen.proj.vertical;
                selectedPlanarScoreRun(kk) = cen.proj.planar;
                selectedVerticalFractionRun(kk) = cen.proj.verticalFraction;

                selectedPositionVerticalFractionRun(kk) = ...
                    cen.proj.positionVerticalFraction;

                selectedVelocityVerticalFractionRun(kk) = ...
                    cen.proj.velocityVerticalFraction;

                verticalModeDominantRun(kk) = ...
                    cen.proj.verticalDominant;

            catch

                centerFoundRun(kk) = false;

            end


            successRun(kk) = true;


        catch ME

            successRun(kk) = false;
            failureMessageRun{kk} = ME.message;

        end

    end


else

    for kk = 1:nRun

        i = rowsToRun(kk);

        if mod(kk,50) == 1 || kk == nRun
            fprintf('  row %d / %d\n',kk,nRun);
        end

        try

            x0 = X0Catalog(i,:).';
            T = periodTU(i);
            mu = muCatalog(i);

            % EXPENSIVE operation: exactly one monodromy integration.
            M_eigs = monodromy_CR3BP([0,T],x0,mu);

            V = M_eigs(:,1:6);
            D = M_eigs(:,7:12);

            spec = analyze_floquet_spectrum(V,D,floquetOpts);

            lambdaRun(kk,:) = spec.lambda.';
            planarScoreRun(kk,:) = spec.planarScore.';
            verticalScoreRun(kk,:) = spec.verticalScore.';

            spectralRadiusRun(kk) = spec.spectralRadius;
            stabilityIndexRun(kk) = spec.stabilityIndex;

            determinantErrorRun(kk) = spec.determinantError;
            reciprocalPairErrorRun(kk) = spec.reciprocalPairError;
            trivialPairDistanceRun(kk) = spec.trivialPairDistance;

            nCenterPairsRun(kk) = spec.nPositiveUnitComplexPairs;

            spectralClassRun(kk) = spec.spectralClass;
            spectralClassCodeRun(kk) = spec.spectralClassCode;

            unstableLambdaRun(kk) = spec.dominantUnstableLambda;
            stableLambdaRun(kk) = spec.dominantStableLambda;


            try

                cen = gmos_center_eigpair_from_eigs( ...
                    V, ...
                    D, ...
                    'vertical', ...
                    centerOpts);

                centerFoundRun(kk) = true;

                selectedLambdaRun(kk) = cen.lambda;
                selectedRhoRun(kk) = cen.rho0;
                selectedRhoCyclesRun(kk) = cen.rho0_cycles;

                selectedVerticalScoreRun(kk) = cen.proj.vertical;
                selectedPlanarScoreRun(kk) = cen.proj.planar;
                selectedVerticalFractionRun(kk) = cen.proj.verticalFraction;

                selectedPositionVerticalFractionRun(kk) = ...
                    cen.proj.positionVerticalFraction;

                selectedVelocityVerticalFractionRun(kk) = ...
                    cen.proj.velocityVerticalFraction;

                verticalModeDominantRun(kk) = ...
                    cen.proj.verticalDominant;

            catch

                centerFoundRun(kk) = false;

            end


            successRun(kk) = true;


        catch ME

            successRun(kk) = false;
            failureMessageRun{kk} = ME.message;

        end

    end

end


runtime_s = toc;


%% ============================================================
%  5. MAP THE WORKER RESULTS BACK TO CATALOG ROWS
%  ============================================================

lambdaAll(rowsToRun,:) = lambdaRun;

planarScoreAll(rowsToRun,:) = planarScoreRun;
verticalScoreAll(rowsToRun,:) = verticalScoreRun;

spectralRadius(rowsToRun) = spectralRadiusRun;
computedStabilityIndex(rowsToRun) = stabilityIndexRun;

determinantError(rowsToRun) = determinantErrorRun;
reciprocalPairError(rowsToRun) = reciprocalPairErrorRun;
trivialPairDistance(rowsToRun) = trivialPairDistanceRun;

nPositiveUnitComplexPairs(rowsToRun) = nCenterPairsRun;

spectralClass(rowsToRun) = spectralClassRun;
spectralClassCode(rowsToRun) = spectralClassCodeRun;

dominantUnstableLambda(rowsToRun) = unstableLambdaRun;
dominantStableLambda(rowsToRun) = stableLambdaRun;

gmosWouldFindCenter(rowsToRun) = centerFoundRun;
verticalModeDominant(rowsToRun) = verticalModeDominantRun;

selectedLambda(rowsToRun) = selectedLambdaRun;
selectedRho_rad(rowsToRun) = selectedRhoRun;
selectedRho_cycles(rowsToRun) = selectedRhoCyclesRun;

selectedVerticalScore(rowsToRun) = selectedVerticalScoreRun;
selectedPlanarScore(rowsToRun) = selectedPlanarScoreRun;
selectedVerticalFraction(rowsToRun) = selectedVerticalFractionRun;

selectedPositionVerticalFraction(rowsToRun) = ...
    selectedPositionVerticalFractionRun;

selectedVelocityVerticalFraction(rowsToRun) = ...
    selectedVelocityVerticalFractionRun;

success(rowsToRun) = successRun;
failureMessage(rowsToRun) = failureMessageRun;


stabilityIndexAbsError = abs( ...
    computedStabilityIndex - ...
    jplStabilityIndex);

stabilityIndexRelError = ...
    stabilityIndexAbsError ./ ...
    max(1,abs(jplStabilityIndex));


fprintf('\nMonodromy survey runtime = %.3f s\n',runtime_s);
fprintf('Successful rows          = %d / %d\n', ...
    nnz(success(rowsToRun)),numel(rowsToRun));

fprintf('Genuine nontrivial center candidates = %d\n', ...
    nnz(gmosWouldFindCenter(rowsToRun)));


%% ============================================================
%  6. BUILD ORBIT-LEVEL SUMMARY TABLE
%  ============================================================

summaryTable = ...
    table( ...
        catalogRow, ...
        jplId, ...
        X0Catalog(:,1), ...
        X0Catalog(:,2), ...
        X0Catalog(:,3), ...
        X0Catalog(:,4), ...
        X0Catalog(:,5), ...
        X0Catalog(:,6), ...
        jacobiCatalog, ...
        periodTU, ...
        periodDays, ...
        jplStabilityIndex, ...
        muCatalog, ...
        success, ...
        spectralRadius, ...
        computedStabilityIndex, ...
        stabilityIndexAbsError, ...
        stabilityIndexRelError, ...
        determinantError, ...
        reciprocalPairError, ...
        trivialPairDistance, ...
        nPositiveUnitComplexPairs, ...
        gmosWouldFindCenter, ...
        verticalModeDominant, ...
        real(selectedLambda), ...
        imag(selectedLambda), ...
        abs(selectedLambda), ...
        selectedRho_rad, ...
        selectedRho_cycles, ...
        selectedVerticalScore, ...
        selectedPlanarScore, ...
        selectedVerticalFraction, ...
        selectedPositionVerticalFraction, ...
        selectedVelocityVerticalFraction, ...
        spectralClass, ...
        spectralClassCode, ...
        real(dominantUnstableLambda), ...
        imag(dominantUnstableLambda), ...
        abs(dominantUnstableLambda), ...
        real(dominantStableLambda), ...
        imag(dominantStableLambda), ...
        abs(dominantStableLambda), ...
        'VariableNames',{ ...
            'CatalogRow', ...
            'JPL_Id', ...
            'x0', ...
            'y0', ...
            'z0', ...
            'vx0', ...
            'vy0', ...
            'vz0', ...
            'JacobiConstant', ...
            'PeriodTU', ...
            'PeriodDays', ...
            'JPL_StabilityIndex', ...
            'MassRatio', ...
            'MonodromySuccess', ...
            'SpectralRadius', ...
            'ComputedStabilityIndex', ...
            'StabilityIndexAbsError', ...
            'StabilityIndexRelError', ...
            'DeterminantError', ...
            'ReciprocalPairError', ...
            'TrivialPairDistance', ...
            'NumPositiveUnitComplexPairs', ...
            'GMOSWouldFindCenter', ...
            'VerticalModeDominant', ...
            'SelectedLambdaReal', ...
            'SelectedLambdaImag', ...
            'SelectedLambdaAbs', ...
            'SelectedRho_rad', ...
            'SelectedRho_cycles', ...
            'SelectedVerticalScore', ...
            'SelectedPlanarScore', ...
            'SelectedVerticalFraction', ...
            'SelectedPositionVerticalFraction', ...
            'SelectedVelocityVerticalFraction', ...
            'SpectralClass', ...
            'SpectralClassCode', ...
            'DominantUnstableLambdaReal', ...
            'DominantUnstableLambdaImag', ...
            'DominantUnstableLambdaAbs', ...
            'DominantStableLambdaReal', ...
            'DominantStableLambdaImag', ...
            'DominantStableLambdaAbs'});


%% ============================================================
%  7. BUILD MULTIPLIER-LEVEL LONG TABLE
%  ============================================================

nLong = ...
    6*nCatalog;


longCatalogRow = ...
    repelem(catalogRow,6);


longMultiplierIndex = ...
    repmat((1:6).',nCatalog,1);


lamVector = ...
    reshape(lambdaAll.',[],1);


planarVector = ...
    reshape(planarScoreAll.',[],1);


verticalVector = ...
    reshape(verticalScoreAll.',[],1);


multiplierTable = ...
    table( ...
        longCatalogRow, ...
        longMultiplierIndex, ...
        real(lamVector), ...
        imag(lamVector), ...
        abs(lamVector), ...
        angle(lamVector), ...
        planarVector, ...
        verticalVector, ...
        'VariableNames',{ ...
            'CatalogRow', ...
            'MultiplierIndex', ...
            'LambdaReal', ...
            'LambdaImag', ...
            'LambdaAbs', ...
            'LambdaAngle_rad', ...
            'PlanarProjectionScore', ...
            'VerticalProjectionScore'});


%% ============================================================
%  8. CONTIGUOUS CATALOG-ROW RUNS (INVENTORY ONLY)
%  ============================================================

candidateRowRuns = ...
    boolean_row_runs( ...
        gmosWouldFindCenter & success);


verticalDominantRowRuns = ...
    boolean_row_runs( ...
        verticalModeDominant & success);


candidateRowRunTable = ...
    row_run_table( ...
        candidateRowRuns, ...
        summaryTable, ...
        'CurrentGMOSCandidate');


verticalDominantRowRunTable = ...
    row_run_table( ...
        verticalDominantRowRuns, ...
        summaryTable, ...
        'VerticalDominantCandidate');


fprintf('\n');
fprintf('============================================================\n');
fprintf('WHOLE-FAMILY FLOQUET SUMMARY\n');
fprintf('============================================================\n');

fprintf([ ...
    '\nIMPORTANT CATALOG-ORDER NOTE:\n', ...
    '  Contiguous catalog-row runs below are inventory summaries only.\n', ...
    '  They are NOT continuous dynamical-family intervals because the\n', ...
    '  JPL table interleaves distinct halo branches in part of the catalog.\n\n']);

fprintf('Successful monodromy rows                : %d / %d\n', ...
    nnz(success),nCatalog);

fprintf('Rows with >=1 current GMOS center pair   : %d\n', ...
    nnz(gmosWouldFindCenter));

fprintf('Rows with vertical-dominant center mode   : %d\n', ...
    nnz(verticalModeDominant));

fprintf('Rows classified saddle-center            : %d\n', ...
    nnz(spectralClass=="saddle-center"));

fprintf('Rows classified elliptic-elliptic        : %d\n', ...
    nnz(spectralClass=="elliptic-elliptic"));

fprintf('Rows classified saddle-saddle            : %d\n', ...
    nnz(spectralClass=="saddle-saddle"));

fprintf('Rows classified complex quartet          : %d\n', ...
    nnz(spectralClass=="complex-quartet"));

fprintf('Rows classified transition/degenerate    : %d\n', ...
    nnz(spectralClass=="transition/degenerate"));


finiteStabilityError = ...
    stabilityIndexAbsError(isfinite(stabilityIndexAbsError));


if ~isempty(finiteStabilityError)

    fprintf('\nJPL stability-index validation:\n');

    fprintf('  max absolute error = %.3e\n', ...
        max(finiteStabilityError));

    fprintf('  RMS absolute error = %.3e\n', ...
        sqrt(mean(finiteStabilityError.^2)));

end


fprintf('\nCurrent GMOS center-candidate catalog-row runs:\n');
disp(candidateRowRunTable);


fprintf('\nVertical-dominant candidate catalog-row runs:\n');
disp(verticalDominantRowRunTable);


%% ============================================================
%  9. FIGURE: CATALOG + STABILITY CROSS-CHECK
%  ============================================================

figCatalog = ...
    figure( ...
        'Name','JPL halo catalog and Floquet stability', ...
        'Color','w', ...
        'Position',[60 60 1500 900]);


tl = ...
    tiledlayout( ...
        figCatalog, ...
        3,2, ...
        'TileSpacing','compact', ...
        'Padding','compact');


%% z0

ax = ...
    nexttile(tl);

plot(ax,catalogRow,X0Catalog(:,3),'k-');
grid(ax,'on');

xlabel(ax,'JPL catalog row');
ylabel(ax,'$z_0$ [LU]');
title(ax,'Halo-family geometry');


%% Jacobi constant

ax = ...
    nexttile(tl);

plot(ax,catalogRow,jacobiCatalog,'k-');
grid(ax,'on');

xlabel(ax,'JPL catalog row');
ylabel(ax,'$C$');
title(ax,'Jacobi constant');


%% Period

ax = ...
    nexttile(tl);

plot(ax,catalogRow,periodDays,'k-');
grid(ax,'on');

xlabel(ax,'JPL catalog row');
ylabel(ax,'Period [days]');
title(ax,'Parent-orbit period');


%% JPL stability index

ax = ...
    nexttile(tl);

semilogy(ax,catalogRow,jplStabilityIndex,'k-');
grid(ax,'on');

xlabel(ax,'JPL catalog row');
ylabel(ax,'JPL stability index');
title(ax,'Database stability index');

yline(ax,1,'k--','$s=1$');


%% Computed spectral radius

ax = ...
    nexttile(tl);

semilogy(ax,catalogRow,spectralRadius,'k-');
grid(ax,'on');

xlabel(ax,'JPL catalog row');
ylabel(ax,'$\rho(M)=\max_i|\lambda_i|$');
title(ax,'Computed monodromy spectral radius');

yline(ax,1,'k--','$|\lambda|=1$');


%% Stability-index cross-check

ax = ...
    nexttile(tl);

semilogy( ...
    ax, ...
    catalogRow, ...
    max(stabilityIndexAbsError,eps), ...
    'k-');

grid(ax,'on');

xlabel(ax,'JPL catalog row');
ylabel(ax,'$|s_{\rm calc}-s_{\rm JPL}|$');
title(ax,'Independent stability-index cross-check');


title( ...
    tl, ...
    'JPL L1 northern-halo catalog: geometry and Floquet stability');


add_highlight_rows_to_figure(figCatalog,highlightRows);


%% ============================================================
%  10. FIGURE: ALL SIX FLOQUET-MULTIPLIER MODULI
%  ============================================================

figModuli = ...
    figure( ...
        'Name','Whole-family Floquet multiplier moduli', ...
        'Color','w', ...
        'Position',[80 80 1450 850]);


tl = ...
    tiledlayout( ...
        figModuli, ...
        3,1, ...
        'TileSpacing','compact', ...
        'Padding','compact');


%% All six multiplier moduli

ax = ...
    nexttile(tl);


hold(ax,'on');
grid(ax,'on');


for j = 1:6

    semilogy( ...
        ax, ...
        catalogRow, ...
        abs(lambdaAll(:,j)), ...
        '.', ...
        'MarkerSize',4);

end


yline(ax,1,'k--');

xlabel(ax,'JPL catalog row');
ylabel(ax,'$|\lambda|$');

title(ax,'All six Floquet-multiplier moduli');


%% Selected GMOS vertical-center rotation number

ax = ...
    nexttile(tl);


scatter( ...
    ax, ...
    catalogRow(gmosWouldFindCenter), ...
    selectedRho_cycles(gmosWouldFindCenter), ...
    18, ...
    selectedVerticalFraction(gmosWouldFindCenter), ...
    'filled');


grid(ax,'on');

xlabel(ax,'JPL catalog row');
ylabel(ax,'$\rho_0/(2\pi)$ [cycles]');

title(ax,'Center rotation number selected by current vertical-mode rule');

cb = ...
    colorbar(ax);

ylabel(cb,'vertical fraction of selected eigenvector');


%% Vertical and planar mode projections

ax = ...
    nexttile(tl);


hold(ax,'on');
grid(ax,'on');


plot( ...
    ax, ...
    catalogRow(gmosWouldFindCenter), ...
    selectedVerticalScore(gmosWouldFindCenter), ...
    '.', ...
    'DisplayName','$S_v$');


plot( ...
    ax, ...
    catalogRow(gmosWouldFindCenter), ...
    selectedPlanarScore(gmosWouldFindCenter), ...
    '.', ...
    'DisplayName','$S_p$');


xlabel(ax,'JPL catalog row');
ylabel(ax,'eigenvector projection norm');

title(ax,'Mode character of selected unit-circle pair');

legend(ax,'Location','best');


title( ...
    tl, ...
    'Floquet structure relevant to GMOS QPO initialization');


add_highlight_rows_to_figure(figModuli,highlightRows);


%% ============================================================
%  10A. FIGURE: VERTICAL CHARACTER OF THE SELECTED CENTER MODE
%  ============================================================

figModeFractions = figure( ...
    'Name','Selected center-mode vertical fractions', ...
    'Color','w', ...
    'Position',[110 110 1450 700]);

tlFrac = tiledlayout(figModeFractions,3,1, ...
    'TileSpacing','compact','Padding','compact');

fractionData = { ...
    selectedVerticalFraction, ...
    selectedPositionVerticalFraction, ...
    selectedVelocityVerticalFraction};

fractionTitles = { ...
    'Full state-space vertical fraction $f_{v,\mathrm{state}}$', ...
    'Position-only vertical fraction $f_{v,\mathrm{pos}}$', ...
    'Velocity-only vertical fraction $f_{v,\mathrm{vel}}$'};

for kFrac = 1:3

    ax = nexttile(tlFrac);

    scatter(ax, ...
        catalogRow(gmosWouldFindCenter), ...
        fractionData{kFrac}(gmosWouldFindCenter), ...
        16, ...
        selectedRho_cycles(gmosWouldFindCenter), ...
        'filled');

    hold(ax,'on');
    grid(ax,'on');

    yline(ax,0.5,'k--','$0.5$');

    ylim(ax,[0 1]);

    xlabel(ax,'JPL catalog row');
    ylabel(ax,'vertical fraction');
    title(ax,fractionTitles{kFrac});

    cb = colorbar(ax);
    ylabel(cb,'$\rho_0/(2\pi)$ [cycles]');

end

title(tlFrac, ...
    ['Selected Floquet-center mode character in JPL catalog order ', ...
     '(inventory only)']);

add_highlight_rows_to_figure(figModeFractions,highlightRows);


%% ============================================================
%  11. FIGURE: SPECTRAL TOPOLOGY + GMOS ADMISSIBILITY
%  ============================================================

figTopology = ...
    figure( ...
        'Name','Floquet topology and GMOS candidate map', ...
        'Color','w', ...
        'Position',[100 100 1450 760]);


tl = ...
    tiledlayout( ...
        figTopology, ...
        3,1, ...
        'TileSpacing','compact', ...
        'Padding','compact');


%% Topology code

ax = ...
    nexttile(tl);


stairs( ...
    ax, ...
    catalogRow, ...
    spectralClassCode, ...
    'k-', ...
    'LineWidth',1);


grid(ax,'on');

xlabel(ax,'JPL catalog row');
ylabel(ax,'spectral topology');

yticks(ax,0:5);

yticklabels( ...
    ax,{ ...
        'failure', ...
        'saddle-center', ...
        'elliptic-elliptic', ...
        'saddle-saddle', ...
        'complex quartet', ...
        'transition'});

ylim(ax,[-0.5 5.5]);

title(ax,'Non-trivial four-multiplier topology');


%% Exact current GMOS center-candidate availability

ax = ...
    nexttile(tl);


stairs( ...
    ax, ...
    catalogRow, ...
    double(gmosWouldFindCenter), ...
    'k-', ...
    'LineWidth',1.2);


hold(ax,'on');


stairs( ...
    ax, ...
    catalogRow, ...
    0.85*double(verticalModeDominant), ...
    'LineWidth',1.2);


grid(ax,'on');

xlabel(ax,'JPL catalog row');
ylabel(ax,'availability');

yticks(ax,[0 0.85 1]);

yticklabels(ax,{ ...
    'none', ...
    'vertical dominant', ...
    'unit-complex pair'});

ylim(ax,[-0.1 1.1]);

title(ax,['Center-mode availability in JPL catalog order ', ...
    '(inventory, not branch continuation)']);


%% Number of unit-complex pairs

ax = ...
    nexttile(tl);


stairs( ...
    ax, ...
    catalogRow, ...
    nPositiveUnitComplexPairs, ...
    'k-', ...
    'LineWidth',1.1);


grid(ax,'on');

xlabel(ax,'JPL catalog row');
ylabel(ax,'number of positive-imaginary unit pairs');

yticks(ax,0:2);

title(ax,'Number of complex conjugate pairs on the unit circle');


title( ...
    tl, ...
    ['Floquet-spectrum inventory in JPL catalog order ', ...
     '(catalog row is not a continuation coordinate)']);


add_highlight_rows_to_figure(figTopology,highlightRows);


%% ============================================================
%  12. FIGURE: FLOQUET MULTIPLIERS NEAR THE UNIT CIRCLE
%  ============================================================

figPlane = ...
    figure( ...
        'Name','Floquet multipliers in complex plane', ...
        'Color','w', ...
        'Position',[120 120 900 800]);


ax = ...
    axes(figPlane);


theta = ...
    linspace(0,2*pi,600);


plot( ...
    ax, ...
    cos(theta), ...
    sin(theta), ...
    'k--', ...
    'LineWidth',1);


hold(ax,'on');
grid(ax,'on');
axis(ax,'equal');


rowColor = ...
    repelem(catalogRow,6);


lambdaVector = ...
    reshape(lambdaAll.',[],1);


planeMask = ...
    isfinite(real(lambdaVector)) & ...
    abs(lambdaVector) <= 3;


scatter( ...
    ax, ...
    real(lambdaVector(planeMask)), ...
    imag(lambdaVector(planeMask)), ...
    9, ...
    rowColor(planeMask), ...
    'filled');


xlabel(ax,'$\mathrm{Re}(\lambda)$');
ylabel(ax,'$\mathrm{Im}(\lambda)$');

title(ax,'Floquet multipliers with $|\lambda|\leq3$');

xlim(ax,[-3 3]);
ylim(ax,[-3 3]);

cb = ...
    colorbar(ax);

ylabel(cb,'JPL catalog row');


%% ============================================================
%  13. FIGURE: SPECTRA AT REPRESENTATIVE MORPHOLOGY-SURVEY ROWS
%  ============================================================

validHighlightRows = ...
    highlightRows( ...
        highlightRows >= 1 & ...
        highlightRows <= nCatalog);


nHighlight = ...
    numel(validHighlightRows);


if nHighlight > 0

    figExamples = ...
        figure( ...
            'Name','Representative Floquet spectra', ...
            'Color','w', ...
            'Position',[80 80 1500 820]);


    nCols = ...
        min(3,nHighlight);


    nRows = ...
        ceil(nHighlight/nCols);


    tl = ...
        tiledlayout( ...
            figExamples, ...
            nRows,nCols, ...
            'TileSpacing','compact', ...
            'Padding','compact');


    for k = 1:nHighlight

        i = ...
            validHighlightRows(k);


        ax = ...
            nexttile(tl);


        plot( ...
            ax, ...
            cos(theta), ...
            sin(theta), ...
            'k--');


        hold(ax,'on');
        grid(ax,'on');
        axis(ax,'equal');


        lam = ...
            lambdaAll(i,:);


        scatter( ...
            ax, ...
            real(lam), ...
            imag(lam), ...
            50, ...
            abs(lam), ...
            'filled');


        xlabel(ax,'$\mathrm{Re}(\lambda)$');
        ylabel(ax,'$\mathrm{Im}(\lambda)$');


        title( ...
            ax, ...
            sprintf( ...
                ['row %d, $z_0=%.3f$\n', ...
                 '%s'], ...
                i, ...
                X0Catalog(i,3), ...
                spectralClass(i)));


        % Keep enough room for moderately unstable multipliers without
        % letting a very large hyperbolic multiplier destroy the unit-circle
        % view.  Large values remain available numerically in the tables.
        xlim(ax,[-4 4]);
        ylim(ax,[-4 4]);

    end


    title( ...
        tl, ...
        'Representative parent-orbit Floquet spectra');

end


%% ============================================================
%  14. FIGURE: JPL STABILITY INDEX VERSUS COMPUTED FORMULA
%  ============================================================

figValidation = ...
    figure( ...
        'Name','JPL stability-index validation', ...
        'Color','w', ...
        'Position',[140 140 700 620]);


valid = ...
    success & ...
    isfinite(jplStabilityIndex) & ...
    isfinite(computedStabilityIndex);


loglog( ...
    jplStabilityIndex(valid), ...
    computedStabilityIndex(valid), ...
    'k.');


hold on;
grid on;
axis equal;


sMin = ...
    min([ ...
        jplStabilityIndex(valid); ...
        computedStabilityIndex(valid)]);


sMax = ...
    max([ ...
        jplStabilityIndex(valid); ...
        computedStabilityIndex(valid)]);


plot( ...
    [sMin sMax], ...
    [sMin sMax], ...
    'k--');


xlabel('JPL stability index');
ylabel('Computed $\frac{1}{2}(\rho+\rho^{-1})$');

title('Validation of the database stability index');


%% ============================================================
%  15. SAVE
%  ============================================================

if saveResults

    summaryCsv = ...
        fullfile( ...
            outputDir, ...
            'JPL_Halo_Floquet_Summary.csv');


    multiplierCsv = ...
        fullfile( ...
            outputDir, ...
            'JPL_Halo_Floquet_Multipliers.csv');


    intervalCsv = ...
        fullfile( ...
            outputDir, ...
            'JPL_Halo_GMOS_Candidate_CatalogRowRuns.csv');


    verticalIntervalCsv = ...
        fullfile( ...
            outputDir, ...
            'JPL_Halo_VerticalDominant_CatalogRowRuns.csv');


    matFile = ...
        fullfile( ...
            outputDir, ...
            'JPL_Halo_Floquet_Survey.mat');


    writetable(summaryTable,summaryCsv);

    writetable(multiplierTable,multiplierCsv);

    writetable(candidateRowRunTable,intervalCsv);

    writetable(verticalDominantRowRunTable,verticalIntervalCsv);


    save( ...
        matFile, ...
        'catalogFile', ...
        'catalog', ...
        'summaryTable', ...
        'multiplierTable', ...
        'candidateRowRunTable', ...
        'verticalDominantRowRunTable', ...
        'selectedVerticalFraction', ...
        'selectedPositionVerticalFraction', ...
        'selectedVelocityVerticalFraction', ...
        'lambdaAll', ...
        'planarScoreAll', ...
        'verticalScoreAll', ...
        'tolUnit', ...
        'tolImag', ...
        'verticalDominanceThreshold', ...
        'runtime_s', ...
        '-v7.3');


    fprintf('\nSaved:\n');
    fprintf('  %s\n',summaryCsv);
    fprintf('  %s\n',multiplierCsv);
    fprintf('  %s\n',intervalCsv);
    fprintf('  %s\n',verticalIntervalCsv);
    fprintf('  %s\n',matFile);

end


%% ============================================================
%  16. SCIENTIFIC INTERPRETATION NOTES
%  ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('INTERPRETATION\n');
fprintf('============================================================\n');

fprintf([ ...
    'GMOSWouldFindCenter is produced by the canonical function\n', ...
    '  gmos_center_eigpair_from_eigs,\n', ...
    'after analyze_floquet_spectrum removes the two trivial +1 multipliers.\n', ...
    'The nontrivial unit-circle tolerances are %.1e and %.1e.\n'], ...
    tolUnit,tolImag);

fprintf([ ...
    '\nSelectedRho_rad is the argument of the genuine nontrivial unit-circle\n', ...
    'candidate with the largest vertical eigenvector projection, matching\n', ...
    'the production mode="vertical" selection logic.\n']);

fprintf([ ...
    '\nVerticalModeDominant means\n', ...
    '  f_v = S_v^2/(S_p^2+S_v^2) > 0.5.\n', ...
    'This replaces the previous arbitrary score > 0.2 interpretation.\n']);

fprintf([ ...
    '\nThe JPL stability index should agree with\n', ...
    '  s = 0.5*(rho(M) + 1/rho(M)),\n', ...
    'where rho(M) is the monodromy spectral radius.\n']);

fprintf([ ...
    '\nCatalogRow is used only as a database inventory coordinate.\n', ...
    'Do NOT interpret adjacent rows as neighboring points on one continuous\n', ...
    'halo branch in the interleaved portion of the JPL catalog.\n', ...
    'Physical bifurcation claims require explicit branch reconstruction.\n']);

fprintf([ ...
    '\nThree vertical-character measures are retained:\n', ...
    '  f_v,state : full 6-state vertical fraction,\n', ...
    '  f_v,pos   : position-only vertical fraction,\n', ...
    '  f_v,vel   : velocity-only vertical fraction.\n', ...
    'VerticalModeDominant still uses f_v,state > 0.5.\n']);

fprintf('============================================================\n');


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================



function runs = boolean_row_runs(mask)
%BOOLEAN_INTERVALS
% Return [startRow,endRow] for every contiguous true interval.

    mask = ...
        logical(mask(:));


    padded = ...
        [false;mask;false];


    d = ...
        diff(padded);


    iStart = ...
        find(d == 1);


    iEnd = ...
        find(d == -1) - 1;


    runs = ...
        [iStart,iEnd];

end


function T = row_run_table(runs,summaryTable,label)
%INTERVAL_TABLE

    if isempty(runs)

        T = ...
            table( ...
                zeros(0,1), ...
                zeros(0,1), ...
                zeros(0,1), ...
                zeros(0,1), ...
                zeros(0,1), ...
                zeros(0,1), ...
                zeros(0,1), ...
                strings(0,1), ...
                'VariableNames',{ ...
                    'StartRow', ...
                    'EndRow', ...
                    'StartJPL_Id', ...
                    'EndJPL_Id', ...
                    'StartZ0', ...
                    'EndZ0', ...
                    'Length', ...
                    'IntervalType'});

        return

    end


    n = ...
        size(runs,1);


    startRow = ...
        runs(:,1);


    endRow = ...
        runs(:,2);


    startId = ...
        summaryTable.JPL_Id(startRow);


    endId = ...
        summaryTable.JPL_Id(endRow);


    startZ0 = ...
        summaryTable.z0(startRow);


    endZ0 = ...
        summaryTable.z0(endRow);


    len = ...
        endRow - ...
        startRow + ...
        1;


    intervalType = ...
        repmat(string(label),n,1);


    T = ...
        table( ...
            startRow, ...
            endRow, ...
            startId, ...
            endId, ...
            startZ0, ...
            endZ0, ...
            len, ...
            intervalType, ...
            'VariableNames',{ ...
                'StartRow', ...
                'EndRow', ...
                'StartJPL_Id', ...
                'EndJPL_Id', ...
                'StartZ0', ...
                'EndZ0', ...
                'Length', ...
                'IntervalType'});

end


function add_highlight_rows_to_figure(fig,rows)
%ADD_HIGHLIGHT_ROWS_TO_FIGURE

    axesList = ...
        findall(fig,'Type','axes');


    for ia = 1:numel(axesList)

        ax = ...
            axesList(ia);


        if strcmp(get(ax,'Tag'),'Colorbar')
            continue
        end


        hold(ax,'on');


        for r = rows(:).'

            xline( ...
                ax, ...
                r, ...
                ':', ...
                sprintf('%d',r), ...
                'LabelVerticalAlignment','bottom', ...
                'HandleVisibility','off');

        end

    end

end
