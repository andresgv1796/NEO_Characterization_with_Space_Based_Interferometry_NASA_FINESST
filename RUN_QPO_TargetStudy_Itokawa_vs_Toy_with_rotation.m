% RUN_QPO_TargetStudy_Itokawa_vs_Toy_with_rotation.m
% Plain-MATLAB conversion of the user's Live Script, including the additive
% rotating-Itokawa diagnostic section. Existing analysis code is preserved.

% RUN_QPO_TargetStudy_Itokawa_vs_Toy_FULL
% Complete standalone FINESST experiment:
% parent periodic orbit
%     -> GMOS QPO
%     -> three equally phased collectors
%     -> propagated physical formation
%     -> all-sky target-suitability study
%     -> automatic selection of the highest-scoring valid target
%     -> accessibility + UVW for that selected target
%     -> Gaskell Itokawa projected into the selected target sky frame
%     -> physical QPO interferometric observation
%     -> ideal toy positive-control observation of the SAME angular truth
%     -> dirty imaging + Högbom CLEAN
%     -> native-resolution and common-resolution objective comparison

% No RA/Dec target is prescribed before the suitability study.
% The target study determines the celestial direction used by every subsequent part of the physical QPO/Itokawa pipeline.
% Uses the existing FINESST helper functions already used by the main pipeline.
clear; clc; close all;

%% Plot defaults
set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
set(groot,'DefaultAxesFontSize',18);
set(groot,'DefaultTextFontSize',20);
set(groot,'DefaultLineLineWidth',1.4);
set(groot,'DefaultFigureColor','w');
set(groot,'DefaultAxesColor','w');
set(groot,'DefaultAxesXColor','k');
set(groot,'DefaultAxesYColor','k');
set(groot,'DefaultTextColor','k');

% ============================================================
%% USER CONFIGURATION
% ============================================================

systemName = 'EarthMoon';
sys = cr3bp_system_parameters(systemName);
mu = sys.mu;

% Parent periodic orbit
poFile = 'L1_northern_halo_orbits_JPL_IC_Earth_Moon.csv';
orbitIndex = 2;

% GMOS QPO
qpoGen = struct();
qpoGen.N = 51;
qpoGen.mode = 'vertical';
qpoGen.initialAmplitude = 1e-5;
qpoGen.corrector = 'fixT';
qpoGen.CtargetMode = 'parent';
qpoGen.maxIter = 10;
qpoGen.tolInf = 1e-10;
qpoGen.verbose = true;
qpoGen.plotInitialGuess = true;
qpoGen.saveGeneratedQPO = false;

% Formation
nCollectors = 3;
phaseStartIndex = 1;
nPrimaryPeriods = 3;
nSavePerPrimaryPeriod = 500;
normalReferenceSynodic = [1 0 0];
timePlotUnit = 'auto';

% Interferometer
lambda_m = 10e-6;

% Accessibility
maxOffNormalDeg = 30.0;
usePlaneNormalSignAmbiguity = true;

% Target suitability study
%
% No target direction is supplied here.  The complete sky is sampled first,
% ranked with the existing FINESST target-suitability metric, and the best
% finite candidate is then used for the rest of the script.
nTargets = 10000;
useParallelScan = true;
nTopTargetsToPrint = 10;

skyNorthReference = [0 0 1];

% Itokawa
itokawaFile = fullfile('ShapeModels','Itokawa','ver64q.tab');
R_IB_Itokawa = eye(3);
nPixItokawa = 401;
marginItokawa = 1.20;
rangeItokawa_km = 1.0e8;
facetRadianceItokawa = 1.0;

% Ideal toy positive control
%
% These are the CURRENT study settings from the supplied script.
rhoMinToy_Glambda = 0.02;
rhoMaxToy_Glambda = 1.00;
nAnglesToy = 6;
nRadiiToy = 40;

% Imaging
addHermitian = true;
weightMode = 'uniform';
phaseAmpThreshold = 1e-5;

% Högbom CLEAN
cleanLoopGain = 0.10;
cleanNIterMax = 3000;
cleanStopFraction = 1e-3;
cleanBeamFitLevel = 0.50;
cleanBeamTruncateSigma = 5.0;

% Objective comparison metrics
correlationSupportLevel = 0.05;
contourLevel = 0.50;

% Integrator
odeOpts = odeset('RelTol',1e-12,'AbsTol',1e-12);

% ============================================================
%% 1. PARENT PERIODIC ORBIT + GMOS QPO
% ============================================================

po = choose_periodic_orbit(poFile,orbitIndex,sys);
xCombiner0 = po.x0(:);
TPrimary = po.T;

fprintf('\n============================================================\n');
fprintf('CR3BP / PARENT PERIODIC ORBIT\n');
fprintf('============================================================\n');
fprintf('  system       = %s\n',sys.name);
fprintf('  mu           = %.16e\n',sys.mu);
fprintf('  L*           = %.6f km\n',sys.Lstar_km);
fprintf('  t*           = %.6f days\n',sys.Tstar_days);
fprintf('  orbit file   = %s\n',po.file);
fprintf('  orbit index  = %d\n',orbitIndex);
fprintf('  TPrimary     = %.16e nondim\n',TPrimary);
fprintf('  TPrimary     = %.6f days\n',po.T_days);

qpo = generate_gmos_qpo_torus_from_po(po,sys,qpoGen,odeOpts);
XCurve = qpo.XCurve;
TQPO = qpo.T;
rhoQPO = qpo.rho;
N = qpo.N;

fprintf('\nGMOS QPO\n');
fprintf('  source       = %s\n',qpo.source);
fprintf('  N            = %d\n',N);
fprintf('  TQPO         = %.16e nondim\n',TQPO);
fprintf('  TQPO         = %.6f days\n',TQPO*sys.Tstar_days);
fprintf('  rhoQPO       = %.16e rad\n',rhoQPO);
fprintf('  rhoQPO/2pi   = %.16e cycles\n',rhoQPO/(2*pi));
fprintf('  period diff  = %.6e\n',abs(TQPO-TPrimary)/TPrimary);
fprintf('============================================================\n');

% ============================================================
%% 2. THREE EQUALLY PHASED COLLECTORS
% ============================================================

jCollectors = qpo_equal_phase_indices(N,nCollectors,phaseStartIndex);
XCollectors0 = XCurve(:,jCollectors);

fprintf('\nCollector QPO curve indices\n');
fprintf('  phaseStartIndex = %d\n',phaseStartIndex);
fprintf('  jCollectors     = [%d %d %d]\n',jCollectors(1),jCollectors(2),jCollectors(3));
fprintf('  phase fractions = [%.6f %.6f %.6f]\n', ...
    (jCollectors(1)-1)/N,(jCollectors(2)-1)/N,(jCollectors(3)-1)/N);

%% QPO invariant-curve / initial formation-section figure
baselineColors = lines(3);

figQPOSection = figure('Color','w','Position',[100 100 1350 600]);
tlQPOSection = tiledlayout(figQPOSection,1,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(tlQPOSection);
plot3(ax,XCurve(1,:)*sys.Lstar_km,XCurve(2,:)*sys.Lstar_km,XCurve(3,:)*sys.Lstar_km,'k-');
hold(ax,'on');
for k = 1:3
    scatter3(ax,XCollectors0(1,k)*sys.Lstar_km,XCollectors0(2,k)*sys.Lstar_km, ...
        XCollectors0(3,k)*sys.Lstar_km,75,baselineColors(k,:),'filled');
end
scatter3(ax,xCombiner0(1)*sys.Lstar_km,xCombiner0(2)*sys.Lstar_km, ...
    xCombiner0(3)*sys.Lstar_km,85,'k','filled');
axis(ax,'equal'); grid(ax,'on'); view(ax,3);
xlabel(ax,'$x$ [km]'); ylabel(ax,'$y$ [km]'); zlabel(ax,'$z$ [km]');
title(ax,'Corrected QPO invariant curve');

ax = nexttile(tlQPOSection);
XCurveRel_km = (XCurve(1:3,:) - xCombiner0(1:3))*sys.Lstar_km;
XCollectorsRel_km = (XCollectors0(1:3,:) - xCombiner0(1:3))*sys.Lstar_km;
plot3(ax,XCurveRel_km(1,:),XCurveRel_km(2,:),XCurveRel_km(3,:),'k-');
hold(ax,'on');
for k = 1:3
    scatter3(ax,XCollectorsRel_km(1,k),XCollectorsRel_km(2,k), ...
        XCollectorsRel_km(3,k),75,baselineColors(k,:),'filled');
end
scatter3(ax,0,0,0,85,'k','filled');
axis(ax,'equal'); grid(ax,'on'); view(ax,3);
xlabel(ax,'$\Delta x$ [km]'); ylabel(ax,'$\Delta y$ [km]'); zlabel(ax,'$\Delta z$ [km]');
title(ax,'Formation section relative to combiner');

title(tlQPOSection,sprintf('GMOS QPO section: orbit %d, $\\rho/2\\pi=%.4g$', ...
    orbitIndex,rhoQPO/(2*pi)));

% ============================================================
%% 3. PROPAGATE PHYSICAL QPO FORMATION
% ============================================================

geomQPO = compute_three_collector_qpo_geometry( ...
    mu,xCombiner0,TPrimary,XCollectors0, ...
    nPrimaryPeriods,nSavePerPrimaryPeriod,odeOpts,normalReferenceSynodic);

geomQPO.config.systemName = sys.name;
geomQPO.config.mu = mu;
geomQPO.config.N = N;
geomQPO.config.poFile = poFile;
geomQPO.config.orbitIndex = orbitIndex;
geomQPO.config.TPrimary = TPrimary;
geomQPO.config.TQPO = TQPO;
geomQPO.config.rhoQPO = rhoQPO;
geomQPO.config.qpoSource = 'generate';
geomQPO.config.qpoGen = qpoGen;
geomQPO.config.jCollectors = jCollectors;
geomQPO.config.phaseStartIndex = phaseStartIndex;
geomQPO.config.lambda_m = lambda_m;
geomQPO.config.normalReferenceSynodic = normalReferenceSynodic;

geomQPO = add_dimensional_geometry(geomQPO,sys,timePlotUnit);

% ============================================================
%% 4. TARGET SUITABILITY SCAN FOR THIS FORMATION
% ============================================================

fprintf('\n============================================================\n');
fprintf('TARGET SUITABILITY SCAN\n');
fprintf('============================================================\n');
fprintf('  sampled directions = %d\n',nTargets);

scanTargets = sample_sphere_targets(nTargets,'fibonacci');

tic
scan = evaluate_target_suitability_scan( ...
    geomQPO,sys,scanTargets,lambda_m, ...
    maxOffNormalDeg,usePlaneNormalSignAmbiguity,useParallelScan);
targetScanTime = toc;

% ------------------------------------------------------------------
% Rank only finite target-study results.
% ------------------------------------------------------------------

validScanMask = ...
    isfinite(scan.score) & ...
    isfinite(scan.raDeg) & ...
    isfinite(scan.decDeg) & ...
    isfinite(scan.accessibleFraction) & ...
    isfinite(scan.uvAnisotropy) & ...
    isfinite(scan.uvMax_Glambda);

validScanIdx = find(validScanMask);

if isempty(validScanIdx)

    error( ...
        'Target suitability study returned no finite candidate directions.');

end


[~,localSort] = ...
    sort( ...
        scan.score(validScanIdx), ...
        'descend');

idxSortScan = ...
    validScanIdx(localSort);

scan.idxSort = ...
    idxSortScan;


% ------------------------------------------------------------------
% Automatically select the highest-scoring candidate.
% ------------------------------------------------------------------

bestIdx = ...
    idxSortScan(1);

raDeg = ...
    scan.raDeg(bestIdx);

decDeg = ...
    scan.decDeg(bestIdx);


fprintf('\nBEST TARGET SELECTED BY THE ALL-SKY STUDY\n');

fprintf('  RA       = %.6f deg\n', ...
    raDeg);

fprintf('  Dec      = %.6f deg\n', ...
    decDeg);

fprintf('  Score    = %.8f\n', ...
    scan.score(bestIdx));

fprintf('  f_acc    = %.8f\n', ...
    scan.accessibleFraction(bestIdx));

fprintf('  eta_uv   = %.8f\n', ...
    scan.uvAnisotropy(bestIdx));

fprintf('  rho_max  = %.8f Glambda\n', ...
    scan.uvMax_Glambda(bestIdx));

fprintf('  runtime  = %.6f s\n', ...
    targetScanTime);

fprintf('============================================================\n');


% ------------------------------------------------------------------
% Print the leading targets so the selected direction can be inspected
% relative to nearby high-scoring alternatives.
% ------------------------------------------------------------------

nTop = ...
    min( ...
        nTopTargetsToPrint, ...
        numel(idxSortScan));

idxTop = ...
    idxSortScan(1:nTop);


targetRank = ...
    (1:nTop).';


topRA_deg = ...
    scan.raDeg(idxTop);
topRA_deg = ...
    topRA_deg(:);


topDec_deg = ...
    scan.decDeg(idxTop);
topDec_deg = ...
    topDec_deg(:);


topScore = ...
    scan.score(idxTop);
topScore = ...
    topScore(:);


topAccessibleFraction = ...
    scan.accessibleFraction(idxTop);
topAccessibleFraction = ...
    topAccessibleFraction(:);


topUVAnisotropy = ...
    scan.uvAnisotropy(idxTop);
topUVAnisotropy = ...
    topUVAnisotropy(:);


topRhoMax_Glambda = ...
    scan.uvMax_Glambda(idxTop);
topRhoMax_Glambda = ...
    topRhoMax_Glambda(:);


topTargetsTable = ...
    table( ...
        targetRank, ...
        topRA_deg, ...
        topDec_deg, ...
        topScore, ...
        topAccessibleFraction, ...
        topUVAnisotropy, ...
        topRhoMax_Glambda, ...
        'VariableNames',{ ...
            'Rank', ...
            'RA_deg', ...
            'Dec_deg', ...
            'Score', ...
            'AccessibleFraction', ...
            'UVAnisotropy', ...
            'RhoMax_Glambda'});


fprintf('\nTOP TARGET-STUDY CANDIDATES\n');

disp(topTargetsTable);


% Existing full-sky diagnostic plots.
figsTargetScan = ...
    plot_target_suitability_scan( ...
        scan);


fprintf('\nThe remainder of the pipeline uses the scan optimum only.\n');

fprintf('  selected RA  = %.6f deg\n', ...
    raDeg);

fprintf('  selected Dec = %.6f deg\n', ...
    decDeg);

% ============================================================
%% 5. STUDY-SELECTED FIXED TARGET: ACCESSIBILITY + UVW
% ============================================================

% From this point onward RA/Dec are no longer free inputs. They are the output of the target-suitability study above.
assert( ...
    bestIdx >= 1 && bestIdx <= numel(scan.score), ...
    'Invalid best-target index returned by the suitability study.');

assert( ...
    isfinite(raDeg) && isfinite(decDeg), ...
    'Selected target direction is not finite.');

target = make_target_definition( ...
    'inertial',true,raDeg,decDeg,[1 0 0],[1 0 0],skyNorthReference);

accessQPO = compute_triangle_accessibility( ...
    geomQPO,target,maxOffNormalDeg,usePlaneNormalSignAmbiguity);

uvQPO = compute_fixed_target_uv(geomQPO,target,sys,lambda_m);

print_geometry_accessibility_uv_summary(geomQPO,accessQPO,uvQPO);

% Full production formation / accessibility / UV plots
figsQPOGeometry = plot_accessibility_and_fixed_target_uv(geomQPO,accessQPO,uvQPO);

% ============================================================
%% 6. LOAD + ORIENT GASKELL ITOKAWA
% ============================================================

if ~isfile(itokawaFile)
    error('Could not find Gaskell Itokawa file:\n%s',itokawaFile);
end

shapeItokawa = load_gaskell_itokawa_shape(itokawaFile);

% Body -> inertial -> target-centered sky coordinates
verticesB = shapeItokawa.vertices;
verticesI = (R_IB_Itokawa*verticesB.').';

eU = target.eUInertial0(:);
eV = target.eVInertial0(:);
eW = target.eWInertial0(:);

verticesSky = [verticesI*eU,verticesI*eV,verticesI*eW];

%% 3-D body-frame + observing-orientation figure
figItokawa3D = figure('Color','w','Position',[100 100 1400 650]);
tl3D = tiledlayout(figItokawa3D,1,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl3D);
patch(ax,'Vertices',verticesB,'Faces',shapeItokawa.faces, ...
    'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
axis(ax,'equal'); grid(ax,'on'); view(ax,3);
xlabel(ax,'$x_B$ [km]'); ylabel(ax,'$y_B$ [km]'); zlabel(ax,'$z_B$ [km]');
title(ax,'Native Gaskell body frame');
camlight(ax,'headlight'); lighting(ax,'gouraud');

ax = nexttile(tl3D);

patch(ax,'Vertices',verticesSky,'Faces',shapeItokawa.faces, ...
    'FaceColor',[0.7 0.7 0.7], ...
    'EdgeColor','none');

axis(ax,'equal');
grid(ax,'on');
hold(ax,'on');

xlabel(ax,'$x_{\rm sky}\;(e_U)$ [km]');
ylabel(ax,'$y_{\rm sky}\;(e_V)$ [km]');
zlabel(ax,'$z_{\rm sky}\;(e_W)$ [km]');

% ------------------------------------------------------------------
% OBSERVING-GEOMETRY VISUALIZATION
%
% Sky-frame convention:
%
%   +e_U = image-plane horizontal direction
%   +e_V = image-plane vertical direction
%   +e_W = observer -> target / phase-center direction
%
% Therefore the physical observer lies on the -e_W side.
%
% This panel is deliberately shown from an OBLIQUE camera angle.
% If the MATLAB camera were placed at the physical observer, the +e_W
% line-of-sight arrow would lie exactly along the camera axis and would
% collapse visually to a point.  The actual observer-facing view is shown
% separately in the 2-D projection figure below.
% ------------------------------------------------------------------

bodyScale = max(vecnorm(verticesSky,2,2));

observerDistance = 2.2*bodyScale;
basisScale = 0.55*bodyScale;

observerPos = [0 0 -observerDistance];
targetPos = [0 0 0];

%% Physical observer-to-target line
plot3( ...
    ax, ...
    [observerPos(1) targetPos(1)], ...
    [observerPos(2) targetPos(2)], ...
    [observerPos(3) targetPos(3)], ...
    'k--', ...
    'LineWidth',1.4);

%% Observer marker
scatter3( ...
    ax, ...
    observerPos(1), ...
    observerPos(2), ...
    observerPos(3), ...
    75, ...
    'k', ...
    'filled');


text( ...
    ax, ...
    observerPos(1), ...
    observerPos(2), ...
    1.08*observerPos(3), ...
    'Observer', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle');

%% Line-of-sight arrow
% Arrow starts at the observer and points toward the target. In this sky frame that direction is exactly +e_W.
losLength = 0.75*observerDistance;

quiver3( ...
    ax, ...
    observerPos(1), ...
    observerPos(2), ...
    observerPos(3), ...
    0, ...
    0, ...
    losLength, ...
    0, ...
    'Color','k', ...
    'LineWidth',2.4, ...
    'MaxHeadSize',0.40);


text( ...
    ax, ...
    0.03*bodyScale, ...
    0.03*bodyScale, ...
    observerPos(3)+0.52*losLength, ...
    'view direction $=+e_W$', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom');

%% Target-centered sky-basis triad
quiver3( ...
    ax, ...
    0,0,0, ...
    basisScale,0,0, ...
    0, ...
    'Color',[0.0000 0.4470 0.7410], ...
    'LineWidth',2.2, ...
    'MaxHeadSize',0.45);

quiver3( ...
    ax, ...
    0,0,0, ...
    0,basisScale,0, ...
    0, ...
    'Color',[0.8500 0.3250 0.0980], ...
    'LineWidth',2.2, ...
    'MaxHeadSize',0.45);

quiver3( ...
    ax, ...
    0,0,0, ...
    0,0,basisScale, ...
    0, ...
    'Color',[0.9290 0.6940 0.1250], ...
    'LineWidth',2.2, ...
    'MaxHeadSize',0.45);


text(ax,1.10*basisScale,0,0,'$+e_U$');
text(ax,0,1.10*basisScale,0,'$+e_V$');
text(ax,0,0,1.10*basisScale,'$+e_W$');

%% Oblique camera used ONLY to visualize the geometry
view(ax,[-38 24]);

camproj(ax,'perspective');

%% Keep the observer and all basis arrows inside the displayed limits
xyzAll = [ ...
    verticesSky; ...
    observerPos; ...
    [ basisScale 0 0]; ...
    [0 basisScale 0]; ...
    [0 0 basisScale]];


spanAll = max(max(xyzAll,[],1)-min(xyzAll,[],1));

plotMargin = 0.10*spanAll;


xlim(ax,[ ...
    min(xyzAll(:,1))-plotMargin, ...
    max(xyzAll(:,1))+plotMargin]);

ylim(ax,[ ...
    min(xyzAll(:,2))-plotMargin, ...
    max(xyzAll(:,2))+plotMargin]);

zlim(ax,[ ...
    min(xyzAll(:,3))-plotMargin, ...
    max(xyzAll(:,3))+plotMargin]);


camlight(ax,'headlight');
lighting(ax,'gouraud');


title( ...
    ax, ...
    sprintf( ...
        ['Sky-frame observing geometry (oblique view), ', ...
         'RA $=%.3f^\circ$, Dec $=%+.3f^\circ$'], ...
        raDeg, ...
        decDeg));


title( ...
    tl3D, ...
    ['Gaskell Itokawa: native body frame and ', ...
     'target-centered observing geometry']);

% ============================================================
%% 7. PROJECT + RASTERIZE THE SAME STATIC ITOKAWA TRUTH ONCE
% ============================================================

projItokawa = project_shape_orthographic(shapeItokawa,R_IB_Itokawa,target);
visItokawa = find_visible_facets(projItokawa);

% Projected visible facets.  Plot the exact (x_sky,y_sky) coordinates used
% by the rasterization rather than relying on a separate 3-D view command.
verticesProj2D = [ ...
    projItokawa.x(:), ...
    projItokawa.y(:), ...
    zeros(numel(projItokawa.x),1)];

figure('Color','w');
patch('Vertices',verticesProj2D, ...
    'Faces',projItokawa.faces(visItokawa.isVisible,:), ...
    'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
view(2); axis equal; grid on;
xlabel('$x_{\\rm sky}$ [km]'); ylabel('$y_{\\rm sky}$ [km]');
title('Observer-facing projection: view along $+e_W$');

projectedExtent_km = max(abs([projItokawa.x(:);projItokawa.y(:)]));
halfWidthItokawa_km = marginItokawa*projectedExtent_km;

imgItokawa = rasterize_visible_facets( ...
    projItokawa,visItokawa,nPixItokawa,halfWidthItokawa_km,facetRadianceItokawa);

simulatedItokawa = make_simulated_target_from_raster( ...
    imgItokawa,rangeItokawa_km,shapeItokawa.type);

% Angular truth image.
%
% simulatedItokawa.I is a DISCRETE FRACTIONAL-FLUX image normalized so
% sum(I(:)) = 1.  Its pixel values are therefore generally far below unity.
% Normalize ONLY the displayed copy by its peak; do not modify the data used
% for visibility calculations.
rawTruthPeak = max(simulatedItokawa.I(:));

if rawTruthPeak <= 0
    error('Rasterized Itokawa truth has no positive pixels.');
end

rawTruthShow = simulatedItokawa.I/rawTruthPeak;

figure('Color','w');
imagesc(simulatedItokawa.l_mas,simulatedItokawa.m_mas,rawTruthShow);
axis image; set(gca,'YDir','normal');
colormap(gca,hot(256)); clim([0 1]);
xlabel('$l$ [mas]'); ylabel('$m$ [mas]');
title('Itokawa angular truth image');
cb = colorbar;
cb.Label.String = 'Relative intensity';

fprintf('  normalized image flux       = %.12f\\n',sum(simulatedItokawa.I(:)));
fprintf('  raw truth peak pixel        = %.12e\\n',rawTruthPeak);

fprintf('\n============================================================\n');
fprintf('ITOKAWA TARGET\n');
fprintf('============================================================\n');
fprintf('  selected RA               = %.6f deg\n',raDeg);
fprintf('  selected Dec              = %.6f deg\n',decDeg);
fprintf('  selection basis           = all-sky suitability maximum\n');
fprintf('  range                     = %.6e km\n',rangeItokawa_km);
fprintf('  mesh vertices             = %d\n',size(shapeItokawa.vertices,1));
fprintf('  mesh facets               = %d\n',size(shapeItokawa.faces,1));
fprintf('  observer-facing facets    = %d\n',nnz(visItokawa.isVisible));
fprintf('  image half-FOV            = %.6f mas\n',simulatedItokawa.fov_mas);
fprintf('============================================================\n');

%% Common baseline bookkeeping
baselineFields = {'uvw12','uvw13','uvw23'};
baselineNames = {'B12','B13','B23'};

% ============================================================

%% ============================================================
%  7A. STARTING-PHASE SENSITIVITY ON THE SAME QPO
%  ============================================================
%
% Hold fixed:
%   - the same parent periodic orbit / QPO family
%   - the same study-selected target (target)
%   - the same static Itokawa truth image (simulatedItokawa)
%
% Change only:
%   - the starting phase index used to pick the three equally spaced
%     collectors on the invariant curve
%
% This isolates how the initial phase on the SAME QPO affects the
% accessible physical UV coverage.

nPhaseCases = min(8,N);

phaseStartCases = unique( ...
    mod(round((0:nPhaseCases-1)*N/nPhaseCases),N) + 1, ...
    'stable');

% Force inclusion of the nominal case already used elsewhere.
phaseStartCases = unique([phaseStartIndex, phaseStartCases], 'stable');
nPhaseCases = numel(phaseStartCases);

geomPhaseCases   = cell(nPhaseCases,1);
accessPhaseCases = cell(nPhaseCases,1);
uvPhaseCases     = cell(nPhaseCases,1);
obsPhaseCases    = cell(nPhaseCases,1);
jCollectorsCases = zeros(nPhaseCases,nCollectors);

nAccessiblePhase = zeros(nPhaseCases,1);
nVisPhase        = zeros(nPhaseCases,1);
rhoMaxPhase      = zeros(nPhaseCases,1);
uvSpanPhase      = zeros(nPhaseCases,1);

fprintf('\n============================================================\n');
fprintf('QPO STARTING-PHASE SENSITIVITY STUDY\n');
fprintf('============================================================\n');
fprintf('  number of tested phase cases = %d\n', nPhaseCases);
fprintf('  tested phaseStartIndex values:\n    ');
fprintf('%d ', phaseStartCases);
fprintf('\n');

for kCase = 1:nPhaseCases

    thisPhaseStart = phaseStartCases(kCase);

    % Three equally spaced collectors on the same invariant curve,
    % shifted only by the starting phase.
    jCollectorsThis = qpo_equal_phase_indices(N,nCollectors,thisPhaseStart);
    XCollectors0This = XCurve(:,jCollectorsThis);

    jCollectorsCases(kCase,:) = jCollectorsThis;

    % Propagate the corresponding physical formation.
    geomThis = compute_three_collector_qpo_geometry( ...
        mu, xCombiner0, TPrimary, XCollectors0This, ...
        nPrimaryPeriods, nSavePerPrimaryPeriod, odeOpts, normalReferenceSynodic);

    geomThis.config = geomQPO.config;
    geomThis.config.jCollectors = jCollectorsThis;
    geomThis.config.phaseStartIndex = thisPhaseStart;

    geomThis = add_dimensional_geometry(geomThis,sys,timePlotUnit);

    % Same target as the nominal study-selected case.
    accessThis = compute_triangle_accessibility( ...
        geomThis, target, maxOffNormalDeg, usePlaneNormalSignAmbiguity);

    uvThis = compute_fixed_target_uv(geomThis,target,sys,lambda_m);

    % Build observations only so that we can visualize the ACTUAL physical
    % UV samples passed into the interferometric model.
    obsThis = make_simulated_interferometric_observations( ...
        simulatedItokawa, uvThis, accessThis, baselineFields, baselineNames);

    geomPhaseCases{kCase}   = geomThis;
    accessPhaseCases{kCase} = accessThis;
    uvPhaseCases{kCase}     = uvThis;
    obsPhaseCases{kCase}    = obsThis;

    nAccessiblePhase(kCase) = numel(obsThis.accessibleIndex);
    nVisPhase(kCase)        = numel(obsThis.V);
    rhoMaxPhase(kCase)      = max(obsThis.rho_Glambda);

    % Determine plot span from the actual physical UV samples.
    thisUVMax = 0;
    for b = 1:numel(obsThis.baselines)
        baseThis = obsThis.baselines(b);

        if isfield(baseThis,'u_Glambda') && isfield(baseThis,'v_Glambda')
            uPlot = baseThis.u_Glambda(:);
            vPlot = baseThis.v_Glambda(:);
        elseif isfield(baseThis,'u_Mlambda') && isfield(baseThis,'v_Mlambda')
            uPlot = baseThis.u_Mlambda(:)/1e3;
            vPlot = baseThis.v_Mlambda(:)/1e3;
        elseif isfield(baseThis,'u') && isfield(baseThis,'v')
            uPlot = baseThis.u(:)/1e9;
            vPlot = baseThis.v(:)/1e9;
        else
            error('Could not identify UV fields in obs.baselines(%d).', b);
        end

        thisUVMax = max(thisUVMax,max(abs([uPlot; vPlot])));
    end

    uvSpanPhase(kCase) = thisUVMax;

    fprintf('\n  case %d\n',kCase);
    fprintf('    phaseStartIndex           = %d\n', thisPhaseStart);
    fprintf('    jCollectors               = [%d %d %d]\n', ...
        jCollectorsThis(1),jCollectorsThis(2),jCollectorsThis(3));
    fprintf('    accessible epochs         = %d\n', nAccessiblePhase(kCase));
    fprintf('    physical visibilities     = %d\n', nVisPhase(kCase));
    fprintf('    max rho                   = %.6f Glambda\n', rhoMaxPhase(kCase));

end

phaseStudyTable = table( ...
    (1:nPhaseCases).', ...
    phaseStartCases(:), ...
    jCollectorsCases(:,1), ...
    jCollectorsCases(:,2), ...
    jCollectorsCases(:,3), ...
    nAccessiblePhase(:), ...
    nVisPhase(:), ...
    rhoMaxPhase(:), ...
    'VariableNames',{ ...
        'Case', ...
        'PhaseStartIndex', ...
        'j1', ...
        'j2', ...
        'j3', ...
        'AccessibleEpochs', ...
        'PhysicalVisibilities', ...
        'RhoMax_Glambda'});

fprintf('\nPHASE-SENSITIVITY SUMMARY TABLE\n');
disp(phaseStudyTable);

%% Initial collector locations on the invariant curve

figPhaseInit = figure('Color','w','Position',[100 100 1400 600]);
tlPhaseInit  = tiledlayout(figPhaseInit,1,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(tlPhaseInit);

plot3(ax, ...
    XCurve(1,:)*sys.Lstar_km, ...
    XCurve(2,:)*sys.Lstar_km, ...
    XCurve(3,:)*sys.Lstar_km, ...
    'k-','LineWidth',1.1);

hold(ax,'on');
grid(ax,'on');
axis(ax,'equal');
view(ax,3);

for kCase = 1:nPhaseCases
    X0This = XCurve(:,jCollectorsCases(kCase,:))*sys.Lstar_km;
    for j = 1:nCollectors
        scatter3(ax, ...
            X0This(1,j), X0This(2,j), X0This(3,j), ...
            48, baselineColors(j,:), 'filled', ...
            'MarkerEdgeColor','k');
    end
end

scatter3(ax, ...
    xCombiner0(1)*sys.Lstar_km, ...
    xCombiner0(2)*sys.Lstar_km, ...
    xCombiner0(3)*sys.Lstar_km, ...
    90,'k','filled');

xlabel(ax,'$x$ [km]');
ylabel(ax,'$y$ [km]');
zlabel(ax,'$z$ [km]');
title(ax,'Collector starting locations on the QPO');

ax = nexttile(tlPhaseInit);

XCurveRel_km = (XCurve(1:3,:) - xCombiner0(1:3))*sys.Lstar_km;
plot3(ax, ...
    XCurveRel_km(1,:), ...
    XCurveRel_km(2,:), ...
    XCurveRel_km(3,:), ...
    'k-','LineWidth',1.1);

hold(ax,'on');
grid(ax,'on');
axis(ax,'equal');
view(ax,3);

for kCase = 1:nPhaseCases
    X0ThisRel = (XCurve(1:3,jCollectorsCases(kCase,:)) - xCombiner0(1:3))*sys.Lstar_km;
    for j = 1:nCollectors
        scatter3(ax, ...
            X0ThisRel(1,j), X0ThisRel(2,j), X0ThisRel(3,j), ...
            48, baselineColors(j,:), 'filled', ...
            'MarkerEdgeColor','k');
    end
end

scatter3(ax,0,0,0,90,'k','filled');

xlabel(ax,'$\Delta x$ [km]');
ylabel(ax,'$\Delta y$ [km]');
zlabel(ax,'$\Delta z$ [km]');
title(ax,'Same starting locations relative to the combiner');

title(tlPhaseInit,'QPO starting-phase cases on the same invariant curve');

%% Side-by-side UV coverage comparison
%
% Use the actual physical UV samples appearing in the observation model.

figPhaseUV = figure('Color','w','Position',[100 100 1700 900]);

nColsPhase = ceil(sqrt(nPhaseCases));
nRowsPhase = ceil(nPhaseCases/nColsPhase);

tlPhaseUV = tiledlayout(figPhaseUV,nRowsPhase,nColsPhase, ...
    'TileSpacing','compact','Padding','compact');

uvPlotMax = 1.05*max(uvSpanPhase);

for kCase = 1:nPhaseCases

    ax = nexttile(tlPhaseUV);
    hold(ax,'on');
    grid(ax,'on');

    obsThis = obsPhaseCases{kCase};

    legendHandles = gobjects(numel(obsThis.baselines),1);

    for b = 1:numel(obsThis.baselines)

        baseThis = obsThis.baselines(b);

        if isfield(baseThis,'u_Glambda') && isfield(baseThis,'v_Glambda')
            uPlot = baseThis.u_Glambda(:);
            vPlot = baseThis.v_Glambda(:);
        elseif isfield(baseThis,'u_Mlambda') && isfield(baseThis,'v_Mlambda')
            uPlot = baseThis.u_Mlambda(:)/1e3;
            vPlot = baseThis.v_Mlambda(:)/1e3;
        elseif isfield(baseThis,'u') && isfield(baseThis,'v')
            uPlot = baseThis.u(:)/1e9;
            vPlot = baseThis.v(:)/1e9;
        else
            error('Could not identify UV fields in obs.baselines(%d).', b);
        end

        legendHandles(b) = scatter(ax, ...
            uPlot, vPlot, ...
            10, baselineColors(b,:), 'filled');
    end

    axis(ax,'equal');
    xlim(ax,[-uvPlotMax uvPlotMax]);
    ylim(ax,[-uvPlotMax uvPlotMax]);

    xlabel(ax,'$u$ [G$\lambda$]');
    ylabel(ax,'$v$ [G$\lambda$]');

    title(ax,sprintf(['start %d\n', ...
        '$N_{\\rm acc}=%d$, $N_V=%d$, $\\rho_{\\max}=%.3f$ G$\\lambda$'], ...
        phaseStartCases(kCase), ...
        nAccessiblePhase(kCase), ...
        nVisPhase(kCase), ...
        rhoMaxPhase(kCase)));

    if kCase == 1
        legend(ax,legendHandles,baselineNames,'Location','best');
    end

end

title(tlPhaseUV,'Effect of QPO starting phase on the physical UV coverage');

%% Optional: full target + UV figure for each phase case
%
% Uncomment if you want the same style as the rest of the pipeline.
%
% for kCase = 1:nPhaseCases
%     figTmp = plot_selected_simulated_target_and_uv_samples( ...
%         simulatedItokawa, ...
%         uvPhaseCases{kCase}, ...
%         accessPhaseCases{kCase}, ...
%         baselineFields, ...
%         baselineNames);
%
%     sgtitle(sprintf('Itokawa target and UV samples: phaseStartIndex = %d', ...
%         phaseStartCases(kCase)));
% end
%% 8. CASE A: PHYSICAL GMOS QPO OBSERVATION
% ============================================================

obsQPO = make_simulated_interferometric_observations( ...
    simulatedItokawa,uvQPO,accessQPO,baselineFields,baselineNames);

print_image_uv_sampling_summary(obsQPO.gridCheck);

fprintf('\nCURRENT QPO OBSERVATION\n');
fprintf('  accessible epochs           = %d\n',numel(obsQPO.accessibleIndex));
fprintf('  physical visibility samples = %d\n',numel(obsQPO.V));
fprintf('  max rho                     = %.6f Glambda\n',max(obsQPO.rho_Glambda));

% Full production figures
figQPOTargetUV = plot_selected_simulated_target_and_uv_samples( ...
    simulatedItokawa,uvQPO,accessQPO,baselineFields,baselineNames);

figsQPOObs = plot_simulated_interferometric_observations(obsQPO,phaseAmpThreshold);

dirtyQPO = make_dirty_image_from_observations( ...
    obsQPO,simulatedItokawa.l_rad,simulatedItokawa.m_rad,addHermitian,weightMode);

figsDirtyQPO = plot_dirty_image_results(simulatedItokawa,dirtyQPO);

cleanParamsQPO = struct();
cleanParamsQPO.loopGain = cleanLoopGain;
cleanParamsQPO.nIterMax = cleanNIterMax;

dirtyBeamPeakQPO = max(real(dirtyQPO.beam(:)));
if dirtyBeamPeakQPO <= 0
    error('QPO dirty beam must have a positive peak.');
end

dirtyImageNormalizedQPO = real(dirtyQPO.image)/dirtyBeamPeakQPO;
cleanParamsQPO.stopThreshold = cleanStopFraction*max(abs(dirtyImageNormalizedQPO(:)));
cleanParamsQPO.cleanMask = true(size(dirtyQPO.image));
cleanParamsQPO.peakMode = 'absolute';
cleanParamsQPO.cleanBeamFitLevel = cleanBeamFitLevel;
cleanParamsQPO.cleanBeamTruncateSigma = cleanBeamTruncateSigma;

cleanQPO = make_clean_image_from_dirty(dirtyQPO,cleanParamsQPO);
figsCleanQPO = plot_clean_image_results(simulatedItokawa,dirtyQPO,cleanQPO);

% ============================================================
%% 9. CASE B: CONFIGURABLE IDEAL TOY POSITIVE CONTROL
% ============================================================

geomToy = make_toy_equilateral_uv_sweep_geometry( ...
    sys,lambda_m,rhoMinToy_Glambda,rhoMaxToy_Glambda,nAnglesToy,nRadiiToy);

% Canonical face-on toy target.  These u,v samples are applied to the SAME
% angular Itokawa image generated above.
targetToy = make_target_definition( ...
    'synodic',false,0,0,[1 0 0],[0 0 1],[0 1 0]);

accessToy = compute_triangle_accessibility(geomToy,targetToy,1e-6,true);
uvToy = compute_fixed_target_uv(geomToy,targetToy,sys,lambda_m);

assert(all(accessToy.accessible),'Every toy epoch should be accessible.');

maxAbsWToy = max(abs([uvToy.uvw12(:,3);uvToy.uvw13(:,3);uvToy.uvw23(:,3)]));
closureToy = uvToy.uvw12 + uvToy.uvw23 - uvToy.uvw13;

fprintf('\nIDEAL TOY FORMATION\n');
fprintf('  radial rings         = %d\n',nRadiiToy);
fprintf('  angular samples      = %d\n',nAnglesToy);
fprintf('  trajectory epochs    = %d\n',numel(geomToy.t));
fprintf('  accessible epochs    = %d\n',nnz(accessToy.accessible));
fprintf('  max |w|              = %.6e wavelengths\n',maxAbsWToy);
fprintf('  max closure error    = %.6e wavelengths\n',max(vecnorm(closureToy,2,2)));

obsToy = make_simulated_interferometric_observations( ...
    simulatedItokawa,uvToy,accessToy,baselineFields,baselineNames);

print_image_uv_sampling_summary(obsToy.gridCheck);

fprintf('\nTOY OBSERVATION\n');
fprintf('  physical visibility samples = %d\n',numel(obsToy.V));
fprintf('  max rho                     = %.6f Glambda\n',max(obsToy.rho_Glambda));

% Full production figures
figToyTargetUV = plot_selected_simulated_target_and_uv_samples( ...
    simulatedItokawa,uvToy,accessToy,baselineFields,baselineNames);

figsToyObs = plot_simulated_interferometric_observations(obsToy,phaseAmpThreshold);

dirtyToy = make_dirty_image_from_observations( ...
    obsToy,simulatedItokawa.l_rad,simulatedItokawa.m_rad,addHermitian,weightMode);

figsDirtyToy = plot_dirty_image_results(simulatedItokawa,dirtyToy);

cleanParamsToy = struct();
cleanParamsToy.loopGain = cleanLoopGain;
cleanParamsToy.nIterMax = cleanNIterMax;

dirtyBeamPeakToy = max(real(dirtyToy.beam(:)));
if dirtyBeamPeakToy <= 0
    error('Toy dirty beam must have a positive peak.');
end

dirtyImageNormalizedToy = real(dirtyToy.image)/dirtyBeamPeakToy;
cleanParamsToy.stopThreshold = cleanStopFraction*max(abs(dirtyImageNormalizedToy(:)));
cleanParamsToy.cleanMask = true(size(dirtyToy.image));
cleanParamsToy.peakMode = 'absolute';
cleanParamsToy.cleanBeamFitLevel = cleanBeamFitLevel;
cleanParamsToy.cleanBeamTruncateSigma = cleanBeamTruncateSigma;

cleanToy = make_clean_image_from_dirty(dirtyToy,cleanParamsToy);
figsCleanToy = plot_clean_image_results(simulatedItokawa,dirtyToy,cleanToy);

% ============================================================
%% 10. NATIVE-RESOLUTION OBJECTIVE METRICS
% ============================================================

% For quantitative morphology metrics use the restored CLEAN component model,
% not the final image with the residual added.  restoredComponents has exactly
% the fitted CLEAN beam as its PSF, so comparison with I_true * B_clean is
% mathematically well-defined.  The full production CLEAN figures still show
% the final restored image including the residual.
imageQPO = real(cleanQPO.restoredComponents);
imageToy = real(cleanToy.restoredComponents);

truthQPO = real(conv2(simulatedItokawa.I,cleanQPO.cleanBeamKernel,'same'));
truthToy = real(conv2(simulatedItokawa.I,cleanToy.cleanBeamKernel,'same'));

% Global amplitude fits for morphology comparison
alphaQPO = sum(truthQPO(:).*imageQPO(:))/sum(imageQPO(:).^2);
alphaToy = sum(truthToy(:).*imageToy(:))/sum(imageToy(:).^2);

imageQPOScaled = alphaQPO*imageQPO;
imageToyScaled = alphaToy*imageToy;

peakTruthQPO = max(truthQPO(:));
peakTruthToy = max(truthToy(:));

nrmseQPONative = sqrt(mean((imageQPOScaled(:)-truthQPO(:)).^2))/peakTruthQPO;
nrmseToyNative = sqrt(mean((imageToyScaled(:)-truthToy(:)).^2))/peakTruthToy;

% CLEAN beam resolution metrics
beamMajorQPO = cleanQPO.params.cleanBeamFwhmMajor_mas;
beamMinorQPO = cleanQPO.params.cleanBeamFwhmMinor_mas;
beamAxisRatioQPO = cleanQPO.params.cleanBeamAxisRatio;
beamPaQPO = cleanQPO.params.cleanBeamPaDeg;

beamMajorToy = cleanToy.params.cleanBeamFwhmMajor_mas;
beamMinorToy = cleanToy.params.cleanBeamFwhmMinor_mas;
beamAxisRatioToy = cleanToy.params.cleanBeamAxisRatio;
beamPaToy = cleanToy.params.cleanBeamPaDeg;

beamAreaQPO_mas2 = pi/(4*log(2))*beamMajorQPO*beamMinorQPO;
beamAreaToy_mas2 = pi/(4*log(2))*beamMajorToy*beamMinorToy;
beamAreaRatioQPOtoToy = beamAreaQPO_mas2/beamAreaToy_mas2;

% ============================================================
%% 11. COMMON-RESOLUTION COMPARISON WITH A SINGLE EXPLICIT BEAM
% ============================================================

% The previous version formed
%      B_common = B_QPO * B_toy,

% which is mathematically valid but unnecessarily blurs both images. It can erase precisely the morphology we are trying to compare.
% Instead define ONE explicit common circular Gaussian beam whose FWHM is only slightly larger than the largest native CLEAN-beam major axis. Since a circular covariance of this size dominates both native Gaussian beam covariance matrices, each reconstruction can be smoothed by an additional Gaussian so that both end with exactly the same common covariance.
% RAW TRUTH:
%      I_true = simulatedItokawa.I

% NATIVE-BEAM TRUTHS:
%      I_true,QPO = I_true * B_QPO
%      I_true,toy = I_true * B_toy

% COMMON TRUTH:
%      I_true,common = I_true * B_common

% Nothing new is inferred about Itokawa in these "truth" images. They are simply the same known rasterized input morphology shown at the angular resolution appropriate to the comparison.
fwhmToSigma = 1/(2*sqrt(2*log(2)));

SigmaQPO_mas2 = clean_beam_covariance_mas2( ...
    beamMajorQPO,beamMinorQPO,beamPaQPO);

SigmaToy_mas2 = clean_beam_covariance_mas2( ...
    beamMajorToy,beamMinorToy,beamPaToy);

% Robust common beam: 2 percent broader than the largest native major FWHM.
% The small margin prevents a nearly singular extra-blur kernel when one
% native major axis is exactly the limiting value.
commonBeamFwhm_mas = 1.02*max([beamMajorQPO,beamMajorToy]);
commonBeamSigma_mas = fwhmToSigma*commonBeamFwhm_mas;
SigmaCommon_mas2 = commonBeamSigma_mas^2*eye(2);

SigmaExtraQPO_mas2 = SigmaCommon_mas2 - SigmaQPO_mas2;
SigmaExtraToy_mas2 = SigmaCommon_mas2 - SigmaToy_mas2;

% Numerical PSD checks.
minEigExtraQPO = min(eig(0.5*(SigmaExtraQPO_mas2+SigmaExtraQPO_mas2.')));
minEigExtraToy = min(eig(0.5*(SigmaExtraToy_mas2+SigmaExtraToy_mas2.')));

if minEigExtraQPO < -1e-12 || minEigExtraToy < -1e-12
    error(['Chosen common beam does not dominate both native beams. ', ...
           'Increase commonBeamFwhm_mas.']);
end

% Angular pixel spacing of the reconstructed image.
dl_mas = mean(diff(simulatedItokawa.l_mas));
dm_mas = mean(diff(simulatedItokawa.m_mas));

commonBeamKernel = gaussian_kernel_from_covariance_mas( ...
    SigmaCommon_mas2,dl_mas,dm_mas,5.0);

extraQPOKernel = gaussian_kernel_from_covariance_mas( ...
    SigmaExtraQPO_mas2,dl_mas,dm_mas,5.0);

extraToyKernel = gaussian_kernel_from_covariance_mas( ...
    SigmaExtraToy_mas2,dl_mas,dm_mas,5.0);

% One common truth, generated directly from the RAW truth.
truthCommon = real(conv2( ...
    simulatedItokawa.I, ...
    commonBeamKernel, ...
    'same'));

% Bring each native CLEAN image to that same common beam.
imageQPOCommon = real(conv2( ...
    imageQPO, ...
    extraQPOKernel, ...
    'same'));

imageToyCommon = real(conv2( ...
    imageToy, ...
    extraToyKernel, ...
    'same'));

% Independent morphology amplitude fits.
alphaQPOCommon = ...
    sum(truthCommon(:).*imageQPOCommon(:))/sum(imageQPOCommon(:).^2);

alphaToyCommon = ...
    sum(truthCommon(:).*imageToyCommon(:))/sum(imageToyCommon(:).^2);

imageQPOCommonScaled = alphaQPOCommon*imageQPOCommon;
imageToyCommonScaled = alphaToyCommon*imageToyCommon;

peakTruthCommon = max(truthCommon(:));

if peakTruthCommon <= 0
    error('Common-resolution truth has no positive peak.');
end

errorQPOCommon = imageQPOCommonScaled-truthCommon;
errorToyCommon = imageToyCommonScaled-truthCommon;

% Global common-resolution NRMSE.
nrmseQPOCommon = sqrt(mean(errorQPOCommon(:).^2))/peakTruthCommon;
nrmseToyCommon = sqrt(mean(errorToyCommon(:).^2))/peakTruthCommon;

% Support-limited NRMSE.  This is often more discriminating than the global
% value because the large empty background no longer dilutes the metric.
metricSupportMask = ...
    truthCommon >= correlationSupportLevel*peakTruthCommon;

nrmseQPOCommonSupport = ...
    sqrt(mean(errorQPOCommon(metricSupportMask).^2))/peakTruthCommon;

nrmseToyCommonSupport = ...
    sqrt(mean(errorToyCommon(metricSupportMask).^2))/peakTruthCommon;

% NCC on the same meaningful target support.
truthVector = truthCommon(metricSupportMask);
qpoVector = imageQPOCommonScaled(metricSupportMask);
toyVector = imageToyCommonScaled(metricSupportMask);

truthVector0 = truthVector-mean(truthVector);
qpoVector0 = qpoVector-mean(qpoVector);
toyVector0 = toyVector-mean(toyVector);

nccQPOCommon = dot(truthVector0,qpoVector0)/(norm(truthVector0)*norm(qpoVector0));
nccToyCommon = dot(truthVector0,toyVector0)/(norm(truthVector0)*norm(toyVector0));

% Contour IoU at the requested common threshold.
truthMaskContour = truthCommon >= contourLevel*peakTruthCommon;
qpoMaskContour = imageQPOCommonScaled >= contourLevel*peakTruthCommon;
toyMaskContour = imageToyCommonScaled >= contourLevel*peakTruthCommon;

iouQPOCommon = ...
    nnz(truthMaskContour & qpoMaskContour)/nnz(truthMaskContour | qpoMaskContour);

iouToyCommon = ...
    nnz(truthMaskContour & toyMaskContour)/nnz(truthMaskContour | toyMaskContour);

fprintf('\nCOMMON RESOLUTION\n');
fprintf('  common circular FWHM        = %.6f mas\n',commonBeamFwhm_mas);
fprintf('  QPO support NRMSE           = %.8e\n',nrmseQPOCommonSupport);
fprintf('  toy support NRMSE           = %.8e\n',nrmseToyCommonSupport);

% ============================================================
%% 12. FINAL NUMERICAL SUMMARY + TABLE
% ============================================================

fprintf('\n============================================================\n');
fprintf('FINAL ITOKAWA RECONSTRUCTION COMPARISON\n');
fprintf('============================================================\n');

fprintf('  scan-selected RA            = %.6f deg\n', ...
    raDeg);

fprintf('  scan-selected Dec           = %.6f deg\n', ...
    decDeg);

fprintf('  target suitability score    = %.8f\n', ...
    scan.score(bestIdx));

fprintf('  target accessible fraction  = %.8f\n', ...
    scan.accessibleFraction(bestIdx));

fprintf('  target UV anisotropy        = %.8f\n', ...
    scan.uvAnisotropy(bestIdx));

fprintf('  target rho max              = %.8f Glambda\n', ...
    scan.uvMax_Glambda(bestIdx));

fprintf('\nCURRENT GMOS QPO\n');
fprintf('  visibility samples          = %d\n',numel(obsQPO.V));
fprintf('  max rho                     = %.6f Glambda\n',max(obsQPO.rho_Glambda));
fprintf('  clean beam                  = %.6f x %.6f mas\n',beamMajorQPO,beamMinorQPO);
fprintf('  beam axis ratio             = %.6f\n',beamAxisRatioQPO);
fprintf('  beam PA                     = %.6f deg\n',beamPaQPO);
fprintf('  beam area                   = %.8e mas^2\n',beamAreaQPO_mas2);
fprintf('  restored-model native NRMSE = %.8e\n',nrmseQPONative);
fprintf('  final residual RMS           = %.8e\n',sqrt(mean(abs(cleanQPO.residual(:)).^2)));
fprintf('  common-resolution NRMSE     = %.8e\n',nrmseQPOCommon);
fprintf('  common support NRMSE        = %.8e\n',nrmseQPOCommonSupport);
fprintf('  common-resolution NCC       = %.8f\n',nccQPOCommon);
fprintf('  common-resolution IoU %.0f%% = %.8f\n',100*contourLevel,iouQPOCommon);

fprintf('\nIDEAL TOY\n');
fprintf('  visibility samples          = %d\n',numel(obsToy.V));
fprintf('  max rho                     = %.6f Glambda\n',max(obsToy.rho_Glambda));
fprintf('  clean beam                  = %.6f x %.6f mas\n',beamMajorToy,beamMinorToy);
fprintf('  beam axis ratio             = %.6f\n',beamAxisRatioToy);
fprintf('  beam PA                     = %.6f deg\n',beamPaToy);
fprintf('  beam area                   = %.8e mas^2\n',beamAreaToy_mas2);
fprintf('  restored-model native NRMSE = %.8e\n',nrmseToyNative);
fprintf('  final residual RMS           = %.8e\n',sqrt(mean(abs(cleanToy.residual(:)).^2)));
fprintf('  common-resolution NRMSE     = %.8e\n',nrmseToyCommon);
fprintf('  common support NRMSE        = %.8e\n',nrmseToyCommonSupport);
fprintf('  common-resolution NCC       = %.8f\n',nccToyCommon);
fprintf('  common-resolution IoU %.0f%% = %.8f\n',100*contourLevel,iouToyCommon);

fprintf('\nDIRECT COMPARISON\n');
fprintf('  beam area QPO / toy         = %.6f\n',beamAreaRatioQPOtoToy);
fprintf('  native NRMSE QPO / toy      = %.6f\n',nrmseQPONative/nrmseToyNative);
fprintf('  common NRMSE QPO / toy      = %.6f\n',nrmseQPOCommon/nrmseToyCommon);
fprintf('  common support NRMSE ratio  = %.6f\n',nrmseQPOCommonSupport/nrmseToyCommonSupport);
fprintf('  common NCC QPO - toy        = %+.8f\n',nccQPOCommon-nccToyCommon);
fprintf('  common IoU QPO - toy        = %+.8f\n',iouQPOCommon-iouToyCommon);
fprintf('============================================================\n');

comparisonTable = table( ...
    ["Current QPO";"Ideal toy"], ...
    [numel(obsQPO.V);numel(obsToy.V)], ...
    [max(obsQPO.rho_Glambda);max(obsToy.rho_Glambda)], ...
    [beamMajorQPO;beamMajorToy], ...
    [beamMinorQPO;beamMinorToy], ...
    [beamAreaQPO_mas2;beamAreaToy_mas2], ...
    [nrmseQPONative;nrmseToyNative], ...
    [nrmseQPOCommon;nrmseToyCommon], ...
    [nrmseQPOCommonSupport;nrmseToyCommonSupport], ...
    [nccQPOCommon;nccToyCommon], ...
    [iouQPOCommon;iouToyCommon], ...
    'VariableNames',{ ...
        'Formation','PhysicalSamples','MaxRho_Glambda', ...
        'BeamMajor_mas','BeamMinor_mas','BeamArea_mas2', ...
        'NativeNRMSE','CommonNRMSE','CommonSupportNRMSE', ...
        'CommonNCC','CommonIoU'});

disp(comparisonTable);

% ============================================================
%% 13. CONSISTENT NATIVE-RESOLUTION SUMMARY FIGURE
% ============================================================

% Sequential positive map
positiveMap = hot(256);

% Diverging blue-white-red map
nColor = 256;
nHalf = floor(nColor/2);
blueToWhite = [linspace(0,1,nHalf).',linspace(0,1,nHalf).',ones(nHalf,1)];
whiteToRed = [ones(nColor-nHalf,1),linspace(1,0,nColor-nHalf).',linspace(1,0,nColor-nHalf).'];
divergingMap = [blueToWhite;whiteToRed];

baselineColors = lines(3);
baselineNamesComparison = {'B12','B13','B23'};
uvwQPOCell = {uvQPO.uvw12,uvQPO.uvw13,uvQPO.uvw23};
uvwToyCell = {uvToy.uvw12,uvToy.uvw13,uvToy.uvw23};
maskQPO = accessQPO.accessible(:);
maskToy = accessToy.accessible(:);

rhoPlotMax = 1.05*max([obsQPO.rho_Glambda(:);obsToy.rho_Glambda(:)]);

figNativeComparison = figure('Color','w','Position',[70 70 1700 850]);
tlNative = tiledlayout(figNativeComparison,2,4,'TileSpacing','compact','Padding','compact');

% QPO UV
ax1 = nexttile(tlNative); hold(ax1,'on'); grid(ax1,'on');
for b = 1:3
    uvwCurrent = uvwQPOCell{b};
    uCurrent = uvwCurrent(maskQPO,1)/1e9;
    vCurrent = uvwCurrent(maskQPO,2)/1e9;
    plot(ax1,uCurrent,vCurrent,'.','Color',baselineColors(b,:), ...
        'MarkerSize',8,'DisplayName',baselineNamesComparison{b});
    plot(ax1,-uCurrent,-vCurrent,'.','Color',baselineColors(b,:), ...
        'MarkerSize',8,'HandleVisibility','off');
end
axis(ax1,'equal'); xlim(ax1,[-rhoPlotMax rhoPlotMax]); ylim(ax1,[-rhoPlotMax rhoPlotMax]);
xlabel(ax1,'$u$ [G$\lambda$]'); ylabel(ax1,'$v$ [G$\lambda$]');
title(ax1,'Current QPO: UV'); legend(ax1,'Location','best');

% QPO beam
ax2 = nexttile(tlNative);
beamQPO = real(dirtyQPO.beam);
imagesc(ax2,dirtyQPO.l_mas,dirtyQPO.m_mas,beamQPO);
axis(ax2,'image'); set(ax2,'YDir','normal');
clim(ax2,[-max(abs(beamQPO(:))) max(abs(beamQPO(:)))]);
colormap(ax2,divergingMap);
xlabel(ax2,'$l$ [mas]'); ylabel(ax2,'$m$ [mas]');
title(ax2,sprintf('QPO beam: %.3f $\\times$ %.3f mas',beamMajorQPO,beamMinorQPO));
colorbar(ax2);

% QPO truth
ax3 = nexttile(tlNative);
imagesc(ax3,simulatedItokawa.l_mas,simulatedItokawa.m_mas,truthQPO/peakTruthQPO);
axis(ax3,'image'); set(ax3,'YDir','normal'); colormap(ax3,positiveMap); clim(ax3,[0 1]);
xlabel(ax3,'$l$ [mas]'); ylabel(ax3,'$m$ [mas]'); title(ax3,'Truth at QPO resolution'); colorbar(ax3);

% QPO CLEAN
ax4 = nexttile(tlNative);
imagesc(ax4,simulatedItokawa.l_mas,simulatedItokawa.m_mas,imageQPOScaled/peakTruthQPO);
axis(ax4,'image'); set(ax4,'YDir','normal'); colormap(ax4,positiveMap); clim(ax4,[0 1]);
xlabel(ax4,'$l$ [mas]'); ylabel(ax4,'$m$ [mas]');
title(ax4,{'QPO restored CLEAN model',sprintf('native NRMSE $=%.3g$',nrmseQPONative)}); colorbar(ax4);

% Toy UV
ax5 = nexttile(tlNative); hold(ax5,'on'); grid(ax5,'on');
for b = 1:3
    uvwCurrent = uvwToyCell{b};
    uCurrent = uvwCurrent(maskToy,1)/1e9;
    vCurrent = uvwCurrent(maskToy,2)/1e9;
    plot(ax5,uCurrent,vCurrent,'.','Color',baselineColors(b,:), ...
        'MarkerSize',8,'DisplayName',baselineNamesComparison{b});
    plot(ax5,-uCurrent,-vCurrent,'.','Color',baselineColors(b,:), ...
        'MarkerSize',8,'HandleVisibility','off');
end
axis(ax5,'equal'); xlim(ax5,[-rhoPlotMax rhoPlotMax]); ylim(ax5,[-rhoPlotMax rhoPlotMax]);
xlabel(ax5,'$u$ [G$\lambda$]'); ylabel(ax5,'$v$ [G$\lambda$]');
title(ax5,'Ideal toy: UV'); legend(ax5,'Location','best');

% Toy beam
ax6 = nexttile(tlNative);
beamToy = real(dirtyToy.beam);
imagesc(ax6,dirtyToy.l_mas,dirtyToy.m_mas,beamToy);
axis(ax6,'image'); set(ax6,'YDir','normal');
clim(ax6,[-max(abs(beamToy(:))) max(abs(beamToy(:)))]);
colormap(ax6,divergingMap);
xlabel(ax6,'$l$ [mas]'); ylabel(ax6,'$m$ [mas]');
title(ax6,sprintf('Toy beam: %.3f $\\times$ %.3f mas',beamMajorToy,beamMinorToy));
colorbar(ax6);

% Toy truth
ax7 = nexttile(tlNative);
imagesc(ax7,simulatedItokawa.l_mas,simulatedItokawa.m_mas,truthToy/peakTruthToy);
axis(ax7,'image'); set(ax7,'YDir','normal'); colormap(ax7,positiveMap); clim(ax7,[0 1]);
xlabel(ax7,'$l$ [mas]'); ylabel(ax7,'$m$ [mas]'); title(ax7,'Truth at toy resolution'); colorbar(ax7);

% Toy CLEAN
ax8 = nexttile(tlNative);
imagesc(ax8,simulatedItokawa.l_mas,simulatedItokawa.m_mas,imageToyScaled/peakTruthToy);
axis(ax8,'image'); set(ax8,'YDir','normal'); colormap(ax8,positiveMap); clim(ax8,[0 1]);
xlabel(ax8,'$l$ [mas]'); ylabel(ax8,'$m$ [mas]');
title(ax8,{'Toy restored CLEAN model',sprintf('native NRMSE $=%.3g$',nrmseToyNative)}); colorbar(ax8);

title(tlNative,sprintf( ...
    ['Scan-selected Gaskell Itokawa target: GMOS QPO versus ideal coverage, ', ...
     'RA $=%.3f^\\circ$, Dec $=%+.3f^\\circ$'],raDeg,decDeg));

% ============================================================
%% 14. COMMON-RESOLUTION OBJECTIVE COMPARISON FIGURE
% ============================================================

% The first panel explicitly shows the RAW rasterized truth used to generate all visibilities. The second panel shows that same truth convolved once with the explicitly chosen common circular beam. The QPO and toy panels can then be compared against exactly that same common-resolution reference.
figCommonComparison = figure( ...
    'Name','Itokawa common-resolution objective comparison', ...
    'Color','w', ...
    'Position',[50 80 1800 850]);

tlCommon = tiledlayout( ...
    figCommonComparison, ...
    2,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

rawTruthShow = simulatedItokawa.I/max(simulatedItokawa.I(:));
truthCommonShow = truthCommon/peakTruthCommon;
imageQPOCommonShow = imageQPOCommonScaled/peakTruthCommon;
imageToyCommonShow = imageToyCommonScaled/peakTruthCommon;
errorQPOCommonShow = abs(errorQPOCommon)/peakTruthCommon;
errorToyCommonShow = abs(errorToyCommon)/peakTruthCommon;

% Raw numerical truth
ax = nexttile(tlCommon);
imagesc(ax,simulatedItokawa.l_mas,simulatedItokawa.m_mas,rawTruthShow);
axis(ax,'image'); set(ax,'YDir','normal'); colormap(ax,positiveMap); clim(ax,[0 1]);
xlabel(ax,'$l$ [mas]'); ylabel(ax,'$m$ [mas]');
title(ax,'Raw rasterized truth'); colorbar(ax);

% Common-resolution truth
ax = nexttile(tlCommon);
imagesc(ax,simulatedItokawa.l_mas,simulatedItokawa.m_mas,truthCommonShow);
axis(ax,'image'); set(ax,'YDir','normal'); colormap(ax,positiveMap); clim(ax,[0 1]);
xlabel(ax,'$l$ [mas]'); ylabel(ax,'$m$ [mas]');
title(ax,sprintf('Truth at common %.3f mas FWHM',commonBeamFwhm_mas)); colorbar(ax);

% QPO common-resolution reconstruction
ax = nexttile(tlCommon);
imagesc(ax,simulatedItokawa.l_mas,simulatedItokawa.m_mas,imageQPOCommonShow);
axis(ax,'image'); set(ax,'YDir','normal'); colormap(ax,positiveMap); clim(ax,[0 1]);
xlabel(ax,'$l$ [mas]'); ylabel(ax,'$m$ [mas]');
title(ax,{ ...
    'QPO at common resolution', ...
    sprintf('NRMSE$_{\\rm supp}=%.3g$, NCC $=%.4f$, IoU$_{%.0f}=%.4f$', ...
    nrmseQPOCommonSupport,nccQPOCommon,100*contourLevel,iouQPOCommon)}); colorbar(ax);

% Toy common-resolution reconstruction
ax = nexttile(tlCommon);
imagesc(ax,simulatedItokawa.l_mas,simulatedItokawa.m_mas,imageToyCommonShow);
axis(ax,'image'); set(ax,'YDir','normal'); colormap(ax,positiveMap); clim(ax,[0 1]);
xlabel(ax,'$l$ [mas]'); ylabel(ax,'$m$ [mas]');
title(ax,{ ...
    'Toy at common resolution', ...
    sprintf('NRMSE$_{\\rm supp}=%.3g$, NCC $=%.4f$, IoU$_{%.0f}=%.4f$', ...
    nrmseToyCommonSupport,nccToyCommon,100*contourLevel,iouToyCommon)}); colorbar(ax);

% QPO absolute error
ax = nexttile(tlCommon);
imagesc(ax,simulatedItokawa.l_mas,simulatedItokawa.m_mas,errorQPOCommonShow);
axis(ax,'image'); set(ax,'YDir','normal'); colormap(ax,positiveMap);
xlabel(ax,'$l$ [mas]'); ylabel(ax,'$m$ [mas]');
title(ax,'QPO absolute error'); colorbar(ax);

% Toy absolute error
ax = nexttile(tlCommon);
imagesc(ax,simulatedItokawa.l_mas,simulatedItokawa.m_mas,errorToyCommonShow);
axis(ax,'image'); set(ax,'YDir','normal'); colormap(ax,positiveMap);
xlabel(ax,'$l$ [mas]'); ylabel(ax,'$m$ [mas]');
title(ax,'Toy absolute error'); colorbar(ax);

title(tlCommon,sprintf( ...
    ['Equal-resolution comparison; common circular FWHM $=%.3f$ mas; ', ...
     'native beam-area ratio QPO/toy $=%.2f$'], ...
    commonBeamFwhm_mas,beamAreaRatioQPOtoToy));

% ============================================================
% ============================================================
%% 15. ROTATING ITOKAWA TARGET — ADDITIVE FORWARD-MODEL DIAGNOSTIC
% This section leaves the complete static target study above unchanged. It adds only the known principal-axis rotation and evaluates the existing Gaskell projection/rasterization chain at multiple rotation phases.
% -------------------------------------------------------------------------
% ADDITIVE ROTATING-TARGET DIAGNOSTIC
%
% This section does NOT overwrite the static target or any of the existing
% QPO/toy observations and reconstructions above. It only adds a rotating
% Gaskell-Itokawa target model for later time-dependent / stroboscopic work.
%
% Convention:
%
%   r_I(t) = R_IB(t) r_B
%
% with the existing static attitude R_IB_Itokawa used as phase zero.
% The Gaskell body +z_B axis is the rotation pole.
% -------------------------------------------------------------------------

rotItokawa = struct();
rotItokawa.period_s   = 12.132381791359627*3600;
rotItokawa.spinAxis_B = [0;0;1];
rotItokawa.R_IB0      = R_IB_Itokawa;
rotItokawa.t0_s       = 0.0;
rotItokawa.spinSense  = +1;

% Relative phase convention only; this does not claim an absolute UTC/PCK
% inertial phase or spin-sense convention for Itokawa.

% The QPO propagation time geomQPO.t is nondimensional CR3BP time.
% Convert elapsed propagation time to physical seconds using t*.
qpoElapsedRot_s = ...
    (geomQPO.t(:) - geomQPO.t(1)) * ...
    sys.Tstar_days * 86400;

nRotEpochs = numel(qpoElapsedRot_s);

rotationPhaseQPO_rad = zeros(nRotEpochs,1);
rotationPhaseQPO_unwrapped_rad = zeros(nRotEpochs,1);

for kRot = 1:nRotEpochs

    [~,rotationPhaseQPO_rad(kRot),rotationPhaseQPO_unwrapped_rad(kRot)] = ...
        asteroid_attitude_at_time( ...
            rotItokawa, ...
            qpoElapsedRot_s(kRot));

end

% Store the inexpensive rotational state at every QPO epoch. These arrays
% will later be used directly for phase-bin / stroboscopic sample selection.
rotatingItokawa = struct();
rotatingItokawa.model = rotItokawa;
rotatingItokawa.qpoElapsed_s = qpoElapsedRot_s;
rotatingItokawa.qpoPhase_rad = rotationPhaseQPO_rad;
rotatingItokawa.qpoPhase_deg = rad2deg(rotationPhaseQPO_rad);
rotatingItokawa.qpoPhaseUnwrapped_rad = rotationPhaseQPO_unwrapped_rad;

% One-period attitude closure check.
[R_IB_rot0,~,~] = asteroid_attitude_at_time(rotItokawa,0.0);
[R_IB_rotT,~,~] = asteroid_attitude_at_time(rotItokawa,rotItokawa.period_s);

rotationClosureError = norm(R_IB_rotT-R_IB_rot0,'fro');
initialAttitudeError = norm(R_IB_rot0-R_IB_Itokawa,'fro');

fprintf('\n============================================================\n');
fprintf('ROTATING GASKELL ITOKAWA TARGET\n');
fprintf('============================================================\n');
fprintf('  rotation period             = %.9f h\n',rotItokawa.period_s/3600);
fprintf('  body-frame spin axis        = [%+.1f %+.1f %+.1f]\n', ...
    rotItokawa.spinAxis_B);
fprintf('  QPO epochs                  = %d\n',nRotEpochs);
fprintf('  QPO elapsed duration        = %.6f days\n',qpoElapsedRot_s(end)/86400);
fprintf('  asteroid rotations covered  = %.6f\n', ...
    qpoElapsedRot_s(end)/rotItokawa.period_s);
fprintf('  initial attitude error      = %.3e\n',initialAttitudeError);
fprintf('  one-period closure error    = %.3e\n',rotationClosureError);
fprintf('============================================================\n');

% Rotation phase sampled at the SAME QPO epochs used by the formation.
figItokawaRotationPhase = figure( ...
    'Name','Itokawa rotation phase over QPO propagation', ...
    'Color','w');

plot( ...
    qpoElapsedRot_s/86400, ...
    rad2deg(rotationPhaseQPO_rad), ...
    'LineWidth',1.4);

grid on;
xlabel('Elapsed physical time [days]');
ylabel('Wrapped Itokawa rotation phase [deg]');
ylim([0 360]);
title('Itokawa rotation phase at QPO propagation epochs');

% -------------------------------------------------------------------------
% Time-dependent Gaskell forward-model snapshots.
%
% Use one FIXED physical raster half-width for every attitude. This is
% essential: changing the raster extent with rotation would change the image
% coordinates themselves and would contaminate later visibility comparisons.
% The maximum body-frame radius safely bounds every orthographic projection.
% -------------------------------------------------------------------------

nRotationSnapshots = 8;
rotationSnapshotTimes_s = ...
    (0:nRotationSnapshots-1).' / nRotationSnapshots * ...
    rotItokawa.period_s;

bodyRadiusItokawa_km = ...
    max(vecnorm(shapeItokawa.vertices,2,2));

halfWidthRotatingItokawa_km = ...
    marginItokawa * bodyRadiusItokawa_km;

rotationSnapshotPhase_rad = zeros(nRotationSnapshots,1);
rotationSnapshotSimulated = cell(nRotationSnapshots,1);
rotationSnapshotVisibleFacets = zeros(nRotationSnapshots,1);

figItokawaRotationSnapshots = figure( ...
    'Name','Rotating Gaskell Itokawa projected truth', ...
    'Color','w', ...
    'Position',[70 80 1500 760]);

tlRot = tiledlayout( ...
    figItokawaRotationSnapshots, ...
    2,4, ...
    'TileSpacing','compact', ...
    'Padding','compact');

for kRot = 1:nRotationSnapshots

    tRot_s = rotationSnapshotTimes_s(kRot);

    [R_IB_rot,phiRot] = ...
        asteroid_attitude_at_time( ...
            rotItokawa, ...
            tRot_s);

    rotationSnapshotPhase_rad(kRot) = phiRot;

    projRot = ...
        project_shape_orthographic( ...
            shapeItokawa, ...
            R_IB_rot, ...
            target);

    visRot = ...
        find_visible_facets( ...
            projRot);

    imgRot = ...
        rasterize_visible_facets( ...
            projRot, ...
            visRot, ...
            nPixItokawa, ...
            halfWidthRotatingItokawa_km, ...
            facetRadianceItokawa);

    simulatedRot = ...
        make_simulated_target_from_raster( ...
            imgRot, ...
            rangeItokawa_km, ...
            shapeItokawa.type);

    rotationSnapshotSimulated{kRot} = simulatedRot;
    rotationSnapshotVisibleFacets(kRot) = nnz(visRot.isVisible);

    axRot = nexttile(tlRot);

    peakRot = max(simulatedRot.I(:));

    if peakRot <= 0
        error('Rotating Itokawa snapshot %d has no positive pixels.',kRot);
    end

    imagesc( ...
        axRot, ...
        simulatedRot.l_mas, ...
        simulatedRot.m_mas, ...
        simulatedRot.I/peakRot);

    axis(axRot,'image');
    set(axRot,'YDir','normal');
    colormap(axRot,hot(256));
    clim(axRot,[0 1]);

    xlabel(axRot,'$l$ [mas]');
    ylabel(axRot,'$m$ [mas]');

    title( ...
        axRot, ...
        sprintf('phase = %.0f deg',rad2deg(phiRot)), ...
        'Interpreter','latex');

end

cbRot = colorbar;
cbRot.Layout.Tile = 'east';
cbRot.Label.String = 'Relative intensity';

title( ...
    tlRot, ...
    'Gaskell Itokawa: same target frame and fixed raster grid through one rotation');

rotatingItokawa.snapshotTimes_s = rotationSnapshotTimes_s;
rotatingItokawa.snapshotPhase_rad = rotationSnapshotPhase_rad;
rotatingItokawa.snapshotPhase_deg = rad2deg(rotationSnapshotPhase_rad);
rotatingItokawa.snapshotSimulated = rotationSnapshotSimulated;
rotatingItokawa.snapshotVisibleFacets = rotationSnapshotVisibleFacets;
rotatingItokawa.halfWidth_km = halfWidthRotatingItokawa_km;

fprintf('\nRotating-target forward-model snapshots\n');
fprintf('  snapshots                   = %d\n',nRotationSnapshots);
fprintf('  fixed raster half-width     = %.6f km\n',halfWidthRotatingItokawa_km);
fprintf('  visible facets range        = %d to %d\n', ...
    min(rotationSnapshotVisibleFacets), ...
    max(rotationSnapshotVisibleFacets));


%% 16. TIME-DEPENDENT ROTATING-ITOKAWA VISIBILITIES AT THE ACTUAL QPO EPOCHS
% ============================================================
%
% For each accessible QPO epoch:
%
%   t_k
%       -> R_IB(t_k)
%       -> projected Gaskell mesh
%       -> visible facets
%       -> fixed-grid raster I_k(l,m)
%       -> V_k(u_k,v_k)
%
% The existing static obsQPO structure is not modified.
% A new structure, obsRotatingQPO, is created.
%
% All rotating images use the same fixed raster grid established
% in Section 15.


%% Accessible QPO epochs

accessibleEpochRot = find(accessQPO.accessible(:));

if ~isequal(accessibleEpochRot(:),obsQPO.accessibleIndex(:))
    error(['Rotating-target epoch bookkeeping disagrees with obsQPO. ', ...
        'The accessible-index ordering must be identical.']);
end

nAccessibleRot = numel(accessibleEpochRot);
nBaselinesRot = numel(baselineFields);

if nBaselinesRot ~= 3
    error('Section 16 expects exactly three baselines: B12, B13, and B23.');
end


%% Extract physical u-v-w samples for each accessible epoch

uRot = zeros(nAccessibleRot,nBaselinesRot);
vRot = zeros(nAccessibleRot,nBaselinesRot);
wRot = zeros(nAccessibleRot,nBaselinesRot);

for bRot = 1:nBaselinesRot

    uvwAllRot = uvQPO.(baselineFields{bRot});

    if size(uvwAllRot,1) ~= numel(accessQPO.accessible) || size(uvwAllRot,2) < 3
        error('Unexpected size for uvQPO.%s.',baselineFields{bRot});
    end

    uvwAccessibleRot = uvwAllRot(accessibleEpochRot,1:3);

    uRot(:,bRot) = uvwAccessibleRot(:,1);
    vRot(:,bRot) = uvwAccessibleRot(:,2);
    wRot(:,bRot) = uvwAccessibleRot(:,3);

end


%% Determine the sample ordering used by obsQPO
%
% Possible orderings:
%
% baseline-major:
%   all B12 epochs, then all B13 epochs, then all B23 epochs
%
% epoch-major:
%   B12, B13, B23 at epoch 1, then epoch 2, etc.

uBaselineMajor = uRot(:);
vBaselineMajor = vRot(:);

uEpochMajor = reshape(uRot.',[],1);
vEpochMajor = reshape(vRot.',[],1);

uObsStatic = obsQPO.u(:);
vObsStatic = obsQPO.v(:);

expectedSampleCount = nAccessibleRot*nBaselinesRot;

if numel(uObsStatic) ~= expectedSampleCount
    error(['obsQPO sample count does not equal accessible epochs ', ...
        'times the number of baselines.']);
end

orderErrorBaselineMajor = max(hypot( ...
    uBaselineMajor-uObsStatic, ...
    vBaselineMajor-vObsStatic));

orderErrorEpochMajor = max(hypot( ...
    uEpochMajor-uObsStatic, ...
    vEpochMajor-vObsStatic));

uvScaleOrder = max(1,max(abs([uObsStatic;vObsStatic])));
orderTolerance = 1e-12*uvScaleOrder;

if orderErrorBaselineMajor <= orderErrorEpochMajor
    rotatingSampleOrder = 'baseline-major';
    orderErrorSelected = orderErrorBaselineMajor;
else
    rotatingSampleOrder = 'epoch-major';
    orderErrorSelected = orderErrorEpochMajor;
end

if orderErrorSelected > orderTolerance
    error(['Could not reproduce the sample ordering in obsQPO from uvQPO. ', ...
        'Best u-v mismatch = %.6e wavelengths.'],orderErrorSelected);
end

fprintf('\n============================================================\n');
fprintf('TIME-DEPENDENT ROTATING ITOKAWA VISIBILITIES\n');
fprintf('============================================================\n');
fprintf('  accessible QPO epochs       = %d\n',nAccessibleRot);
fprintf('  baselines per epoch         = %d\n',nBaselinesRot);
fprintf('  total physical samples      = %d\n',expectedSampleCount);
fprintf('  detected sample ordering    = %s\n',rotatingSampleOrder);
fprintf('  ordering u-v mismatch       = %.3e wavelengths\n',orderErrorSelected);


%% Audit direct NUDFT against the existing static production result
%
% This uses the original static simulatedItokawa image.
% If this does not reproduce obsQPO.V, stop before doing any
% rotating-target calculations.

[jStaticPix,iStaticPix,fluxStaticPix] = find(simulatedItokawa.I);

lStaticPix_rad = simulatedItokawa.l_rad(iStaticPix);
mStaticPix_rad = simulatedItokawa.m_rad(jStaticPix);

lStaticPix_rad = lStaticPix_rad(:);
mStaticPix_rad = mStaticPix_rad(:);
fluxStaticPix = fluxStaticPix(:);

assert(iscolumn(lStaticPix_rad), ...
    'lStaticPix_rad must be a column vector.');

assert(iscolumn(mStaticPix_rad), ...
    'mStaticPix_rad must be a column vector.');

assert(iscolumn(fluxStaticPix), ...
    'fluxStaticPix must be a column vector.');

assert(numel(lStaticPix_rad) == numel(fluxStaticPix), ...
    'Static l-coordinate and flux vectors have inconsistent lengths.');

assert(numel(mStaticPix_rad) == numel(fluxStaticPix), ...
    'Static m-coordinate and flux vectors have inconsistent lengths.');

VStaticManual = complex(zeros(nAccessibleRot,nBaselinesRot));

for bRot = 1:nBaselinesRot

    for kAcc = 1:nAccessibleRot

        phaseStatic = ...
            uRot(kAcc,bRot)*lStaticPix_rad + ...
            vRot(kAcc,bRot)*mStaticPix_rad;

        VStaticManual(kAcc,bRot) = ...
            sum(fluxStaticPix .* exp(-2*pi*1i*phaseStatic));

    end

end


%% Flatten the manually computed static visibility in production order

switch rotatingSampleOrder

    case 'baseline-major'
        VStaticManualVector = VStaticManual(:);

    case 'epoch-major'
        VStaticManualVector = reshape(VStaticManual.',[],1);

    otherwise
        error('Unexpected rotatingSampleOrder.');

end

staticNUDFTError = VStaticManualVector-obsQPO.V(:);

staticNUDFTMaxError = max(abs(staticNUDFTError));
staticNUDFTRMSError = sqrt(mean(abs(staticNUDFTError).^2));

fprintf('  static NUDFT max error      = %.3e\n',staticNUDFTMaxError);
fprintf('  static NUDFT RMS error      = %.3e\n',staticNUDFTRMSError);

if staticNUDFTMaxError > 1e-10
    error(['Direct NUDFT does not reproduce the existing static ', ...
        'obsQPO visibility data. Stop before generating rotating data.']);
end


%% Frozen phase-zero target on the fixed rotating-target image grid
%
% This control uses the same fixed raster grid as every rotating image.
% Therefore differences between this target and the rotating target are
% caused by attitude change rather than changing raster dimensions.

[R_IB_frozenRot,~,~] = ...
    asteroid_attitude_at_time(rotItokawa,rotItokawa.t0_s);

projFrozenRot = ...
    project_shape_orthographic( ...
        shapeItokawa, ...
        R_IB_frozenRot, ...
        target);

visFrozenRot = find_visible_facets(projFrozenRot);

imgFrozenRot = ...
    rasterize_visible_facets( ...
        projFrozenRot, ...
        visFrozenRot, ...
        nPixItokawa, ...
        halfWidthRotatingItokawa_km, ...
        facetRadianceItokawa);

simulatedFrozenRot = ...
    make_simulated_target_from_raster( ...
        imgFrozenRot, ...
        rangeItokawa_km, ...
        shapeItokawa.type);


%% Nonzero pixels of the frozen phase-zero image

[jFrozenPix,iFrozenPix,fluxFrozenPix] = find(simulatedFrozenRot.I);

lFrozenPix_rad = simulatedFrozenRot.l_rad(iFrozenPix);
mFrozenPix_rad = simulatedFrozenRot.m_rad(jFrozenPix);

lFrozenPix_rad = lFrozenPix_rad(:);
mFrozenPix_rad = mFrozenPix_rad(:);
fluxFrozenPix = fluxFrozenPix(:);

assert(iscolumn(lFrozenPix_rad), ...
    'lFrozenPix_rad must be a column vector.');

assert(iscolumn(mFrozenPix_rad), ...
    'mFrozenPix_rad must be a column vector.');

assert(iscolumn(fluxFrozenPix), ...
    'fluxFrozenPix must be a column vector.');

assert(numel(lFrozenPix_rad) == numel(fluxFrozenPix), ...
    'Frozen l-coordinate and flux vectors have inconsistent lengths.');

assert(numel(mFrozenPix_rad) == numel(fluxFrozenPix), ...
    'Frozen m-coordinate and flux vectors have inconsistent lengths.');


%% Frozen-target visibilities on the actual QPO u-v samples

VFrozenFixedGrid = complex(zeros(nAccessibleRot,nBaselinesRot));

for bRot = 1:nBaselinesRot

    for kAcc = 1:nAccessibleRot

        phaseFrozen = ...
            uRot(kAcc,bRot)*lFrozenPix_rad + ...
            vRot(kAcc,bRot)*mFrozenPix_rad;

        VFrozenFixedGrid(kAcc,bRot) = ...
            sum(fluxFrozenPix .* exp(-2*pi*1i*phaseFrozen));

    end

end


%% Allocate time-dependent rotating-target output

VRotatingQPO = complex(zeros(nAccessibleRot,nBaselinesRot));

rotationPhaseAccessible_rad = zeros(nAccessibleRot,1);
rotationPhaseAccessible_deg = zeros(nAccessibleRot,1);

elapsedAccessible_s = qpoElapsedRot_s(accessibleEpochRot);

visibleFacetCountRot = zeros(nAccessibleRot,1);
nonzeroPixelCountRot = zeros(nAccessibleRot,1);
fluxNormalizationErrorRot = zeros(nAccessibleRot,1);

lRotatingGrid_rad = [];
mRotatingGrid_rad = [];

lRotatingGrid_mas = [];
mRotatingGrid_mas = [];


%% Generate rotating target at every accessible QPO epoch

fprintf('\nGenerating rotating target visibilities...\n');

ticRotatingVis = tic;

for kAcc = 1:nAccessibleRot

    epochIndex = accessibleEpochRot(kAcc);
    tRot_s = qpoElapsedRot_s(epochIndex);

    [R_IB_rot,phiRot,~] = ...
        asteroid_attitude_at_time(rotItokawa,tRot_s);

    rotationPhaseAccessible_rad(kAcc) = phiRot;
    rotationPhaseAccessible_deg(kAcc) = rad2deg(phiRot);


    %% Project the rotating Gaskell mesh

    projRotEpoch = ...
        project_shape_orthographic( ...
            shapeItokawa, ...
            R_IB_rot, ...
            target);

    visRotEpoch = find_visible_facets(projRotEpoch);


    %% Rasterize on the same fixed grid for every epoch

    imgRotEpoch = ...
        rasterize_visible_facets( ...
            projRotEpoch, ...
            visRotEpoch, ...
            nPixItokawa, ...
            halfWidthRotatingItokawa_km, ...
            facetRadianceItokawa);

    simulatedRotEpoch = ...
        make_simulated_target_from_raster( ...
            imgRotEpoch, ...
            rangeItokawa_km, ...
            shapeItokawa.type);


    %% Verify that the angular image grid is fixed

    if kAcc == 1

        lRotatingGrid_rad = simulatedRotEpoch.l_rad(:);
        mRotatingGrid_rad = simulatedRotEpoch.m_rad(:);

        lRotatingGrid_mas = simulatedRotEpoch.l_mas(:);
        mRotatingGrid_mas = simulatedRotEpoch.m_mas(:);

    else

        lGridError = max(abs( ...
            simulatedRotEpoch.l_rad(:)-lRotatingGrid_rad));

        mGridError = max(abs( ...
            simulatedRotEpoch.m_rad(:)-mRotatingGrid_rad));

        gridErrorRot = max(lGridError,mGridError);

        if gridErrorRot > 1e-15
            error('Rotating target angular grid changed at QPO epoch %d.', ...
                epochIndex);
        end

    end


    %% Extract nonzero image pixels

    [jPixRot,iPixRot,fluxPixRot] = find(simulatedRotEpoch.I);

    if isempty(fluxPixRot)
        error('Rotating target has no nonzero pixels at QPO epoch %d.', ...
            epochIndex);
    end

    lPixRot_rad = simulatedRotEpoch.l_rad(iPixRot);
    mPixRot_rad = simulatedRotEpoch.m_rad(jPixRot);

    lPixRot_rad = lPixRot_rad(:);
    mPixRot_rad = mPixRot_rad(:);
    fluxPixRot = fluxPixRot(:);

    assert(iscolumn(lPixRot_rad), ...
        'lPixRot_rad must be a column vector.');

    assert(iscolumn(mPixRot_rad), ...
        'mPixRot_rad must be a column vector.');

    assert(iscolumn(fluxPixRot), ...
        'fluxPixRot must be a column vector.');

    assert(numel(lPixRot_rad) == numel(fluxPixRot), ...
        'Rotating l-coordinate and flux vectors have inconsistent lengths.');

    assert(numel(mPixRot_rad) == numel(fluxPixRot), ...
        'Rotating m-coordinate and flux vectors have inconsistent lengths.');


    %% Epoch diagnostics

    visibleFacetCountRot(kAcc) = nnz(visRotEpoch.isVisible);
    nonzeroPixelCountRot(kAcc) = numel(fluxPixRot);
    fluxNormalizationErrorRot(kAcc) = abs(sum(fluxPixRot)-1);


    %% Evaluate B12, B13, and B23 at this epoch

    for bRot = 1:nBaselinesRot

        phaseRot = ...
            uRot(kAcc,bRot)*lPixRot_rad + ...
            vRot(kAcc,bRot)*mPixRot_rad;

        VRotatingQPO(kAcc,bRot) = ...
            sum(fluxPixRot .* exp(-2*pi*1i*phaseRot));

    end


    %% Progress output

    if kAcc == 1 || mod(kAcc,25) == 0 || kAcc == nAccessibleRot

        fprintf(['  rotating epoch %4d / %4d: ', ...
            'phase = %7.2f deg, visible facets = %d\n'], ...
            kAcc, ...
            nAccessibleRot, ...
            rotationPhaseAccessible_deg(kAcc), ...
            visibleFacetCountRot(kAcc));

    end

end

rotatingVisibilityRuntime_s = toc(ticRotatingVis);


%% Flatten rotating and frozen visibility arrays in production order

switch rotatingSampleOrder

    case 'baseline-major'

        VRotatingVector = VRotatingQPO(:);
        VFrozenVector = VFrozenFixedGrid(:);

        sampleEpochIndexRot = ...
            repmat(accessibleEpochRot,nBaselinesRot,1);

        sampleBaselineIndexRot = [ ...
            ones(nAccessibleRot,1); ...
            2*ones(nAccessibleRot,1); ...
            3*ones(nAccessibleRot,1)];

        sampleElapsed_sRot = ...
            repmat(elapsedAccessible_s,nBaselinesRot,1);

        sampleRotationPhase_radRot = ...
            repmat(rotationPhaseAccessible_rad,nBaselinesRot,1);


    case 'epoch-major'

        VRotatingVector = reshape(VRotatingQPO.',[],1);
        VFrozenVector = reshape(VFrozenFixedGrid.',[],1);

        sampleEpochIndexRot = ...
            repelem(accessibleEpochRot,nBaselinesRot);

        sampleBaselineIndexRot = ...
            repmat((1:nBaselinesRot).',nAccessibleRot,1);

        sampleElapsed_sRot = ...
            repelem(elapsedAccessible_s,nBaselinesRot);

        sampleRotationPhase_radRot = ...
            repelem(rotationPhaseAccessible_rad,nBaselinesRot);


    otherwise

        error('Unexpected rotatingSampleOrder.');

end


%% Final bookkeeping checks

assert(numel(VRotatingVector) == numel(obsQPO.V), ...
    'Rotating visibility vector has incorrect length.');

assert(numel(VFrozenVector) == numel(obsQPO.V), ...
    'Frozen visibility vector has incorrect length.');

assert(numel(sampleEpochIndexRot) == numel(obsQPO.V), ...
    'Sample epoch-index vector has incorrect length.');

assert(numel(sampleBaselineIndexRot) == numel(obsQPO.V), ...
    'Sample baseline-index vector has incorrect length.');

assert(numel(sampleElapsed_sRot) == numel(obsQPO.V), ...
    'Sample time vector has incorrect length.');

assert(numel(sampleRotationPhase_radRot) == numel(obsQPO.V), ...
    'Sample rotation-phase vector has incorrect length.');


%% Create observation structure compatible with the existing imaging code

obsRotatingQPO = obsQPO;

obsRotatingQPO.V = ...
    reshape(VRotatingVector,size(obsQPO.V));

obsRotatingQPO.VFrozenPhaseZero = ...
    reshape(VFrozenVector,size(obsQPO.V));


%% Update top-level derived visibility fields if they exist

if isfield(obsRotatingQPO,'amplitude') && ...
        numel(obsRotatingQPO.amplitude) == numel(VRotatingVector)

    obsRotatingQPO.amplitude = ...
        reshape(abs(VRotatingVector),size(obsRotatingQPO.amplitude));

end

if isfield(obsRotatingQPO,'phase') && ...
        numel(obsRotatingQPO.phase) == numel(VRotatingVector)

    obsRotatingQPO.phase = ...
        reshape(angle(VRotatingVector),size(obsRotatingQPO.phase));

end

if isfield(obsRotatingQPO,'phase_rad') && ...
        numel(obsRotatingQPO.phase_rad) == numel(VRotatingVector)

    obsRotatingQPO.phase_rad = ...
        reshape(angle(VRotatingVector),size(obsRotatingQPO.phase_rad));

end

if isfield(obsRotatingQPO,'phase_deg') && ...
        numel(obsRotatingQPO.phase_deg) == numel(VRotatingVector)

    obsRotatingQPO.phase_deg = ...
        reshape(rad2deg(angle(VRotatingVector)), ...
        size(obsRotatingQPO.phase_deg));

end


%% Update baseline-resolved visibility fields if they exist

if isfield(obsRotatingQPO,'baselines')

    for bRot = 1:min(nBaselinesRot,numel(obsRotatingQPO.baselines))

        if isfield(obsRotatingQPO.baselines(bRot),'V') && ...
                numel(obsRotatingQPO.baselines(bRot).V) == nAccessibleRot

            oldSize = size(obsRotatingQPO.baselines(bRot).V);

            obsRotatingQPO.baselines(bRot).V = ...
                reshape(VRotatingQPO(:,bRot),oldSize);

        end

        if isfield(obsRotatingQPO.baselines(bRot),'amplitude') && ...
                numel(obsRotatingQPO.baselines(bRot).amplitude) == nAccessibleRot

            oldSize = size(obsRotatingQPO.baselines(bRot).amplitude);

            obsRotatingQPO.baselines(bRot).amplitude = ...
                reshape(abs(VRotatingQPO(:,bRot)),oldSize);

        end

        if isfield(obsRotatingQPO.baselines(bRot),'phase') && ...
                numel(obsRotatingQPO.baselines(bRot).phase) == nAccessibleRot

            oldSize = size(obsRotatingQPO.baselines(bRot).phase);

            obsRotatingQPO.baselines(bRot).phase = ...
                reshape(angle(VRotatingQPO(:,bRot)),oldSize);

        end

        if isfield(obsRotatingQPO.baselines(bRot),'phase_rad') && ...
                numel(obsRotatingQPO.baselines(bRot).phase_rad) == nAccessibleRot

            oldSize = size(obsRotatingQPO.baselines(bRot).phase_rad);

            obsRotatingQPO.baselines(bRot).phase_rad = ...
                reshape(angle(VRotatingQPO(:,bRot)),oldSize);

        end

        if isfield(obsRotatingQPO.baselines(bRot),'phase_deg') && ...
                numel(obsRotatingQPO.baselines(bRot).phase_deg) == nAccessibleRot

            oldSize = size(obsRotatingQPO.baselines(bRot).phase_deg);

            obsRotatingQPO.baselines(bRot).phase_deg = ...
                reshape(rad2deg(angle(VRotatingQPO(:,bRot))),oldSize);

        end

    end

end


%% Add explicit time and rotation metadata

obsRotatingQPO.epochIndex = sampleEpochIndexRot;
obsRotatingQPO.baselineIndex = sampleBaselineIndexRot;

obsRotatingQPO.elapsed_s = sampleElapsed_sRot;

obsRotatingQPO.rotationPhase_rad = sampleRotationPhase_radRot;
obsRotatingQPO.rotationPhase_deg = ...
    rad2deg(sampleRotationPhase_radRot);

obsRotatingQPO.rotationPeriod_s = rotItokawa.period_s;
obsRotatingQPO.sampleOrder = rotatingSampleOrder;


%% Store transparent epoch-resolved rotating-target data

rotatingQPOData = struct();

rotatingQPOData.accessibleEpochIndex = accessibleEpochRot;

rotatingQPOData.elapsed_s = elapsedAccessible_s;
rotatingQPOData.elapsed_days = elapsedAccessible_s/86400;

rotatingQPOData.rotationPhase_rad = rotationPhaseAccessible_rad;
rotatingQPOData.rotationPhase_deg = rotationPhaseAccessible_deg;

rotatingQPOData.u = uRot;
rotatingQPOData.v = vRot;
rotatingQPOData.w = wRot;

rotatingQPOData.V = VRotatingQPO;
rotatingQPOData.VFrozenPhaseZero = VFrozenFixedGrid;

rotatingQPOData.visibleFacetCount = visibleFacetCountRot;
rotatingQPOData.nonzeroPixelCount = nonzeroPixelCountRot;
rotatingQPOData.fluxNormalizationError = fluxNormalizationErrorRot;

rotatingQPOData.l_rad = lRotatingGrid_rad;
rotatingQPOData.m_rad = mRotatingGrid_rad;

rotatingQPOData.l_mas = lRotatingGrid_mas;
rotatingQPOData.m_mas = mRotatingGrid_mas;

rotatingQPOData.simulatedPhaseZero = simulatedFrozenRot;
rotatingQPOData.sampleOrder = rotatingSampleOrder;


%% Summary diagnostics

maxFluxNormalizationErrorRot = ...
    max(fluxNormalizationErrorRot);

visibilityDifferenceRot = ...
    VRotatingVector-VFrozenVector;

maxMotionVisibilityDifference = ...
    max(abs(visibilityDifferenceRot));

rmsMotionVisibilityDifference = ...
    sqrt(mean(abs(visibilityDifferenceRot).^2));

fprintf('\nROTATING-QPO DATA SET COMPLETE\n');
fprintf('  runtime                      = %.3f s\n', ...
    rotatingVisibilityRuntime_s);
fprintf('  max flux normalization error = %.3e\n', ...
    maxFluxNormalizationErrorRot);
fprintf('  visible facets range         = %d to %d\n', ...
    min(visibleFacetCountRot),max(visibleFacetCountRot));
fprintf('  nonzero pixel range          = %d to %d\n', ...
    min(nonzeroPixelCountRot),max(nonzeroPixelCountRot));
fprintf('  max |Vrot-Vfrozen|           = %.6e\n', ...
    maxMotionVisibilityDifference);
fprintf('  RMS |Vrot-Vfrozen|           = %.6e\n', ...
    rmsMotionVisibilityDifference);
fprintf('============================================================\n');


%% Diagnostic figure: visibility versus asteroid rotation phase

phaseAmpThresholdRot = 1e-5;

figRotatingVisibility = figure( ...
    'Name','Time-dependent Itokawa visibility versus rotation phase', ...
    'Color','w', ...
    'Position',[100 100 1200 800]);

tlRotatingVisibility = tiledlayout( ...
    figRotatingVisibility, ...
    2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% Visibility amplitude

axAmpRot = nexttile(tlRotatingVisibility);

hold(axAmpRot,'on');
grid(axAmpRot,'on');
box(axAmpRot,'on');

for bRot = 1:nBaselinesRot

    plot( ...
        axAmpRot, ...
        rotationPhaseAccessible_deg, ...
        abs(VRotatingQPO(:,bRot)), ...
        '.', ...
        'MarkerSize',9, ...
        'DisplayName',baselineNames{bRot});

end

xlabel(axAmpRot,'Itokawa rotation phase [deg]');
ylabel(axAmpRot,'|V|');

xlim(axAmpRot,[0 360]);

legend(axAmpRot,'Location','best');

title(axAmpRot,'Visibility amplitude at the actual QPO samples');


%% Visibility phase

axPhaseRot = nexttile(tlRotatingVisibility);

hold(axPhaseRot,'on');
grid(axPhaseRot,'on');
box(axPhaseRot,'on');

for bRot = 1:nBaselinesRot

    phaseMaskRot = ...
        abs(VRotatingQPO(:,bRot)) > phaseAmpThresholdRot;

    plot( ...
        axPhaseRot, ...
        rotationPhaseAccessible_deg(phaseMaskRot), ...
        rad2deg(angle(VRotatingQPO(phaseMaskRot,bRot))), ...
        '.', ...
        'MarkerSize',9, ...
        'DisplayName',baselineNames{bRot});

end

xlabel(axPhaseRot,'Itokawa rotation phase [deg]');
ylabel(axPhaseRot,'Visibility phase [deg]');

xlim(axPhaseRot,[0 360]);
ylim(axPhaseRot,[-180 180]);

legend(axPhaseRot,'Location','best');

title(axPhaseRot,'Visibility phase at the actual QPO samples');

title(tlRotatingVisibility, ...
    'Time-dependent visibility of the rotating Gaskell Itokawa target');


%% Diagnostic figure: rotating versus frozen target

figRotatingDifference = figure( ...
    'Name','Rotating versus frozen Itokawa visibility difference', ...
    'Color','w', ...
    'Position',[120 120 1000 650]);

axDifferenceRot = axes(figRotatingDifference);

hold(axDifferenceRot,'on');
grid(axDifferenceRot,'on');
box(axDifferenceRot,'on');

for bRot = 1:nBaselinesRot

    plot( ...
        axDifferenceRot, ...
        rotationPhaseAccessible_deg, ...
        abs(VRotatingQPO(:,bRot)-VFrozenFixedGrid(:,bRot)), ...
        '.', ...
        'MarkerSize',9, ...
        'DisplayName',baselineNames{bRot});

end

xlabel(axDifferenceRot,'Itokawa rotation phase [deg]');
ylabel(axDifferenceRot,'|Vrot - Vfrozen|');

xlim(axDifferenceRot,[0 360]);

legend(axDifferenceRot,'Location','best');

title(axDifferenceRot, ...
    'Visibility inconsistency caused by treating the rotating target as static');
    %% 


 %% 17. DISENTANGLE ROTATION PHASE FROM BASELINE SPATIAL FREQUENCY
% ============================================================
%
% The quantity
%
%   |V_rot - V_frozen|
%
% depends on both the asteroid attitude and the sampled spatial frequency.
%
% This section examines the same data in two complementary ways:
%
%   1. visibility mismatch versus rho, colored by asteroid phase;
%   2. visibility mismatch versus asteroid phase, colored by rho.
%
% This does not yet remove the effect of baseline orientation in the u-v
% plane. It isolates baseline LENGTH first.


%% Spatial frequency magnitude

rhoRot = hypot(uRot,vRot);
rhoRot_Glambda = rhoRot/1e9;

deltaVComplex = VRotatingQPO-VFrozenFixedGrid;
deltaVAbs = abs(deltaVComplex);

deltaAmplitude = ...
    abs(abs(VRotatingQPO)-abs(VFrozenFixedGrid));

phaseDifference_deg = ...
    rad2deg(angle(VRotatingQPO.*conj(VFrozenFixedGrid)));


%% Basic numerical summary by baseline

fprintf('\n');
fprintf('============================================================\n');
fprintf('ROTATION SENSITIVITY VERSUS SPATIAL FREQUENCY\n');
fprintf('============================================================\n');

for bRot = 1:nBaselinesRot

    rhoThis = rhoRot_Glambda(:,bRot);
    deltaThis = deltaVAbs(:,bRot);

    C = corrcoef(rhoThis,deltaThis);

    if numel(C) == 4
        rhoDeltaCorrelation = C(1,2);
    else
        rhoDeltaCorrelation = NaN;
    end

    fprintf('\n%s\n',baselineNames{bRot});
    fprintf('  rho minimum              = %.6f Glambda\n',min(rhoThis));
    fprintf('  rho median               = %.6f Glambda\n',median(rhoThis));
    fprintf('  rho maximum              = %.6f Glambda\n',max(rhoThis));
    fprintf('  median |Delta V|         = %.6f\n',median(deltaThis));
    fprintf('  maximum |Delta V|        = %.6f\n',max(deltaThis));
    fprintf('  corr(rho,|Delta V|)      = %+.6f\n',rhoDeltaCorrelation);

end

fprintf('\n============================================================\n');


%% Figure 1: mismatch versus spatial frequency, colored by rotation phase

figDeltaVvsRho = figure( ...
    'Name','Rotating-target mismatch versus spatial frequency', ...
    'Color','w', ...
    'Position',[80 100 1500 520]);

tlDeltaVvsRho = tiledlayout( ...
    figDeltaVvsRho, ...
    1,nBaselinesRot, ...
    'TileSpacing','compact', ...
    'Padding','compact');

nRhoBins = 12;

for bRot = 1:nBaselinesRot

    ax = nexttile(tlDeltaVvsRho);

    hold(ax,'on');
    grid(ax,'on');
    box(ax,'on');

    rhoThis = rhoRot_Glambda(:,bRot);
    deltaThis = deltaVAbs(:,bRot);
    phaseThis = rotationPhaseAccessible_deg;

    scatter( ...
        ax, ...
        rhoThis, ...
        deltaThis, ...
        22, ...
        phaseThis, ...
        'filled');

    clim(ax,[0 360]);

    xlabel(ax,'\rho [G\lambda]');
    ylabel(ax,'|V_{rot}-V_{frozen}|');

    title(ax,baselineNames{bRot});


    %% Binned median trend versus rho

    rhoMinThis = min(rhoThis);
    rhoMaxThis = max(rhoThis);

    if rhoMaxThis > rhoMinThis

        rhoEdges = linspace( ...
            rhoMinThis, ...
            rhoMaxThis, ...
            nRhoBins+1);

        rhoBin = discretize(rhoThis,rhoEdges);

        rhoBinCenter = NaN(nRhoBins,1);
        deltaBinMedian = NaN(nRhoBins,1);

        for jBin = 1:nRhoBins

            useBin = rhoBin == jBin;

            if any(useBin)

                rhoBinCenter(jBin) = ...
                    median(rhoThis(useBin));

                deltaBinMedian(jBin) = ...
                    median(deltaThis(useBin));

            end

        end

        validBins = ...
            isfinite(rhoBinCenter) & ...
            isfinite(deltaBinMedian);

        plot( ...
            ax, ...
            rhoBinCenter(validBins), ...
            deltaBinMedian(validBins), ...
            'k-', ...
            'LineWidth',2.0);

    end

end

colormap(figDeltaVvsRho,turbo(256));

cb = colorbar(ax);
cb.Label.String = 'Itokawa rotation phase [deg]';

title( ...
    tlDeltaVvsRho, ...
    'Visibility mismatch versus sampled spatial frequency');


%% Figure 2: mismatch versus asteroid phase, colored by spatial frequency

figDeltaVvsPhase = figure( ...
    'Name','Rotating-target mismatch versus asteroid phase', ...
    'Color','w', ...
    'Position',[80 100 1500 520]);

tlDeltaVvsPhase = tiledlayout( ...
    figDeltaVvsPhase, ...
    1,nBaselinesRot, ...
    'TileSpacing','compact', ...
    'Padding','compact');

rhoGlobalMin = min(rhoRot_Glambda(:));
rhoGlobalMax = max(rhoRot_Glambda(:));

for bRot = 1:nBaselinesRot

    ax = nexttile(tlDeltaVvsPhase);

    hold(ax,'on');
    grid(ax,'on');
    box(ax,'on');

    scatter( ...
        ax, ...
        rotationPhaseAccessible_deg, ...
        deltaVAbs(:,bRot), ...
        22, ...
        rhoRot_Glambda(:,bRot), ...
        'filled');

    clim(ax,[rhoGlobalMin rhoGlobalMax]);

    xlabel(ax,'Itokawa rotation phase [deg]');
    ylabel(ax,'|V_{rot}-V_{frozen}|');

    xlim(ax,[0 360]);

    title(ax,baselineNames{bRot});

end

colormap(figDeltaVvsPhase,turbo(256));

cb = colorbar(ax);
cb.Label.String = '\rho [G\lambda]';

title( ...
    tlDeltaVvsPhase, ...
    'Visibility mismatch versus asteroid phase');


%% Figure 3: separate amplitude and phase effects

figAmplitudePhaseMismatch = figure( ...
    'Name','Amplitude and phase effects of target rotation', ...
    'Color','w', ...
    'Position',[100 80 1300 800]);

tlAmplitudePhaseMismatch = tiledlayout( ...
    figAmplitudePhaseMismatch, ...
    2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% Amplitude-only difference

axAmpDifference = nexttile(tlAmplitudePhaseMismatch);

hold(axAmpDifference,'on');
grid(axAmpDifference,'on');
box(axAmpDifference,'on');

for bRot = 1:nBaselinesRot

    scatter( ...
        axAmpDifference, ...
        rhoRot_Glambda(:,bRot), ...
        deltaAmplitude(:,bRot), ...
        18, ...
        'filled', ...
        'DisplayName',baselineNames{bRot});

end

xlabel(axAmpDifference,'\rho [G\lambda]');
ylabel(axAmpDifference,'||V_{rot}|-|V_{frozen}||');

legend(axAmpDifference,'Location','best');

title(axAmpDifference, ...
    'Change in visibility amplitude caused by asteroid rotation');


%% Complex phase difference
%
% Only show points where both rotating and frozen visibility amplitudes
% are sufficiently large that the phase comparison is meaningful.

phaseComparisonThreshold = 0.05;

axPhaseDifference = nexttile(tlAmplitudePhaseMismatch);

hold(axPhaseDifference,'on');
grid(axPhaseDifference,'on');
box(axPhaseDifference,'on');

for bRot = 1:nBaselinesRot

    validPhaseDifference = ...
        abs(VRotatingQPO(:,bRot)) > phaseComparisonThreshold & ...
        abs(VFrozenFixedGrid(:,bRot)) > phaseComparisonThreshold;

    scatter( ...
        axPhaseDifference, ...
        rhoRot_Glambda(validPhaseDifference,bRot), ...
        phaseDifference_deg(validPhaseDifference,bRot), ...
        18, ...
        'filled', ...
        'DisplayName',baselineNames{bRot});

end

xlabel(axPhaseDifference,'\rho [G\lambda]');
ylabel(axPhaseDifference,'arg(V_{rot} V_{frozen}^*) [deg]');

ylim(axPhaseDifference,[-180 180]);

legend(axPhaseDifference,'Location','best');

title(axPhaseDifference, ...
    sprintf('Rotation-induced phase change for |V| > %.2f', ...
    phaseComparisonThreshold));

title( ...
    tlAmplitudePhaseMismatch, ...
    'How target rotation changes the complex visibility');
%% 18. CONTROLLED SINGLE-UV ROTATION EXPERIMENT
% ============================================================
%
% Hold one Fourier-plane sample fixed while Itokawa rotates.
%
% This removes the QPO evolution completely. The only changing variable is
% asteroid rotation phase.
%
% Fixed controlled sample:
%
%   rho = 0.20 Glambda
%   psi = 0 deg
%
% where
%
%   u = rho*cos(psi)
%   v = rho*sin(psi)


%% Controlled UV point

rhoControl_Glambda = 0.20;
psiControl_deg = 0.0;

rhoControl = rhoControl_Glambda*1e9;

uControl = rhoControl*cosd(psiControl_deg);
vControl = rhoControl*sind(psiControl_deg);

fprintf('\n');
fprintf('============================================================\n');
fprintf('CONTROLLED SINGLE-UV ROTATION EXPERIMENT\n');
fprintf('============================================================\n');
fprintf('  rho                = %.6f Glambda\n',rhoControl_Glambda);
fprintf('  psi                = %.3f deg\n',psiControl_deg);
fprintf('  u                  = %.6e wavelengths\n',uControl);
fprintf('  v                  = %.6e wavelengths\n',vControl);


%% Rotation phases
%
% Include both 0 and 360 deg so that one-period closure is visible.

phaseControl_deg = (0:5:360).';
nPhaseControl = numel(phaseControl_deg);

timeControl_s = ...
    rotItokawa.t0_s + ...
    (phaseControl_deg/360)*rotItokawa.period_s;


%% Allocate output

VControl = complex(zeros(nPhaseControl,1));

visibleFacetCountControl = zeros(nPhaseControl,1);
nonzeroPixelCountControl = zeros(nPhaseControl,1);
fluxErrorControl = zeros(nPhaseControl,1);

actualPhaseControl_deg = zeros(nPhaseControl,1);


%% Sweep asteroid rotation phase

fprintf('  phase samples       = %d\n',nPhaseControl);
fprintf('------------------------------------------------------------\n');

ticControl = tic;

for kPhase = 1:nPhaseControl

    tControl_s = timeControl_s(kPhase);

    [R_IB_control,phiControl,~] = ...
        asteroid_attitude_at_time(rotItokawa,tControl_s);

    actualPhaseControl_deg(kPhase) = rad2deg(phiControl);


    %% Instantaneous projected asteroid

    projControl = ...
        project_shape_orthographic( ...
            shapeItokawa, ...
            R_IB_control, ...
            target);

    visControl = ...
        find_visible_facets(projControl);


    %% Fixed-grid raster

    imgControl = ...
        rasterize_visible_facets( ...
            projControl, ...
            visControl, ...
            nPixItokawa, ...
            halfWidthRotatingItokawa_km, ...
            facetRadianceItokawa);

    simulatedControl = ...
        make_simulated_target_from_raster( ...
            imgControl, ...
            rangeItokawa_km, ...
            shapeItokawa.type);


    %% Extract nonzero pixels

    [jPixControl,iPixControl,fluxPixControl] = ...
        find(simulatedControl.I);

    if isempty(fluxPixControl)
        error('Controlled target has no nonzero pixels at phase %.3f deg.', ...
            phaseControl_deg(kPhase));
    end

    lPixControl_rad = ...
        simulatedControl.l_rad(iPixControl);

    mPixControl_rad = ...
        simulatedControl.m_rad(jPixControl);

    lPixControl_rad = lPixControl_rad(:);
    mPixControl_rad = mPixControl_rad(:);
    fluxPixControl = fluxPixControl(:);

    assert(iscolumn(lPixControl_rad));
    assert(iscolumn(mPixControl_rad));
    assert(iscolumn(fluxPixControl));

    assert(numel(lPixControl_rad) == numel(fluxPixControl));
    assert(numel(mPixControl_rad) == numel(fluxPixControl));


    %% Visibility at the single fixed UV point

    fourierPhaseControl = ...
        uControl*lPixControl_rad + ...
        vControl*mPixControl_rad;

    VControl(kPhase) = ...
        sum( ...
            fluxPixControl .* ...
            exp(-2*pi*1i*fourierPhaseControl));


    %% Diagnostics

    visibleFacetCountControl(kPhase) = ...
        nnz(visControl.isVisible);

    nonzeroPixelCountControl(kPhase) = ...
        numel(fluxPixControl);

    fluxErrorControl(kPhase) = ...
        abs(sum(fluxPixControl)-1);

end

runtimeControl_s = toc(ticControl);


%% Phase-zero reference

VControl0 = VControl(1);

deltaVControl = ...
    VControl-VControl0;

deltaVControlAbs = ...
    abs(deltaVControl);

deltaAmplitudeControl = ...
    abs(abs(VControl)-abs(VControl0));

deltaPhaseControl_deg = ...
    rad2deg( ...
        angle(VControl.*conj(VControl0)));


%% Closure check

visibilityClosureError = ...
    abs(VControl(end)-VControl(1));

attitudePhaseClosureError_deg = ...
    abs(actualPhaseControl_deg(end)-actualPhaseControl_deg(1));


%% Find maximum complex visibility departure

[maxDeltaVControl,indexMaxDeltaVControl] = ...
    max(deltaVControlAbs);

phaseMaxDeltaVControl_deg = ...
    phaseControl_deg(indexMaxDeltaVControl);


%% Summary

fprintf('\nCONTROLLED ROTATION RESULT\n');
fprintf('  runtime                      = %.3f s\n',runtimeControl_s);
fprintf('  |V| at phase 0 deg          = %.6f\n',abs(VControl0));
fprintf('  arg(V) at phase 0 deg       = %.6f deg\n', ...
    rad2deg(angle(VControl0)));
fprintf('  maximum |V(phi)-V(0)|       = %.6f\n',maxDeltaVControl);
fprintf('  phase of maximum mismatch   = %.3f deg\n', ...
    phaseMaxDeltaVControl_deg);
fprintf('  one-period visibility error = %.3e\n', ...
    visibilityClosureError);
fprintf('  max flux normalization err  = %.3e\n', ...
    max(fluxErrorControl));
fprintf('============================================================\n');

if abs(VControl0) < 0.05
    warning(['The phase-zero visibility amplitude is small. ', ...
        'Phase differences relative to V(0) should be interpreted carefully.']);
end


%% Figure 1: visibility amplitude and phase versus asteroid phase

figControlVisibility = figure( ...
    'Name','Controlled visibility versus Itokawa rotation phase', ...
    'Color','w', ...
    'Position',[100 100 1100 800]);

tlControlVisibility = tiledlayout( ...
    figControlVisibility, ...
    2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% Visibility amplitude

axControlAmp = nexttile(tlControlVisibility);

plot( ...
    axControlAmp, ...
    phaseControl_deg, ...
    abs(VControl), ...
    'o-', ...
    'LineWidth',1.4, ...
    'MarkerSize',4);

grid(axControlAmp,'on');
box(axControlAmp,'on');

xlim(axControlAmp,[0 360]);

xlabel(axControlAmp,'Itokawa rotation phase [deg]');
ylabel(axControlAmp,'|V|');

title(axControlAmp, ...
    sprintf('Fixed sample: rho = %.2f Glambda, psi = %.0f deg', ...
    rhoControl_Glambda,psiControl_deg));


%% Visibility phase

axControlPhase = nexttile(tlControlVisibility);

plot( ...
    axControlPhase, ...
    phaseControl_deg, ...
    rad2deg(angle(VControl)), ...
    'o-', ...
    'LineWidth',1.4, ...
    'MarkerSize',4);

grid(axControlPhase,'on');
box(axControlPhase,'on');

xlim(axControlPhase,[0 360]);
ylim(axControlPhase,[-180 180]);

xlabel(axControlPhase,'Itokawa rotation phase [deg]');
ylabel(axControlPhase,'arg(V) [deg]');

title(axControlPhase,'Complex visibility phase');


%% Figure 2: mismatch relative to phase zero

figControlMismatch = figure( ...
    'Name','Controlled visibility mismatch from phase zero', ...
    'Color','w', ...
    'Position',[120 120 1100 800]);

tlControlMismatch = tiledlayout( ...
    figControlMismatch, ...
    2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% Complex mismatch

axControlDeltaV = nexttile(tlControlMismatch);

plot( ...
    axControlDeltaV, ...
    phaseControl_deg, ...
    deltaVControlAbs, ...
    'o-', ...
    'LineWidth',1.4, ...
    'MarkerSize',4);

grid(axControlDeltaV,'on');
box(axControlDeltaV,'on');

xlim(axControlDeltaV,[0 360]);

xlabel(axControlDeltaV,'Itokawa rotation phase [deg]');
ylabel(axControlDeltaV,'|V(phi)-V(0)|');

title(axControlDeltaV,'Complex visibility departure from phase zero');


%% Wrapped phase difference

axControlDeltaPhase = nexttile(tlControlMismatch);

plot( ...
    axControlDeltaPhase, ...
    phaseControl_deg, ...
    deltaPhaseControl_deg, ...
    'o-', ...
    'LineWidth',1.4, ...
    'MarkerSize',4);

grid(axControlDeltaPhase,'on');
box(axControlDeltaPhase,'on');

xlim(axControlDeltaPhase,[0 360]);
ylim(axControlDeltaPhase,[-180 180]);

xlabel(axControlDeltaPhase,'Itokawa rotation phase [deg]');
ylabel(axControlDeltaPhase,'$\arg[V(\phi)V(0)^*]$ [deg]', ...
    'Interpreter','latex');

title(axControlDeltaPhase,'Wrapped phase change relative to phase zero');


%% Figure 3: trajectory of the visibility in the complex plane

figControlComplexPlane = figure( ...
    'Name','Controlled visibility trajectory in the complex plane', ...
    'Color','w', ...
    'Position',[150 150 750 700]);

axControlComplex = axes(figControlComplexPlane);

scatter( ...
    axControlComplex, ...
    real(VControl), ...
    imag(VControl), ...
    35, ...
    phaseControl_deg, ...
    'filled');

hold(axControlComplex,'on');

plot( ...
    axControlComplex, ...
    real(VControl), ...
    imag(VControl), ...
    '-', ...
    'LineWidth',1.0);

plot( ...
    axControlComplex, ...
    real(VControl0), ...
    imag(VControl0), ...
    'kp', ...
    'MarkerSize',13, ...
    'MarkerFaceColor','w');

grid(axControlComplex,'on');
box(axControlComplex,'on');
axis(axControlComplex,'equal');

xlabel(axControlComplex,'Re(V)');
ylabel(axControlComplex,'Im(V)');

cbControl = colorbar(axControlComplex);
cbControl.Label.String = 'Itokawa rotation phase [deg]';

title(axControlComplex, ...
    sprintf('Visibility trajectory at rho = %.2f Glambda, psi = %.0f deg', ...
    rhoControl_Glambda,psiControl_deg));
 


%% 19. LOCAL PHASE-BIN WIDTH AROUND THE REFERENCE ATTITUDE
% ============================================================
%
% Zoom in near phi0 = 0 deg for the same controlled UV sample used in
% Section 18. This estimates how wide a stroboscopic phase bin can be
% before the static-target approximation becomes inaccurate.
%
% Relative complex-visibility error:
%
%   epsilon_V(dphi) = |V(phi0+dphi)-V(phi0)| / |V(phi0)|
%
% We also track:
%
%   amplitude error  = ||V(phi)|-|V(phi0)|| / |V(phi0)|
%   phase error      = arg( V(phi) V(phi0)^* )
%
% The result is a local criterion for choosing a phase-bin half-width.


%% Use the same controlled UV sample as Section 18

if ~exist('rhoControl_Glambda','var') || ~exist('psiControl_deg','var')
    error(['Section 19 expects the controlled UV sample from Section 18 ', ...
        'to have been defined already.']);
end

rhoControl = rhoControl_Glambda*1e9;
uControl = rhoControl*cosd(psiControl_deg);
vControl = rhoControl*sind(psiControl_deg);


%% Reference phase and local phase offsets

phi0Local_deg = 0.0;

deltaPhaseLocal_deg = (-20:0.5:20).';
nDeltaLocal = numel(deltaPhaseLocal_deg);

phaseLocalWrapped_deg = mod(phi0Local_deg + deltaPhaseLocal_deg,360);
timeLocal_s = rotItokawa.t0_s + (phaseLocalWrapped_deg/360)*rotItokawa.period_s;


%% Allocate output

VLocal = complex(zeros(nDeltaLocal,1));

actualLocalPhase_deg = zeros(nDeltaLocal,1);
visibleFacetCountLocal = zeros(nDeltaLocal,1);
nonzeroPixelCountLocal = zeros(nDeltaLocal,1);
fluxErrorLocal = zeros(nDeltaLocal,1);


%% Evaluate the rotating target near the reference phase

fprintf('\n');
fprintf('============================================================\n');
fprintf('LOCAL STROBOSCOPIC BIN-WIDTH EXPERIMENT\n');
fprintf('============================================================\n');
fprintf('  controlled sample rho       = %.6f Glambda\n',rhoControl_Glambda);
fprintf('  controlled sample psi       = %.3f deg\n',psiControl_deg);
fprintf('  reference phase phi0        = %.3f deg\n',phi0Local_deg);
fprintf('  local phase range           = [%.3f, %.3f] deg\n', ...
    min(deltaPhaseLocal_deg),max(deltaPhaseLocal_deg));
fprintf('  local phase step            = %.3f deg\n', ...
    deltaPhaseLocal_deg(2)-deltaPhaseLocal_deg(1));

ticLocal = tic;

for kLocal = 1:nDeltaLocal

    tLocal_s = timeLocal_s(kLocal);

    [R_IB_local,phiLocal,~] = ...
        asteroid_attitude_at_time(rotItokawa,tLocal_s);

    actualLocalPhase_deg(kLocal) = rad2deg(phiLocal);


    %% Project and rasterize the asteroid at this local phase

    projLocal = ...
        project_shape_orthographic( ...
            shapeItokawa, ...
            R_IB_local, ...
            target);

    visLocal = find_visible_facets(projLocal);

    imgLocal = ...
        rasterize_visible_facets( ...
            projLocal, ...
            visLocal, ...
            nPixItokawa, ...
            halfWidthRotatingItokawa_km, ...
            facetRadianceItokawa);

    simulatedLocal = ...
        make_simulated_target_from_raster( ...
            imgLocal, ...
            rangeItokawa_km, ...
            shapeItokawa.type);


    %% Nonzero pixels

    [jPixLocal,iPixLocal,fluxPixLocal] = find(simulatedLocal.I);

    if isempty(fluxPixLocal)
        error('No nonzero pixels found at local phase offset %.3f deg.', ...
            deltaPhaseLocal_deg(kLocal));
    end

    lPixLocal_rad = simulatedLocal.l_rad(iPixLocal);
    mPixLocal_rad = simulatedLocal.m_rad(jPixLocal);

    lPixLocal_rad = lPixLocal_rad(:);
    mPixLocal_rad = mPixLocal_rad(:);
    fluxPixLocal = fluxPixLocal(:);

    assert(iscolumn(lPixLocal_rad));
    assert(iscolumn(mPixLocal_rad));
    assert(iscolumn(fluxPixLocal));

    assert(numel(lPixLocal_rad) == numel(fluxPixLocal));
    assert(numel(mPixLocal_rad) == numel(fluxPixLocal));


    %% Visibility at the fixed UV sample

    phaseFourierLocal = ...
        uControl*lPixLocal_rad + ...
        vControl*mPixLocal_rad;

    VLocal(kLocal) = ...
        sum(fluxPixLocal .* exp(-2*pi*1i*phaseFourierLocal));


    %% Diagnostics

    visibleFacetCountLocal(kLocal) = nnz(visLocal.isVisible);
    nonzeroPixelCountLocal(kLocal) = numel(fluxPixLocal);
    fluxErrorLocal(kLocal) = abs(sum(fluxPixLocal)-1);

end

runtimeLocal_s = toc(ticLocal);


%% Reference visibility at delta phase = 0

indexPhi0Local = find(abs(deltaPhaseLocal_deg) < 1e-12,1,'first');

if isempty(indexPhi0Local)
    error('Could not find the local reference sample at delta phase = 0.');
end

V0Local = VLocal(indexPhi0Local);

if abs(V0Local) < 1e-14
    error('Reference visibility magnitude is too small for relative-error analysis.');
end


%% Relative error metrics

epsilonComplexLocal = abs(VLocal-V0Local)/abs(V0Local);

epsilonAmplitudeLocal = ...
    abs(abs(VLocal)-abs(V0Local))/abs(V0Local);

deltaPhaseWrappedLocal_deg = ...
    rad2deg(angle(VLocal.*conj(V0Local)));


%% Symmetric half-widths for selected tolerances
%
% For each tolerance tau, we find the largest half-width H such that
% epsilonComplexLocal <= tau for every sampled point with |dphi| <= H.

toleranceList = [0.01; 0.05; 0.10];
binHalfWidth_deg = NaN(size(toleranceList));

absDeltaLocal_deg = abs(deltaPhaseLocal_deg);
uniqueAbsDelta_deg = unique(absDeltaLocal_deg);

for iTol = 1:numel(toleranceList)

    tau = toleranceList(iTol);
    validHalfWidth = NaN;

    for kWidth = 1:numel(uniqueAbsDelta_deg)

        H = uniqueAbsDelta_deg(kWidth);
        useWindow = absDeltaLocal_deg <= H + 1e-12;

        if all(epsilonComplexLocal(useWindow) <= tau + 1e-14)
            validHalfWidth = H;
        else
            break;
        end

    end

    binHalfWidth_deg(iTol) = validHalfWidth;

end


%% Report summary

maxComplexErrorLocal = max(epsilonComplexLocal);
maxAmplitudeErrorLocal = max(epsilonAmplitudeLocal);
maxWrappedPhaseErrorLocal = max(abs(deltaPhaseWrappedLocal_deg));

fprintf('\nLOCAL BIN-WIDTH RESULT\n');
fprintf('  runtime                      = %.3f s\n',runtimeLocal_s);
fprintf('  |V(phi0)|                    = %.6f\n',abs(V0Local));
fprintf('  arg(V(phi0))                 = %.6f deg\n',rad2deg(angle(V0Local)));
fprintf('  max relative complex error   = %.6f\n',maxComplexErrorLocal);
fprintf('  max relative amplitude error = %.6f\n',maxAmplitudeErrorLocal);
fprintf('  max wrapped phase change     = %.6f deg\n',maxWrappedPhaseErrorLocal);
fprintf('  max flux normalization err   = %.3e\n',max(fluxErrorLocal));
fprintf('\n');

for iTol = 1:numel(toleranceList)

    if isnan(binHalfWidth_deg(iTol))
        fprintf('  half-width for %.1f%% error   = below sampled resolution\n', ...
            100*toleranceList(iTol));
    else
        fprintf('  half-width for %.1f%% error   = +/- %.3f deg\n', ...
            100*toleranceList(iTol),binHalfWidth_deg(iTol));
    end

end

fprintf('============================================================\n');


%% Figure 1: local error metrics versus phase offset

figLocalErrors = figure( ...
    'Name','Local stroboscopic bin-width metrics', ...
    'Color','w', ...
    'Position',[100 100 1100 850]);

tlLocalErrors = tiledlayout( ...
    figLocalErrors, ...
    3,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% Relative complex-visibility error

axLocalComplex = nexttile(tlLocalErrors);

plot( ...
    axLocalComplex, ...
    deltaPhaseLocal_deg, ...
    epsilonComplexLocal, ...
    'o-', ...
    'LineWidth',1.4, ...
    'MarkerSize',4);

hold(axLocalComplex,'on');
grid(axLocalComplex,'on');
box(axLocalComplex,'on');

for iTol = 1:numel(toleranceList)
    yline(axLocalComplex,toleranceList(iTol),'--');
end

xlabel(axLocalComplex,'Phase offset from \phi_0 [deg]');
ylabel(axLocalComplex,'|V(\phi)-V(\phi_0)| / |V(\phi_0)|');

title(axLocalComplex, ...
    sprintf('Relative complex-visibility error near \\phi_0 = %.0f deg', ...
    phi0Local_deg));


%% Relative amplitude error

axLocalAmplitude = nexttile(tlLocalErrors);

plot( ...
    axLocalAmplitude, ...
    deltaPhaseLocal_deg, ...
    epsilonAmplitudeLocal, ...
    'o-', ...
    'LineWidth',1.4, ...
    'MarkerSize',4);

grid(axLocalAmplitude,'on');
box(axLocalAmplitude,'on');

xlabel(axLocalAmplitude,'Phase offset from \phi_0 [deg]');
ylabel(axLocalAmplitude,'||V(\phi)|-|V(\phi_0)|| / |V(\phi_0)|');

title(axLocalAmplitude,'Relative visibility-amplitude error');


%% Wrapped phase change

axLocalPhase = nexttile(tlLocalErrors);

plot( ...
    axLocalPhase, ...
    deltaPhaseLocal_deg, ...
    deltaPhaseWrappedLocal_deg, ...
    'o-', ...
    'LineWidth',1.4, ...
    'MarkerSize',4);

grid(axLocalPhase,'on');
box(axLocalPhase,'on');

xlabel(axLocalPhase,'Phase offset from \phi_0 [deg]');
ylabel(axLocalPhase,'$\arg[V(\phi)V(\phi_0)^*]$ [deg]', ...
    'Interpreter','latex');

title(axLocalPhase,'Wrapped visibility-phase change');


%% Figure 2: local complex-plane zoom

figLocalComplexPlane = figure( ...
    'Name','Local complex-plane trajectory near the reference phase', ...
    'Color','w', ...
    'Position',[150 150 800 700]);

axLocalComplexPlane = axes(figLocalComplexPlane);

scatter( ...
    axLocalComplexPlane, ...
    real(VLocal), ...
    imag(VLocal), ...
    45, ...
    deltaPhaseLocal_deg, ...
    'filled');

hold(axLocalComplexPlane,'on');

plot( ...
    axLocalComplexPlane, ...
    real(VLocal), ...
    imag(VLocal), ...
    '-', ...
    'LineWidth',1.0);

plot( ...
    axLocalComplexPlane, ...
    real(V0Local), ...
    imag(V0Local), ...
    'kp', ...
    'MarkerSize',13, ...
    'MarkerFaceColor','w');

grid(axLocalComplexPlane,'on');
box(axLocalComplexPlane,'on');
axis(axLocalComplexPlane,'equal');

xlabel(axLocalComplexPlane,'Re(V)');
ylabel(axLocalComplexPlane,'Im(V)');

cbLocal = colorbar(axLocalComplexPlane);
cbLocal.Label.String = 'Phase offset from \phi_0 [deg]';

title(axLocalComplexPlane, ...
    sprintf(['Local complex-plane trajectory near \\phi_0 = %.0f deg ', ...
    '(rho = %.2f Glambda, psi = %.0f deg)'], ...
    phi0Local_deg,rhoControl_Glambda,psiControl_deg));

%% 20. STROBOSCOPIC BIN WIDTH VERSUS SPATIAL FREQUENCY
% ============================================================
%
% Controlled experiment:
%
%   psi = 0 deg
%
% while varying
%
%   rho = 0.10, 0.20, 0.30, 0.40 Glambda.
%
% At every asteroid phase the target image is generated only once.
% The four visibility samples are then evaluated from that same image.
%
% We compute the largest symmetric phase-bin half-width H such that
%
%   |V(phi)-V(phi0)| / |V(phi0)| <= tolerance
%
% for every sampled phase satisfying
%
%   |phi-phi0| <= H.


%% Controlled Fourier samples

rhoMulti_Glambda = [0.10 0.20 0.30 0.40];
psiMulti_deg = 0.0;

nRhoMulti = numel(rhoMulti_Glambda);

rhoMulti = rhoMulti_Glambda*1e9;

uMulti = rhoMulti*cosd(psiMulti_deg);
vMulti = rhoMulti*sind(psiMulti_deg);


%% Local phase grid
%
% 0.25 deg resolution is four times finer than Section 19.

phi0Multi_deg = 0.0;

deltaPhaseMulti_deg = (-20:0.25:20).';
nPhaseMulti = numel(deltaPhaseMulti_deg);

phaseMultiWrapped_deg = ...
    mod(phi0Multi_deg + deltaPhaseMulti_deg,360);

timeMulti_s = ...
    rotItokawa.t0_s + ...
    (phaseMultiWrapped_deg/360)*rotItokawa.period_s;


%% Allocate visibility matrix
%
% Rows    = asteroid phase offsets
% Columns = controlled rho values

VMulti = complex(zeros(nPhaseMulti,nRhoMulti));

fluxErrorMulti = zeros(nPhaseMulti,1);


%% Run controlled asteroid sweep

fprintf('\n');
fprintf('============================================================\n');
fprintf('STROBOSCOPIC BIN WIDTH VERSUS SPATIAL FREQUENCY\n');
fprintf('============================================================\n');
fprintf('  psi                         = %.3f deg\n',psiMulti_deg);
fprintf('  rho values                  = ');

for iRho = 1:nRhoMulti
    fprintf('%.2f ',rhoMulti_Glambda(iRho));
end

fprintf('Glambda\n');
fprintf('  reference phase             = %.3f deg\n',phi0Multi_deg);
fprintf('  local phase range           = [%.3f, %.3f] deg\n', ...
    min(deltaPhaseMulti_deg),max(deltaPhaseMulti_deg));
fprintf('  phase resolution            = %.3f deg\n', ...
    deltaPhaseMulti_deg(2)-deltaPhaseMulti_deg(1));
fprintf('  target rasters              = %d\n',nPhaseMulti);
fprintf('------------------------------------------------------------\n');

ticMulti = tic;

for kPhase = 1:nPhaseMulti

    tMulti_s = timeMulti_s(kPhase);

    [R_IB_multi,~,~] = ...
        asteroid_attitude_at_time(rotItokawa,tMulti_s);


    %% Instantaneous asteroid image

    projMulti = ...
        project_shape_orthographic( ...
            shapeItokawa, ...
            R_IB_multi, ...
            target);

    visMulti = find_visible_facets(projMulti);

    imgMulti = ...
        rasterize_visible_facets( ...
            projMulti, ...
            visMulti, ...
            nPixItokawa, ...
            halfWidthRotatingItokawa_km, ...
            facetRadianceItokawa);

    simulatedMulti = ...
        make_simulated_target_from_raster( ...
            imgMulti, ...
            rangeItokawa_km, ...
            shapeItokawa.type);


    %% Nonzero image pixels

    [jPixMulti,iPixMulti,fluxPixMulti] = find(simulatedMulti.I);

    if isempty(fluxPixMulti)
        error('No nonzero pixels at phase offset %.3f deg.', ...
            deltaPhaseMulti_deg(kPhase));
    end

    lPixMulti_rad = simulatedMulti.l_rad(iPixMulti);
    mPixMulti_rad = simulatedMulti.m_rad(jPixMulti);

    lPixMulti_rad = lPixMulti_rad(:);
    mPixMulti_rad = mPixMulti_rad(:);
    fluxPixMulti = fluxPixMulti(:);

    assert(iscolumn(lPixMulti_rad));
    assert(iscolumn(mPixMulti_rad));
    assert(iscolumn(fluxPixMulti));

    assert(numel(lPixMulti_rad) == numel(fluxPixMulti));
    assert(numel(mPixMulti_rad) == numel(fluxPixMulti));

    fluxErrorMulti(kPhase) = abs(sum(fluxPixMulti)-1);


    %% Evaluate all controlled rho values from the same image

    for iRho = 1:nRhoMulti

        fourierPhaseMulti = ...
            uMulti(iRho)*lPixMulti_rad + ...
            vMulti(iRho)*mPixMulti_rad;

        VMulti(kPhase,iRho) = ...
            sum(fluxPixMulti .* ...
            exp(-2*pi*1i*fourierPhaseMulti));

    end


    %% Progress

    if kPhase == 1 || mod(kPhase,40) == 0 || kPhase == nPhaseMulti

        fprintf('  phase sample %4d / %4d: delta phi = %+7.2f deg\n', ...
            kPhase,nPhaseMulti,deltaPhaseMulti_deg(kPhase));

    end

end

runtimeMulti_s = toc(ticMulti);


%% Reference visibility at delta phi = 0

indexPhi0Multi = ...
    find(abs(deltaPhaseMulti_deg) < 1e-12,1,'first');

if isempty(indexPhi0Multi)
    error('Could not locate delta phi = 0 in the multi-rho phase grid.');
end

V0Multi = VMulti(indexPhi0Multi,:);

if any(abs(V0Multi) < 1e-14)
    error('At least one reference visibility is too close to zero.');
end


%% Relative complex-visibility error

epsilonMulti = ...
    abs(VMulti - V0Multi) ./ abs(V0Multi);


%% Relative amplitude error

epsilonAmplitudeMulti = ...
    abs(abs(VMulti) - abs(V0Multi)) ./ abs(V0Multi);


%% Wrapped phase change

deltaVisibilityPhaseMulti_deg = ...
    rad2deg(angle(VMulti .* conj(V0Multi)));


%% Determine symmetric phase-bin half-widths

toleranceMulti = [0.01 0.05 0.10];
nToleranceMulti = numel(toleranceMulti);

halfWidthMulti_deg = ...
    NaN(nToleranceMulti,nRhoMulti);

absDeltaMulti_deg = abs(deltaPhaseMulti_deg);
uniqueAbsDeltaMulti_deg = unique(absDeltaMulti_deg);

for iRho = 1:nRhoMulti

    for iTol = 1:nToleranceMulti

        tau = toleranceMulti(iTol);
        validHalfWidth = 0.0;

        for iWidth = 1:numel(uniqueAbsDeltaMulti_deg)

            H = uniqueAbsDeltaMulti_deg(iWidth);

            useWindow = ...
                absDeltaMulti_deg <= H + 1e-12;

            if all(epsilonMulti(useWindow,iRho) <= tau + 1e-14)

                validHalfWidth = H;

            else

                break;

            end

        end

        halfWidthMulti_deg(iTol,iRho) = validHalfWidth;

    end

end


%% Convert phase half-width to observing-time half-width

minutesPerDegreeItokawa = ...
    rotItokawa.period_s/(360*60);

halfWidthMulti_min = ...
    halfWidthMulti_deg*minutesPerDegreeItokawa;

fullWindowMulti_min = ...
    2*halfWidthMulti_min;


%% Numerical summary

smallestPositiveOffsetMulti_deg = ...
    min(uniqueAbsDeltaMulti_deg(uniqueAbsDeltaMulti_deg > 0));

fprintf('\n');
fprintf('MULTI-RHO BIN-WIDTH RESULT\n');
fprintf('  runtime                      = %.3f s\n',runtimeMulti_s);
fprintf('  phase-grid resolution        = %.3f deg\n', ...
    smallestPositiveOffsetMulti_deg);
fprintf('  time per asteroid degree     = %.6f min/deg\n', ...
    minutesPerDegreeItokawa);
fprintf('  max flux normalization err   = %.3e\n', ...
    max(fluxErrorMulti));
fprintf('\n');


for iRho = 1:nRhoMulti

    fprintf('rho = %.2f Glambda\n',rhoMulti_Glambda(iRho));
    fprintf('  |V(phi0)|                  = %.6f\n', ...
        abs(V0Multi(iRho)));
    fprintf('  arg(V(phi0))               = %+.6f deg\n', ...
        rad2deg(angle(V0Multi(iRho))));

    for iTol = 1:nToleranceMulti

        Hdeg = halfWidthMulti_deg(iTol,iRho);
        Hmin = halfWidthMulti_min(iTol,iRho);

        if Hdeg == 0

            fprintf( ...
                '  %.1f%% half-width             < %.3f deg (unresolved)\n', ...
                100*toleranceMulti(iTol), ...
                smallestPositiveOffsetMulti_deg);

        else

            fprintf( ...
                '  %.1f%% half-width             = +/- %.3f deg', ...
                100*toleranceMulti(iTol), ...
                Hdeg);

            fprintf( ...
                ' = +/- %.3f min\n', ...
                Hmin);

        end

    end

    fprintf('\n');

end

fprintf('============================================================\n');


%% Store controlled-study results

stroboscopicRhoStudy = struct();

stroboscopicRhoStudy.rho_Glambda = rhoMulti_Glambda;
stroboscopicRhoStudy.psi_deg = psiMulti_deg;

stroboscopicRhoStudy.deltaPhase_deg = deltaPhaseMulti_deg;

stroboscopicRhoStudy.V = VMulti;
stroboscopicRhoStudy.V0 = V0Multi;

stroboscopicRhoStudy.relativeComplexError = epsilonMulti;
stroboscopicRhoStudy.relativeAmplitudeError = epsilonAmplitudeMulti;

stroboscopicRhoStudy.visibilityPhaseChange_deg = ...
    deltaVisibilityPhaseMulti_deg;

stroboscopicRhoStudy.tolerance = toleranceMulti;

stroboscopicRhoStudy.halfWidth_deg = halfWidthMulti_deg;
stroboscopicRhoStudy.halfWidth_min = halfWidthMulti_min;
stroboscopicRhoStudy.fullWindow_min = fullWindowMulti_min;


%% Figure 1: relative complex error for each rho

figRhoError = figure( ...
    'Name','Stroboscopic error versus spatial frequency', ...
    'Color','w', ...
    'Position',[100 100 1100 750]);

axRhoError = axes(figRhoError);

hold(axRhoError,'on');
grid(axRhoError,'on');
box(axRhoError,'on');


legendRho = cell(1,nRhoMulti);

for iRho = 1:nRhoMulti

    plot( ...
        axRhoError, ...
        deltaPhaseMulti_deg, ...
        epsilonMulti(:,iRho), ...
        'LineWidth',1.6);

    legendRho{iRho} = [ ...
        '$\rho=', ...
        num2str(rhoMulti_Glambda(iRho),'%.2f'), ...
        '\,\mathrm{G}\lambda$'];

end


for iTol = 1:nToleranceMulti

    yline( ...
        axRhoError, ...
        toleranceMulti(iTol), ...
        '--');

end


xlabel( ...
    axRhoError, ...
    '$\Delta\phi$ [deg]', ...
    'Interpreter','latex');

ylabel( ...
    axRhoError, ...
    '$|V(\phi)-V(\phi_0)|/|V(\phi_0)|$', ...
    'Interpreter','latex');

title( ...
    axRhoError, ...
    'Local complex-visibility error at fixed baseline direction', ...
    'Interpreter','none');

legend( ...
    axRhoError, ...
    legendRho, ...
    'Interpreter','latex', ...
    'Location','best');


%% Figure 2: allowable phase-bin half-width versus rho

figRhoHalfWidth = figure( ...
    'Name','Allowable stroboscopic phase-bin width versus spatial frequency', ...
    'Color','w', ...
    'Position',[120 120 1000 700]);

axRhoHalfWidth = axes(figRhoHalfWidth);

hold(axRhoHalfWidth,'on');
grid(axRhoHalfWidth,'on');
box(axRhoHalfWidth,'on');


legendTolerance = cell(1,nToleranceMulti);

for iTol = 1:nToleranceMulti

    plot( ...
        axRhoHalfWidth, ...
        rhoMulti_Glambda, ...
        halfWidthMulti_deg(iTol,:), ...
        'o-', ...
        'LineWidth',1.6, ...
        'MarkerSize',7);

    legendTolerance{iTol} = [ ...
        num2str(100*toleranceMulti(iTol),'%.0f'), ...
        '\% visibility-error tolerance'];

end


xlabel( ...
    axRhoHalfWidth, ...
    '$\rho$ [G$\lambda$]', ...
    'Interpreter','latex');

ylabel( ...
    axRhoHalfWidth, ...
    'Symmetric phase-bin half-width [deg]', ...
    'Interpreter','none');

title( ...
    axRhoHalfWidth, ...
    'Allowable phase-bin width versus spatial frequency', ...
    'Interpreter','none');

legend( ...
    axRhoHalfWidth, ...
    legendTolerance, ...
    'Interpreter','latex', ...
    'Location','best');


%% Figure 3: same half-width expressed as observing time

figRhoTimeWindow = figure( ...
    'Name','Stroboscopic observing-time window versus spatial frequency', ...
    'Color','w', ...
    'Position',[140 140 1000 700]);

axRhoTimeWindow = axes(figRhoTimeWindow);

hold(axRhoTimeWindow,'on');
grid(axRhoTimeWindow,'on');
box(axRhoTimeWindow,'on');


for iTol = 1:nToleranceMulti

    plot( ...
        axRhoTimeWindow, ...
        rhoMulti_Glambda, ...
        halfWidthMulti_min(iTol,:), ...
        'o-', ...
        'LineWidth',1.6, ...
        'MarkerSize',7);

end


xlabel( ...
    axRhoTimeWindow, ...
    '$\rho$ [G$\lambda$]', ...
    'Interpreter','latex');

ylabel( ...
    axRhoTimeWindow, ...
    'Observing-window half-width [min]', ...
    'Interpreter','none');

title( ...
    axRhoTimeWindow, ...
    'Allowed observing time around each asteroid-phase recurrence', ...
    'Interpreter','none');

legend( ...
    axRhoTimeWindow, ...
    legendTolerance, ...
    'Interpreter','latex', ...
    'Location','best');

%% ============================================================
%  21. PHASE-SYNCHRONOUS NATURAL-QPO OBSERVATION GRID
%  ============================================================
%
% Purpose
% -------
% Build ONE master data set from which every stroboscopic binning case is
% drawn.
%
% The formation is never retimed or controlled here.  The three collectors
% follow the same natural QPO used above.  Observation times are prescribed
% by asteroid rotational phase:
%
%   t(n,DeltaPhi) = t0 + n*Tasteroid
%                 + (DeltaPhi/360 deg)*Tasteroid/spinSense.
%
% DeltaPhi is represented on [-180,180).  A finite FULL phase-bin width W
% later accepts samples satisfying
%
%   |DeltaPhi| <= W/2.
%
% Thus W = 0 deg is the exact-phase stroboscopic limit and W = 360 deg is
% the unbinned full-rotation endpoint on the same controlled cadence.
%
% IMPORTANT:
%   phaseSampleStep_deg controls the temporal sampling inside each phase
%   bin.  The final threshold should be refined with a smaller value after
%   the first broad sweep has located the transition region.

phaseSampleStep_deg = 5.0;

phaseBinWidth_deg = [ ...
      0 ...
     10 ...
     20 ...
     30 ...
     40 ...
     60 ...
     90 ...
    120 ...
    180 ...
    240 ...
    360];

campaignFractions = [ ...
    0.25 ...
    0.50 ...
    1.00 ...
    2.00 ...
    3.00];

% Minimum number of distinct accessible observing epochs required before
% attempting a CLEAN reconstruction.
minAccessibleEpochsPerCase = 5;

% Block size used only for the direct rotating-target NUDFT below.
% It bounds memory while evaluating several visibilities at once.
visibilityBlockSize = 16;

% Reconstruction-fidelity criteria used only in Section 25.
% These are STUDY CRITERIA, not fundamental information-theory limits.
fidelitySupportNRMSEMax = 0.20;
fidelityNCCMin = 0.90;
maxBlurNRMSEPenalty = 0.05;
maxBlurNCCPenalty = 0.05;

% Numerical tolerances for this study.
phaseGridTolerance_deg = 1e-10;
phaseClosureTolerance_deg = 1e-8;
zeroPhaseVisibilityTolerance = 1e-10;
beamEqualityTolerance = 1e-11;


%% Validate the phase grid

if ~isscalar(phaseSampleStep_deg) || ...
        ~isfinite(phaseSampleStep_deg) || ...
        phaseSampleStep_deg <= 0 || ...
        phaseSampleStep_deg > 180

    error('phaseSampleStep_deg must lie in (0,180] deg.');

end

nPhaseOffsets = round(360/phaseSampleStep_deg);

if abs(nPhaseOffsets*phaseSampleStep_deg - 360) > 1e-10
    error('phaseSampleStep_deg must divide 360 deg exactly.');
end

phaseOffsetGrid_deg = ...
    -180 + (0:nPhaseOffsets-1).' * phaseSampleStep_deg;

[minimumZeroOffset_deg,indexPhaseZero] = ...
    min(abs(phaseOffsetGrid_deg));

if minimumZeroOffset_deg > phaseGridTolerance_deg
    error('The phase grid must contain DeltaPhi = 0 exactly.');
end

phaseOffsetGrid_deg(indexPhaseZero) = 0;

phaseBinWidth_deg = unique(sort(phaseBinWidth_deg(:).'));

if isempty(phaseBinWidth_deg) || ...
        phaseBinWidth_deg(1) < 0 || ...
        phaseBinWidth_deg(end) > 360 || ...
        ~any(abs(phaseBinWidth_deg) <= phaseGridTolerance_deg)

    error(['phaseBinWidth_deg must lie in [0,360] deg and ', ...
        'must include the exact-phase W = 0 case.']);

end

if ~ismember(rotItokawa.spinSense,[-1 1])
    error('rotItokawa.spinSense must be +1 or -1.');
end


%% Time scales and campaign durations

tQPO_s = qpoElapsedRot_s(:);

if any(diff(tQPO_s) <= 0)
    error('qpoElapsedRot_s must be strictly increasing.');
end

tQPOStart_s = tQPO_s(1);
tQPOEnd_s = tQPO_s(end);
qpoDuration_s = tQPOEnd_s-tQPOStart_s;

Tstar_s = sys.Tstar_days*86400;
TPrimary_s = TPrimary*Tstar_s;
TItokawa_s = rotItokawa.period_s;

campaignFractions = unique(sort(campaignFractions(:).'));
campaignDuration_s = campaignFractions*TPrimary_s;

keepCampaign = ...
    campaignDuration_s <= qpoDuration_s + ...
    100*eps(max(1,qpoDuration_s));

campaignFractions = campaignFractions(keepCampaign);
campaignDuration_s = campaignDuration_s(keepCampaign);

if isempty(campaignFractions)
    error('No requested campaign duration fits inside the propagated QPO interval.');
end

campaignDuration_days = campaignDuration_s/86400;

nCampaignCases = numel(campaignFractions);
nBinCases = numel(phaseBinWidth_deg);


%% Build the phase-synchronous master observation schedule
%
% Include one extra recurrence center on each side, then retain only times
% that actually fall inside the propagated QPO interval.

nFirstCenter = floor( ...
    (tQPOStart_s-rotItokawa.t0_s)/TItokawa_s) - 1;

nLastCenter = ceil( ...
    (tQPOEnd_s-rotItokawa.t0_s)/TItokawa_s) + 1;

recurrenceCenterNumber = (nFirstCenter:nLastCenter).';
nRecurrenceCenters = numel(recurrenceCenterNumber);

nScheduleMax = nRecurrenceCenters*nPhaseOffsets;

masterTime_s = zeros(nScheduleMax,1);
masterPhaseOffset_deg = zeros(nScheduleMax,1);
masterPhaseIndex = zeros(nScheduleMax,1);
masterRecurrenceNumber = zeros(nScheduleMax,1);

scheduleCount = 0;

for iCenter = 1:nRecurrenceCenters

    recurrenceNumber = recurrenceCenterNumber(iCenter);

    recurrenceCenterTime_s = ...
        rotItokawa.t0_s + recurrenceNumber*TItokawa_s;

    for iPhase = 1:nPhaseOffsets

        phaseOffset_deg = phaseOffsetGrid_deg(iPhase);

        timeOffset_s = ...
            (phaseOffset_deg/360) * ...
            TItokawa_s / ...
            rotItokawa.spinSense;

        candidateTime_s = ...
            recurrenceCenterTime_s + timeOffset_s;

        insideQPO = ...
            candidateTime_s >= tQPOStart_s - 1e-9 && ...
            candidateTime_s <= tQPOEnd_s + 1e-9;

        if ~insideQPO
            continue
        end

        scheduleCount = scheduleCount + 1;

        masterTime_s(scheduleCount) = candidateTime_s;
        masterPhaseOffset_deg(scheduleCount) = phaseOffset_deg;
        masterPhaseIndex(scheduleCount) = iPhase;
        masterRecurrenceNumber(scheduleCount) = recurrenceNumber;

    end

end

masterTime_s = masterTime_s(1:scheduleCount);
masterPhaseOffset_deg = masterPhaseOffset_deg(1:scheduleCount);
masterPhaseIndex = masterPhaseIndex(1:scheduleCount);
masterRecurrenceNumber = masterRecurrenceNumber(1:scheduleCount);

[masterTime_s,sortMaster] = sort(masterTime_s);

masterPhaseOffset_deg = masterPhaseOffset_deg(sortMaster);
masterPhaseIndex = masterPhaseIndex(sortMaster);
masterRecurrenceNumber = masterRecurrenceNumber(sortMaster);

nMasterEpochs = numel(masterTime_s);

if nMasterEpochs < 2
    error('Phase-synchronous master schedule contains fewer than two epochs.');
end

minimumMasterTimeSpacing_s = min(diff(masterTime_s));

if minimumMasterTimeSpacing_s <= 0
    error('Duplicate or nonmonotonic times were generated in the master schedule.');
end


%% Verify that each scheduled time has the requested asteroid phase

actualMasterPhase_rad = zeros(nMasterEpochs,1);

for kMaster = 1:nMasterEpochs

    [~,actualMasterPhase_rad(kMaster),~] = ...
        asteroid_attitude_at_time( ...
            rotItokawa, ...
            masterTime_s(kMaster));

end

requestedMasterPhase_rad = ...
    deg2rad(masterPhaseOffset_deg);

phaseScheduleError_rad = ...
    atan2( ...
        sin(actualMasterPhase_rad-requestedMasterPhase_rad), ...
        cos(actualMasterPhase_rad-requestedMasterPhase_rad));

maxPhaseScheduleError_deg = ...
    max(abs(rad2deg(phaseScheduleError_rad)));

if maxPhaseScheduleError_deg > phaseClosureTolerance_deg
    error(['Master schedule does not reproduce the requested asteroid ', ...
        'phases. Max phase error = %.3e deg.'], ...
        maxPhaseScheduleError_deg);
end


%% Dense natural-QPO collector solutions
%
% These reproduce the already propagated collector histories on geomQPO.t,
% but also permit exact evaluation at the new observation times.

solCollectorsStrobe = cell(1,3);
maxDenseGridDifference = zeros(3,1);

for iCollector = 1:3

    [YDenseGrid,solCollectorsStrobe{iCollector}] = ...
        propagate_state_on_grid( ...
            XCollectors0(:,iCollector), ...
            geomQPO.t, ...
            mu, ...
            odeOpts);

    existingGridHistory = ...
        geomQPO.YCollectors(:,:,iCollector);

    if ~isequal(size(YDenseGrid),size(existingGridHistory))
        error('Dense collector history has an unexpected size.');
    end

    maxDenseGridDifference(iCollector) = ...
        max(abs(YDenseGrid-existingGridHistory),[],'all');

end


%% Evaluate the natural QPO at every master observation epoch

tMaster_nd = ...
    geomQPO.t(1) + masterTime_s/Tstar_s;

if any(tMaster_nd < geomQPO.t(1)-1e-13) || ...
        any(tMaster_nd > geomQPO.t(end)+1e-13)

    error('At least one master observation epoch lies outside the dense QPO solution.');

end

Y1Master = deval(solCollectorsStrobe{1},tMaster_nd).';
Y2Master = deval(solCollectorsStrobe{2},tMaster_nd).';
Y3Master = deval(solCollectorsStrobe{3},tMaster_nd).';

r1Master = Y1Master(:,1:3);
r2Master = Y2Master(:,1:3);
r3Master = Y3Master(:,1:3);


%% Physical baselines in the synodic frame

b12SynodicMaster = r2Master-r1Master;
b13SynodicMaster = r3Master-r1Master;
b23SynodicMaster = r3Master-r2Master;

baselineClosureSynodicMaster = ...
    b12SynodicMaster + ...
    b23SynodicMaster - ...
    b13SynodicMaster;

maxBaselineClosureSynodicMaster = ...
    max(vecnorm(baselineClosureSynodicMaster,2,2));


%% Rotate to inertial and project into the fixed target basis

b12InertialMaster = ...
    rotate_synodic_to_inertial_rows( ...
        b12SynodicMaster, ...
        tMaster_nd);

b13InertialMaster = ...
    rotate_synodic_to_inertial_rows( ...
        b13SynodicMaster, ...
        tMaster_nd);

b23InertialMaster = ...
    rotate_synodic_to_inertial_rows( ...
        b23SynodicMaster, ...
        tMaster_nd);

eUMaster = target.eUInertial0(:);
eVMaster = target.eVInertial0(:);
eWMaster = target.eWInertial0(:);

lengthScaleWavelengths = ...
    sys.Lstar_km*1e3/lambda_m;

uvw12Master = ...
    lengthScaleWavelengths * [ ...
        b12InertialMaster*eUMaster, ...
        b12InertialMaster*eVMaster, ...
        b12InertialMaster*eWMaster];

uvw13Master = ...
    lengthScaleWavelengths * [ ...
        b13InertialMaster*eUMaster, ...
        b13InertialMaster*eVMaster, ...
        b13InertialMaster*eWMaster];

uvw23Master = ...
    lengthScaleWavelengths * [ ...
        b23InertialMaster*eUMaster, ...
        b23InertialMaster*eVMaster, ...
        b23InertialMaster*eWMaster];

uvwClosureMaster = ...
    uvw12Master + ...
    uvw23Master - ...
    uvw13Master;

maxUVWClosureMaster = ...
    max(vecnorm(uvwClosureMaster,2,2));


%% Accessibility at every exact master epoch

triangleNormalMaster = ...
    cross(uvw12Master,uvw13Master,2);

triangleNormalNormMaster = ...
    vecnorm(triangleNormalMaster,2,2);

if any(triangleNormalNormMaster <= 0)
    error('Degenerate collector triangle encountered in the master schedule.');
end

cosOffNormalMaster = ...
    triangleNormalMaster(:,3) ./ ...
    triangleNormalNormMaster;

if usePlaneNormalSignAmbiguity
    cosOffNormalMaster = abs(cosOffNormalMaster);
end

cosOffNormalMaster = ...
    min(1,max(-1,cosOffNormalMaster));

offNormalMaster_deg = ...
    acosd(cosOffNormalMaster);

accessibleMaster = ...
    offNormalMaster_deg <= maxOffNormalDeg;

accessibleMasterIndex = ...
    find(accessibleMaster);

nAccessibleMasterEpochs = ...
    numel(accessibleMasterIndex);

if nAccessibleMasterEpochs < minAccessibleEpochsPerCase
    error(['The complete master schedule contains too few accessible ', ...
        'epochs for the reconstruction study.']);
end


%% Build production-compatible UV and accessibility structures

uvMaster = uvQPO;

uvMaster.uvw12 = uvw12Master;
uvMaster.uvw13 = uvw13Master;
uvMaster.uvw23 = uvw23Master;

accessMaster = accessQPO;
accessMaster.accessible = accessibleMaster;

if isfield(accessMaster,'offNormalDeg')
    accessMaster.offNormalDeg = offNormalMaster_deg;
end


%% Frozen phase-zero control over the COMPLETE accessible master schedule
%
% This is the production observation constructor and therefore also provides
% a direct audit of the sample ordering used downstream.

obsFrozenMaster = ...
    make_simulated_interferometric_observations( ...
        simulatedFrozenRot, ...
        uvMaster, ...
        accessMaster, ...
        baselineFields, ...
        baselineNames);

if ~isequal(obsFrozenMaster.accessibleIndex(:),accessibleMasterIndex(:))
    error(['Master observation bookkeeping disagrees with the ', ...
        'new accessibility mask.']);
end


%% Determine physical-sample ordering used by the production constructor

nBaselinesMaster = numel(baselineFields);

if nBaselinesMaster ~= 3
    error('The stroboscopic study currently expects exactly three baselines.');
end

uAccessibleMaster = [ ...
    uvw12Master(accessibleMaster,1), ...
    uvw13Master(accessibleMaster,1), ...
    uvw23Master(accessibleMaster,1)];

vAccessibleMaster = [ ...
    uvw12Master(accessibleMaster,2), ...
    uvw13Master(accessibleMaster,2), ...
    uvw23Master(accessibleMaster,2)];

wAccessibleMaster = [ ...
    uvw12Master(accessibleMaster,3), ...
    uvw13Master(accessibleMaster,3), ...
    uvw23Master(accessibleMaster,3)];

uBaselineMajor = uAccessibleMaster(:);
vBaselineMajor = vAccessibleMaster(:);

uEpochMajor = reshape(uAccessibleMaster.',[],1);
vEpochMajor = reshape(vAccessibleMaster.',[],1);

uObsMaster = obsFrozenMaster.u(:);
vObsMaster = obsFrozenMaster.v(:);

expectedMasterSamples = ...
    nAccessibleMasterEpochs*nBaselinesMaster;

if numel(uObsMaster) ~= expectedMasterSamples
    error(['Master observation sample count does not equal accessible ', ...
        'epochs times three baselines.']);
end

sampleOrderScale = ...
    max(1,max(abs([uObsMaster;vObsMaster])));

sampleOrderTolerance = ...
    1e-12*sampleOrderScale;

baselineMajorError = ...
    max(abs([ ...
        uBaselineMajor-uObsMaster; ...
        vBaselineMajor-vObsMaster]));

epochMajorError = ...
    max(abs([ ...
        uEpochMajor-uObsMaster; ...
        vEpochMajor-vObsMaster]));

if baselineMajorError <= epochMajorError && ...
        baselineMajorError <= sampleOrderTolerance

    stroboscopicSampleOrder = 'baseline-major';

elseif epochMajorError < baselineMajorError && ...
        epochMajorError <= sampleOrderTolerance

    stroboscopicSampleOrder = 'epoch-major';

else

    error(['Could not identify the physical-sample ordering used by ', ...
        'make_simulated_interferometric_observations.']);

end


%% Sample-level metadata aligned exactly with obsFrozenMaster.V

accessibleMasterTime_s = ...
    masterTime_s(accessibleMasterIndex);

accessibleMasterPhaseOffset_deg = ...
    masterPhaseOffset_deg(accessibleMasterIndex);

accessibleMasterPhaseIndex = ...
    masterPhaseIndex(accessibleMasterIndex);

switch stroboscopicSampleOrder

    case 'baseline-major'

        sampleMasterEpoch = ...
            repmat(accessibleMasterIndex,nBaselinesMaster,1);

        sampleBaselineIndex = [ ...
            ones(nAccessibleMasterEpochs,1); ...
            2*ones(nAccessibleMasterEpochs,1); ...
            3*ones(nAccessibleMasterEpochs,1)];

        sampleTime_s = ...
            repmat(accessibleMasterTime_s,nBaselinesMaster,1);

        samplePhaseOffset_deg = ...
            repmat(accessibleMasterPhaseOffset_deg,nBaselinesMaster,1);

        samplePhaseIndex = ...
            repmat(accessibleMasterPhaseIndex,nBaselinesMaster,1);

    case 'epoch-major'

        sampleMasterEpoch = ...
            repelem(accessibleMasterIndex,nBaselinesMaster);

        sampleBaselineIndex = ...
            repmat((1:nBaselinesMaster).',nAccessibleMasterEpochs,1);

        sampleTime_s = ...
            repelem(accessibleMasterTime_s,nBaselinesMaster);

        samplePhaseOffset_deg = ...
            repelem(accessibleMasterPhaseOffset_deg,nBaselinesMaster);

        samplePhaseIndex = ...
            repelem(accessibleMasterPhaseIndex,nBaselinesMaster);

    otherwise

        error('Unexpected stroboscopicSampleOrder.');

end

if numel(sampleTime_s) ~= numel(obsFrozenMaster.V)
    error('Sample-level master metadata has an incorrect length.');
end


%% Cache one asteroid image per unique rotational phase
%
% The target image depends on rotational phase, not on recurrence number.
% Reusing each phase image prevents repeated rasterization of the same
% physical orientation during a long campaign.

phaseLPix_rad = cell(nPhaseOffsets,1);
phaseMPix_rad = cell(nPhaseOffsets,1);
phaseFluxPix = cell(nPhaseOffsets,1);

phaseVisibleFacetCount = zeros(nPhaseOffsets,1);
phaseNonzeroPixelCount = zeros(nPhaseOffsets,1);
phaseFluxError = zeros(nPhaseOffsets,1);

fprintf('\n============================================================\n');
fprintf('PHASE-SYNCHRONOUS ROTATING-TARGET LIBRARY\n');
fprintf('============================================================\n');
fprintf('  phase spacing               = %.6f deg\n',phaseSampleStep_deg);
fprintf('  unique phase images         = %d\n',nPhaseOffsets);

for iPhase = 1:nPhaseOffsets

    phaseOffset_deg = phaseOffsetGrid_deg(iPhase);

    if iPhase == indexPhaseZero

        simulatedPhase = simulatedFrozenRot;
        phaseVisibleFacetCount(iPhase) = nnz(visFrozenRot.isVisible);

    else

        representativeTime_s = ...
            rotItokawa.t0_s + ...
            (phaseOffset_deg/360) * ...
            TItokawa_s / ...
            rotItokawa.spinSense;

        [R_IB_phase,~,~] = ...
            asteroid_attitude_at_time( ...
                rotItokawa, ...
                representativeTime_s);

        projPhase = ...
            project_shape_orthographic( ...
                shapeItokawa, ...
                R_IB_phase, ...
                target);

        visPhase = ...
            find_visible_facets(projPhase);

        imgPhase = ...
            rasterize_visible_facets( ...
                projPhase, ...
                visPhase, ...
                nPixItokawa, ...
                halfWidthRotatingItokawa_km, ...
                facetRadianceItokawa);

        simulatedPhase = ...
            make_simulated_target_from_raster( ...
                imgPhase, ...
                rangeItokawa_km, ...
                shapeItokawa.type);

        phaseVisibleFacetCount(iPhase) = ...
            nnz(visPhase.isVisible);

    end

    if numel(simulatedPhase.l_rad) ~= numel(simulatedFrozenRot.l_rad) || ...
            numel(simulatedPhase.m_rad) ~= numel(simulatedFrozenRot.m_rad)

        error('Rotating-target phase image changed the fixed image-grid dimensions.');

    end

    gridDifferencePhase = max([ ...
        max(abs(simulatedPhase.l_rad(:)-simulatedFrozenRot.l_rad(:))); ...
        max(abs(simulatedPhase.m_rad(:)-simulatedFrozenRot.m_rad(:)))]);

    if gridDifferencePhase > 1e-15
        error('Rotating-target phase image changed the fixed angular grid.');
    end

    [jPixPhase,iPixPhase,fluxPixPhase] = ...
        find(simulatedPhase.I);

    if isempty(fluxPixPhase)
        error('No nonzero target pixels at phase offset %.6f deg.',phaseOffset_deg);
    end

    lPixPhase_rad = ...
        simulatedPhase.l_rad(iPixPhase);

    mPixPhase_rad = ...
        simulatedPhase.m_rad(jPixPhase);

    phaseLPix_rad{iPhase} = lPixPhase_rad(:);
    phaseMPix_rad{iPhase} = mPixPhase_rad(:);
    phaseFluxPix{iPhase} = fluxPixPhase(:);

    phaseNonzeroPixelCount(iPhase) = numel(fluxPixPhase);
    phaseFluxError(iPhase) = abs(sum(fluxPixPhase)-1);

end

fprintf('  max flux-normalization error = %.3e\n',max(phaseFluxError));
fprintf('  visible-facet range          = %d to %d\n', ...
    min(phaseVisibleFacetCount),max(phaseVisibleFacetCount));
fprintf('  nonzero-pixel range          = %d to %d\n', ...
    min(phaseNonzeroPixelCount),max(phaseNonzeroPixelCount));
fprintf('============================================================\n');


%% Rotating-target visibilities on exactly the same master UV samples
%
% VRotAccessibleMaster is stored as [accessible epoch x baseline].
% Each phase image is reused for every recurrence having that same phase.

VRotAccessibleMaster = ...
    complex(nan(nAccessibleMasterEpochs,nBaselinesMaster));

fprintf('\nGenerating phase-synchronous rotating visibilities...\n');

ticStroboscopicVisibility = tic;

for iPhase = 1:nPhaseOffsets

    accessibleRowsThisPhase = ...
        find(accessibleMasterPhaseIndex == iPhase);

    if isempty(accessibleRowsThisPhase)
        continue
    end

    lPix = phaseLPix_rad{iPhase};
    mPix = phaseMPix_rad{iPhase};
    fluxPix = phaseFluxPix{iPhase};

    for b = 1:nBaselinesMaster

        uThis = ...
            uAccessibleMaster(accessibleRowsThisPhase,b);

        vThis = ...
            vAccessibleMaster(accessibleRowsThisPhase,b);

        nRowsThisPhase = numel(accessibleRowsThisPhase);

        for iBlockStart = 1:visibilityBlockSize:nRowsThisPhase

            iBlockEnd = ...
                min(iBlockStart+visibilityBlockSize-1,nRowsThisPhase);

            localRows = ...
                iBlockStart:iBlockEnd;

            physicalRows = ...
                accessibleRowsThisPhase(localRows);

            fourierPhase = ...
                uThis(localRows)*lPix.' + ...
                vThis(localRows)*mPix.';

            VRotAccessibleMaster(physicalRows,b) = ...
                exp(-2*pi*1i*fourierPhase) * fluxPix;

        end

    end

end

stroboscopicVisibilityRuntime_s = ...
    toc(ticStroboscopicVisibility);

if any(~isfinite(real(VRotAccessibleMaster)),'all') || ...
        any(~isfinite(imag(VRotAccessibleMaster)),'all')

    error('At least one master rotating visibility was not evaluated.');

end


%% Flatten rotating visibilities into the production sample order

switch stroboscopicSampleOrder

    case 'baseline-major'

        VRotMasterVector = ...
            VRotAccessibleMaster(:);

    case 'epoch-major'

        VRotMasterVector = ...
            reshape(VRotAccessibleMaster.',[],1);

    otherwise

        error('Unexpected stroboscopicSampleOrder.');

end

VFrozenMasterVector = ...
    obsFrozenMaster.V(:);

if numel(VRotMasterVector) ~= numel(VFrozenMasterVector)
    error('Rotating and frozen master visibility vectors have different lengths.');
end


%% Exact-phase forward-model audit
%
% At DeltaPhi = 0 the rotating image is deliberately the same fixed-grid
% phase-zero image used by the frozen control.  Therefore the two
% visibility vectors must agree to numerical precision at the same UV points.

zeroPhaseSampleMask = ...
    abs(samplePhaseOffset_deg) <= phaseGridTolerance_deg;

if ~any(zeroPhaseSampleMask)
    error('The accessible master data contain no exact-phase samples.');
end

zeroPhaseVisibilityDifference = ...
    VRotMasterVector(zeroPhaseSampleMask) - ...
    VFrozenMasterVector(zeroPhaseSampleMask);

maxZeroPhaseVisibilityError = ...
    max(abs(zeroPhaseVisibilityDifference));

rmsZeroPhaseVisibilityError = ...
    sqrt(mean(abs(zeroPhaseVisibilityDifference).^2));

if maxZeroPhaseVisibilityError > zeroPhaseVisibilityTolerance
    error(['Exact-phase rotating and frozen visibilities disagree. ', ...
        'Max |DeltaV| = %.3e.'], ...
        maxZeroPhaseVisibilityError);
end


%% Master-study summary

fprintf('\n============================================================\n');
fprintf('PHASE-SYNCHRONOUS NATURAL-QPO MASTER DATA SET\n');
fprintf('============================================================\n');
fprintf('  QPO duration                 = %.6f days\n',qpoDuration_s/86400);
fprintf('  asteroid period              = %.6f h\n',TItokawa_s/3600);
fprintf('  phase-grid spacing           = %.6f deg\n',phaseSampleStep_deg);
fprintf('  cadence represented          = %.6f min per phase step\n', ...
    phaseSampleStep_deg/360*TItokawa_s/60);
fprintf('  master epochs                = %d\n',nMasterEpochs);
fprintf('  accessible master epochs     = %d\n',nAccessibleMasterEpochs);
fprintf('  physical master samples      = %d\n',numel(VFrozenMasterVector));
fprintf('  sample ordering              = %s\n',stroboscopicSampleOrder);
fprintf('  max requested-phase error    = %.3e deg\n',maxPhaseScheduleError_deg);
fprintf('  min master time spacing      = %.6f min\n',minimumMasterTimeSpacing_s/60);
fprintf('  collector dense-grid errors  = %.3e %.3e %.3e\n', ...
    maxDenseGridDifference(1), ...
    maxDenseGridDifference(2), ...
    maxDenseGridDifference(3));
fprintf('  synodic baseline closure     = %.3e nondimensional\n', ...
    maxBaselineClosureSynodicMaster);
fprintf('  UVW closure                  = %.3e wavelengths\n', ...
    maxUVWClosureMaster);
fprintf('  exact-phase max |Vrot-V0|    = %.3e\n', ...
    maxZeroPhaseVisibilityError);
fprintf('  exact-phase RMS |Vrot-V0|    = %.3e\n', ...
    rmsZeroPhaseVisibilityError);
fprintf('  rotating-visibility runtime  = %.3f s\n', ...
    stroboscopicVisibilityRuntime_s);
fprintf('============================================================\n');


%% Exact-phase UV coverage figure

figExactPhaseUV = figure( ...
    'Name','Exact-phase stroboscopic natural-QPO UV coverage', ...
    'Color','w', ...
    'Position',[120 100 850 750]);

axExactPhaseUV = axes(figExactPhaseUV);
hold(axExactPhaseUV,'on');
grid(axExactPhaseUV,'on');
box(axExactPhaseUV,'on');

baselineColorsStrobe = lines(3);

exactAccessibleMaster = ...
    accessibleMaster & ...
    abs(masterPhaseOffset_deg) <= phaseGridTolerance_deg;

for b = 1:3

    uvwThis = uvMaster.(baselineFields{b});

    scatter( ...
        axExactPhaseUV, ...
        uvwThis(exactAccessibleMaster,1)/1e9, ...
        uvwThis(exactAccessibleMaster,2)/1e9, ...
        36, ...
        baselineColorsStrobe(b,:), ...
        'filled', ...
        'DisplayName',baselineNames{b});

end

axis(axExactPhaseUV,'equal');

xlabel(axExactPhaseUV,'$u$ [G$\lambda$]');
ylabel(axExactPhaseUV,'$v$ [G$\lambda$]');

title( ...
    axExactPhaseUV, ...
    'Natural QPO sampled at repeated phase-zero asteroid aspect', ...
    'Interpreter','none');

legend( ...
    axExactPhaseUV, ...
    'Location','best', ...
    'Interpreter','none');


%% ============================================================
%  22. END-TO-END RECONSTRUCTION SWEEP
%  ============================================================
%
% Every case uses:
%
%   1. the SAME physical timestamps and UV samples for the rotating target
%      and the frozen phase-zero control;
%
%   2. the SAME weighting rule and Högbom CLEAN implementation;
%
%   3. a full phase-bin width W centered on phase zero.
%
% The native CLEAN beam is stored here but is NOT used as the final
% comparison resolution.  Section 23 converts every comparable case to one
% common analysis beam.

caseSuccess = false(nCampaignCases,nBinCases);
caseAccessibleEpochs = zeros(nCampaignCases,nBinCases);
casePhysicalSamples = zeros(nCampaignCases,nBinCases);

caseBeamMajor_mas = nan(nCampaignCases,nBinCases);
caseBeamMinor_mas = nan(nCampaignCases,nBinCases);
caseBeamPA_deg = nan(nCampaignCases,nBinCases);

caseDirtyBeamDifference = nan(nCampaignCases,nBinCases);
caseCleanBeamKernelDifference = nan(nCampaignCases,nBinCases);

caseFrozenIterations = nan(nCampaignCases,nBinCases);
caseRotatingIterations = nan(nCampaignCases,nBinCases);

caseImageFrozenNative = cell(nCampaignCases,nBinCases);
caseImageRotatingNative = cell(nCampaignCases,nBinCases);

fprintf('\n============================================================\n');
fprintf('STROBOSCOPIC RECONSTRUCTION SWEEP\n');
fprintf('============================================================\n');
fprintf('  campaign cases = %d\n',nCampaignCases);
fprintf('  bin-width cases = %d\n',nBinCases);
fprintf('  minimum accessible epochs per case = %d\n', ...
    minAccessibleEpochsPerCase);
fprintf('============================================================\n');

for iCampaign = 1:nCampaignCases

    campaignStart_s = tQPOStart_s;
    campaignEnd_s = ...
        campaignStart_s + campaignDuration_s(iCampaign);

    for iBin = 1:nBinCases

        fullBinWidth_deg = phaseBinWidth_deg(iBin);
        halfBinWidth_deg = 0.5*fullBinWidth_deg;

        inCampaignSample = ...
            sampleTime_s >= campaignStart_s - 1e-9 & ...
            sampleTime_s <= campaignEnd_s + 1e-9;

        if fullBinWidth_deg == 0

            inPhaseBinSample = ...
                abs(samplePhaseOffset_deg) <= ...
                phaseGridTolerance_deg;

        else

            inPhaseBinSample = ...
                abs(samplePhaseOffset_deg) <= ...
                halfBinWidth_deg + phaseGridTolerance_deg;

        end

        keepSample = ...
            inCampaignSample & ...
            inPhaseBinSample;

        selectedMasterEpochs = ...
            unique(sampleMasterEpoch(keepSample),'stable');

        nSelectedEpochs = ...
            numel(selectedMasterEpochs);

        nSelectedSamples = ...
            nnz(keepSample);

        caseAccessibleEpochs(iCampaign,iBin) = ...
            nSelectedEpochs;

        casePhysicalSamples(iCampaign,iBin) = ...
            nSelectedSamples;

        fprintf('\n  campaign %.3f TPrimary (%.3f d), W = %.1f deg\n', ...
            campaignFractions(iCampaign), ...
            campaignDuration_days(iCampaign), ...
            fullBinWidth_deg);

        fprintf('    accessible epochs     = %d\n',nSelectedEpochs);
        fprintf('    physical samples      = %d\n',nSelectedSamples);

        if nSelectedEpochs < minAccessibleEpochsPerCase

            fprintf('    skipped: too few accessible epochs\n');
            continue

        end

        try

            %% Frozen-control observation structure
            %
            % The production master structure is copied, then every
            % sample-level top-level field used by imaging is restricted to
            % the exact same physical sample mask.

            obsFrozenCase = obsFrozenMaster;

            obsFrozenCase.u = ...
                obsFrozenMaster.u(keepSample);

            obsFrozenCase.v = ...
                obsFrozenMaster.v(keepSample);

            if isfield(obsFrozenMaster,'w') && ...
                    numel(obsFrozenMaster.w) == numel(VFrozenMasterVector)

                obsFrozenCase.w = ...
                    obsFrozenMaster.w(keepSample);

            end

            obsFrozenCase.V = ...
                VFrozenMasterVector(keepSample);

            obsFrozenCase.rho = ...
                hypot(obsFrozenCase.u,obsFrozenCase.v);

            obsFrozenCase.rho_Glambda = ...
                obsFrozenCase.rho/1e9;

            obsFrozenCase.accessibleIndex = ...
                selectedMasterEpochs(:);

            if isfield(obsFrozenCase,'amplitude')
                obsFrozenCase.amplitude = abs(obsFrozenCase.V);
            end

            if isfield(obsFrozenCase,'phase')
                obsFrozenCase.phase = angle(obsFrozenCase.V);
            end

            if isfield(obsFrozenCase,'phase_rad')
                obsFrozenCase.phase_rad = angle(obsFrozenCase.V);
            end

            if isfield(obsFrozenCase,'phase_deg')
                obsFrozenCase.phase_deg = rad2deg(angle(obsFrozenCase.V));
            end


            %% Rotating-target observation with IDENTICAL physical sampling

            obsRotatingCase = obsFrozenCase;

            obsRotatingCase.V = ...
                VRotMasterVector(keepSample);

            if isfield(obsRotatingCase,'amplitude')
                obsRotatingCase.amplitude = abs(obsRotatingCase.V);
            end

            if isfield(obsRotatingCase,'phase')
                obsRotatingCase.phase = angle(obsRotatingCase.V);
            end

            if isfield(obsRotatingCase,'phase_rad')
                obsRotatingCase.phase_rad = angle(obsRotatingCase.V);
            end

            if isfield(obsRotatingCase,'phase_deg')
                obsRotatingCase.phase_deg = ...
                    rad2deg(angle(obsRotatingCase.V));
            end


            %% Frozen dirty image

            dirtyFrozenCase = ...
                make_dirty_image_from_observations( ...
                    obsFrozenCase, ...
                    simulatedFrozenRot.l_rad, ...
                    simulatedFrozenRot.m_rad, ...
                    addHermitian, ...
                    weightMode);


            %% Rotating dirty image

            dirtyRotatingCase = ...
                make_dirty_image_from_observations( ...
                    obsRotatingCase, ...
                    simulatedFrozenRot.l_rad, ...
                    simulatedFrozenRot.m_rad, ...
                    addHermitian, ...
                    weightMode);


            %% Same-sampling dirty-beam audit
            %
            % The dirty beam depends on UV coordinates and weights, not on
            % the target visibility values.  It therefore must be the same
            % for these two data sets.

            dirtyBeamDifference = ...
                max(abs( ...
                    dirtyFrozenCase.beam(:) - ...
                    dirtyRotatingCase.beam(:)));

            caseDirtyBeamDifference(iCampaign,iBin) = ...
                dirtyBeamDifference;

            if dirtyBeamDifference > beamEqualityTolerance
                error(['Frozen and rotating controls produced different ', ...
                    'dirty beams under identical sampling. Max difference ', ...
                    '= %.3e.'], ...
                    dirtyBeamDifference);
            end


            %% CLEAN parameter rule for the frozen control

            cleanParamsFrozen = struct();

            cleanParamsFrozen.loopGain = cleanLoopGain;
            cleanParamsFrozen.nIterMax = cleanNIterMax;
            cleanParamsFrozen.cleanBeamFitLevel = cleanBeamFitLevel;
            cleanParamsFrozen.cleanBeamTruncateSigma = ...
                cleanBeamTruncateSigma;

            dirtyBeamPeakFrozen = ...
                max(real(dirtyFrozenCase.beam(:)));

            if dirtyBeamPeakFrozen <= 0
                error('Frozen-control dirty beam has a nonpositive peak.');
            end

            dirtyImageNormalizedFrozen = ...
                real(dirtyFrozenCase.image)/dirtyBeamPeakFrozen;

            cleanParamsFrozen.stopThreshold = ...
                cleanStopFraction * ...
                max(abs(dirtyImageNormalizedFrozen(:)));


            %% Frozen CLEAN

            cleanFrozenCase = ...
                make_clean_image_from_dirty( ...
                    dirtyFrozenCase, ...
                    cleanParamsFrozen);


            %% CLEAN parameter rule for the rotating data

            cleanParamsRotating = cleanParamsFrozen;

            dirtyBeamPeakRotating = ...
                max(real(dirtyRotatingCase.beam(:)));

            if dirtyBeamPeakRotating <= 0
                error('Rotating-target dirty beam has a nonpositive peak.');
            end

            dirtyImageNormalizedRotating = ...
                real(dirtyRotatingCase.image)/dirtyBeamPeakRotating;

            cleanParamsRotating.stopThreshold = ...
                cleanStopFraction * ...
                max(abs(dirtyImageNormalizedRotating(:)));


            %% Rotating CLEAN

            cleanRotatingCase = ...
                make_clean_image_from_dirty( ...
                    dirtyRotatingCase, ...
                    cleanParamsRotating);


            %% Same-sampling fitted-beam audit

            if ~isequal( ...
                    size(cleanFrozenCase.cleanBeamKernel), ...
                    size(cleanRotatingCase.cleanBeamKernel))

                error(['Frozen and rotating CLEAN kernels have different ', ...
                    'dimensions under identical UV sampling.']);

            end

            cleanBeamKernelDifference = ...
                max(abs( ...
                    cleanFrozenCase.cleanBeamKernel(:) - ...
                    cleanRotatingCase.cleanBeamKernel(:)));

            caseCleanBeamKernelDifference(iCampaign,iBin) = ...
                cleanBeamKernelDifference;

            if cleanBeamKernelDifference > beamEqualityTolerance
                error(['Frozen and rotating controls produced different ', ...
                    'CLEAN restoring kernels under identical sampling. ', ...
                    'Max difference = %.3e.'], ...
                    cleanBeamKernelDifference);
            end


            %% Store native reconstructions and beam
            %
            % Use restoredComponents consistently with the earlier
            % morphology-validation sections.  Section 23 will bring these
            % images to one fixed common analysis resolution.

            caseImageFrozenNative{iCampaign,iBin} = ...
                real(cleanFrozenCase.restoredComponents);

            caseImageRotatingNative{iCampaign,iBin} = ...
                real(cleanRotatingCase.restoredComponents);

            caseBeamMajor_mas(iCampaign,iBin) = ...
                cleanFrozenCase.params.cleanBeamFwhmMajor_mas;

            caseBeamMinor_mas(iCampaign,iBin) = ...
                cleanFrozenCase.params.cleanBeamFwhmMinor_mas;

            caseBeamPA_deg(iCampaign,iBin) = ...
                cleanFrozenCase.params.cleanBeamPaDeg;

            caseFrozenIterations(iCampaign,iBin) = ...
                cleanFrozenCase.nIter;

            caseRotatingIterations(iCampaign,iBin) = ...
                cleanRotatingCase.nIter;

            caseSuccess(iCampaign,iBin) = true;

            fprintf('    native beam          = %.4f x %.4f mas\n', ...
                caseBeamMajor_mas(iCampaign,iBin), ...
                caseBeamMinor_mas(iCampaign,iBin));

            fprintf('    dirty-beam match     = %.3e\n', ...
                dirtyBeamDifference);

            fprintf('    CLEAN-beam match     = %.3e\n', ...
                cleanBeamKernelDifference);

        catch ME

            warning( ...
                ['Stroboscopic reconstruction case failed at ', ...
                 'campaign %.3f TPrimary and W = %.1f deg:\n%s'], ...
                campaignFractions(iCampaign), ...
                fullBinWidth_deg, ...
                ME.message);

        end

    end

end


%% ============================================================
%  23. FIXED ANALYSIS RESOLUTION + TARGET-RESOLUTION SWEEP SETUP
%  ============================================================
%
% This section replaces the previous single-range common-resolution study.
%
% Scientific purpose
% ------------------
% The previous run established that exact-phase stroboscopic synthesis works
% for the current noiseless geometry, but the phase-zero Itokawa truth was
% only weakly resolved at the common analysis resolution.  As a result, even
% a full 360-deg rotational bin could still satisfy loose morphology metrics.
%
% We now vary asteroid range while keeping fixed:
%
%   - the natural QPO and all physical UV samples;
%   - wavelength;
%   - the physical Gaskell shape and phase-zero attitude;
%   - the phase-bin definitions;
%   - CLEAN settings;
%   - one angular analysis beam for the entire study.
%
% Changing range therefore changes the number of resolved target elements
% without changing the formation dynamics.
%
% Assessment states are explicitly separated into:
%
%   PASS       : reconstructed and satisfies all selected fidelity criteria;
%   FAIL       : reconstructed at the common resolution but violates at least
%                one selected criterion;
%   UNASSESSED : cannot be compared at the chosen common resolution, or the
%                reconstruction failed numerically.
%
% IMPORTANT:
% The thresholds below remain user-selected STUDY criteria.  They are not
% fundamental information-theory recoverability limits.


%% Locate the exact-phase bin

[minimumZeroBinWidth_deg,iZeroBin] = ...
    min(abs(phaseBinWidth_deg));

if minimumZeroBinWidth_deg > phaseGridTolerance_deg
    error('Could not locate the W = 0 exact-phase bin.');
end


%% FAST PILOT SUBSET
%
% First locate the blur transition using only the longest campaign and a
% reduced phase-bin grid.  This cuts the range sweep from
%
%   4 ranges x 5 campaigns x 11 widths = 220 cases
%
% to
%
%   4 ranges x 1 campaign x 8 widths = 32 cases.
%
% Once a transition is found, refine only the relevant range / width region.

pilotCampaignIndices = ...
    nCampaignCases;

pilotBinWidths_deg = [ ...
      0 ...
     30 ...
     60 ...
     90 ...
    120 ...
    180 ...
    240 ...
    360];

pilotBinIndices = ...
    nan(size(pilotBinWidths_deg));

for iPilotBin = 1:numel(pilotBinWidths_deg)

    [distancePilotBin_deg,indexPilotBin] = ...
        min(abs(phaseBinWidth_deg-pilotBinWidths_deg(iPilotBin)));

    if distancePilotBin_deg > phaseGridTolerance_deg

        error( ...
            'Requested pilot bin width %.6f deg is not in phaseBinWidth_deg.', ...
            pilotBinWidths_deg(iPilotBin));

    end

    pilotBinIndices(iPilotBin) = ...
        indexPilotBin;

end

pilotBinIndices = ...
    unique(pilotBinIndices,'stable');

pilotBinWidths_deg = ...
    phaseBinWidth_deg(pilotBinIndices);

if pilotBinIndices(1) ~= iZeroBin
    error('The pilot bin grid must begin with W = 0 deg.');
end

nPilotCampaigns = ...
    numel(pilotCampaignIndices);

nPilotBins = ...
    numel(pilotBinIndices);




%% Choose one angular analysis beam for every range and every tested case
%
% Recommended default:
%
%   'referenceExactPhase'
%
% uses the longest-campaign exact-phase native beam.  Short campaigns whose
% native beam is broader are then correctly labeled UNASSESSED rather than
% FAIL.
%
% Alternative:
%
%   'allSuccessfulCases'
%
% deliberately chooses a broader circular beam that dominates every
% successful native beam from Section 22.  This permits more cases to be
% compared, at the cost of lower common angular resolution.

analysisBeamMode = 'referenceExactPhase';

analysisBeamSafetyFactor = 1.02;

iReferenceCampaign = nCampaignCases;

switch lower(analysisBeamMode)

    case 'referenceexactphase'

        if ~caseSuccess(iReferenceCampaign,iZeroBin)
            error(['The longest-campaign exact-phase case did not ', ...
                'reconstruct successfully in Section 22.']);
        end

        analysisBeamFwhm_mas = ...
            analysisBeamSafetyFactor * ...
            caseBeamMajor_mas(iReferenceCampaign,iZeroBin);


    case 'allsuccessfulcases'

        successfulMajorAxes = ...
            caseBeamMajor_mas(caseSuccess & isfinite(caseBeamMajor_mas));

        if isempty(successfulMajorAxes)
            error('No successful native CLEAN beams are available.');
        end

        analysisBeamFwhm_mas = ...
            analysisBeamSafetyFactor * ...
            max(successfulMajorAxes);


    otherwise

        error('Unknown analysisBeamMode: %s',analysisBeamMode);

end

fwhmToSigma = ...
    1/(2*sqrt(2*log(2)));

analysisBeamSigma_mas = ...
    fwhmToSigma*analysisBeamFwhm_mas;

SigmaAnalysis_mas2 = ...
    analysisBeamSigma_mas^2*eye(2);

covarianceTolerance = ...
    1e-12*max(1,norm(SigmaAnalysis_mas2,'fro'));


%% Resolution comparability is assessed inside each range case
%
% The underlying dirty beam is set by the UV samples, but the fitted CLEAN
% restoring beam is estimated from that beam on the reconstruction image
% grid.  Because the angular pixel spacing changes with target range, the
% fitted FWHM can differ slightly from the Section-22 value even though the
% physical UV coverage is unchanged.  Therefore we do NOT require exact
% equality with the base-range fitted beam.  Instead, each range/campaign/bin
% case computes its own fitted native-beam covariance and asks only whether
% that covariance can be blurred to the fixed analysis covariance.


%% Target-range sweep
%
% The current range is included first, followed by progressively closer
% ranges.  The physical raster is unchanged.  Only the angular scale changes.
%
% This broad logarithmic sweep is intended to locate where rotational blur
% starts to become reconstruction-limiting.  Refine the range list later if
% a transition appears between two neighboring values.

rangeScaleStudy = [ ...
    1.000 ...
    0.500 ...
    0.250 ...
    0.125];

rangeStudy_km = ...
    rangeItokawa_km*rangeScaleStudy;

nRangeCases = ...
    numel(rangeStudy_km);


fprintf('\n============================================================\n');
fprintf('FAST TARGET-RESOLUTION PILOT\n');
fprintf('============================================================\n');
fprintf('  campaign                  = %.3f TPrimary (%.3f d)\n', ...
    campaignFractions(pilotCampaignIndices(1)), ...
    campaignDuration_days(pilotCampaignIndices(1)));
fprintf('  tested full bin widths    = ');
fprintf('%.0f ',pilotBinWidths_deg);
fprintf('deg\n');
fprintf('  range cases               = %d\n',numel(rangeScaleStudy));
fprintf('  reconstruction cases      = %d\n', ...
    numel(rangeScaleStudy)*nPilotCampaigns*nPilotBins);
fprintf('  CLEAN runs                = %d\n', ...
    2*numel(rangeScaleStudy)*nPilotCampaigns*nPilotBins);
fprintf('============================================================\n');


%% Representative finite-bin widths to retain for image inspection

representativeWidths_deg = [0 60 180 360];

representativeBinIndex = ...
    nan(size(representativeWidths_deg));

for iRep = 1:numel(representativeWidths_deg)

    [distanceRep,idxRep] = ...
        min(abs(phaseBinWidth_deg-representativeWidths_deg(iRep)));

    if distanceRep <= phaseGridTolerance_deg
        representativeBinIndex(iRep) = idxRep;
    end

end

representativeBinIndex = ...
    representativeBinIndex(isfinite(representativeBinIndex));

representativeBinIndex = ...
    unique(representativeBinIndex,'stable');

nRepresentativeBins = ...
    numel(representativeBinIndex);


%% Allocate target-resolution diagnostics

projectedMajorSpan_mas = ...
    nan(nRangeCases,1);

projectedMinorSpan_mas = ...
    nan(nRangeCases,1);

resolvedElementsMajor = ...
    nan(nRangeCases,1);

resolvedElementsMinor = ...
    nan(nRangeCases,1);

rangePixelScaleL_mas = ...
    nan(nRangeCases,1);

rangePixelScaleM_mas = ...
    nan(nRangeCases,1);

rangeExactPhaseVisibilityError = ...
    nan(nRangeCases,1);


%% Allocate range x campaign x bin reconstruction metrics

rangeCaseSuccess = ...
    false(nRangeCases,nCampaignCases,nBinCases);

rangeCaseAssessed = ...
    false(nRangeCases,nCampaignCases,nBinCases);

rangeCaseResolutionComparable = ...
    false(nRangeCases,nCampaignCases,nBinCases);

rangeNativeBeamMajor_mas = ...
    nan(nRangeCases,nCampaignCases,nBinCases);

rangeNativeBeamMinor_mas = ...
    nan(nRangeCases,nCampaignCases,nBinCases);

rangeNativeBeamPA_deg = ...
    nan(nRangeCases,nCampaignCases,nBinCases);

rangeNativeBeamMajorShift_mas = ...
    nan(nRangeCases,nCampaignCases,nBinCases);

rangeNativeBeamMinorShift_mas = ...
    nan(nRangeCases,nCampaignCases,nBinCases);

% Assessment code:
%   0 = UNASSESSED
%   1 = FAIL
%   2 = PASS
rangeCaseAssessmentCode = ...
    zeros(nRangeCases,nCampaignCases,nBinCases,'uint8');

rangeNRMSEFrozen = ...
    nan(nRangeCases,nCampaignCases,nBinCases);

rangeNRMSERotating = ...
    nan(nRangeCases,nCampaignCases,nBinCases);

rangeNCCFrozen = ...
    nan(nRangeCases,nCampaignCases,nBinCases);

rangeNCCRotating = ...
    nan(nRangeCases,nCampaignCases,nBinCases);

rangeIoUFrozen = ...
    nan(nRangeCases,nCampaignCases,nBinCases);

rangeIoURotating = ...
    nan(nRangeCases,nCampaignCases,nBinCases);

rangeDeltaNRMSEBlur = ...
    nan(nRangeCases,nCampaignCases,nBinCases);

rangeDeltaNCCBlur = ...
    nan(nRangeCases,nCampaignCases,nBinCases);

rangeDeltaIoUBlur = ...
    nan(nRangeCases,nCampaignCases,nBinCases);


%% Store only a small number of representative images
%
% Avoid retaining hundreds of 401 x 401 reconstruction arrays.

rangeTruthAnalysis = ...
    cell(nRangeCases,1);

representativeFrozenAnalysis = ...
    cell(nRangeCases,nRepresentativeBins);

representativeRotatingAnalysis = ...
    cell(nRangeCases,nRepresentativeBins);


%% Section-22 native images are no longer needed
%
% Keep the beam metadata, success masks, and sample counts, but free the
% large native-image cell arrays before the range sweep.

if exist('caseImageFrozenNative','var')
    clear caseImageFrozenNative
end

if exist('caseImageRotatingNative','var')
    clear caseImageRotatingNative
end


fprintf('\n============================================================\n');
fprintf('TARGET-RESOLUTION / RANGE SWEEP SETUP\n');
fprintf('============================================================\n');
fprintf('  analysis beam mode           = %s\n',analysisBeamMode);
fprintf('  fixed analysis FWHM          = %.6f mas\n',analysisBeamFwhm_mas);
fprintf('  tested ranges [km]           = ');

for iRange = 1:nRangeCases
    fprintf('%.6e ',rangeStudy_km(iRange));
end

fprintf('\n');
fprintf('  range scale relative to base = ');

for iRange = 1:nRangeCases
    fprintf('%.3f ',rangeScaleStudy(iRange));
end

fprintf('\n');
fprintf('============================================================\n');


%% ============================================================
%  24. END-TO-END RANGE x CAMPAIGN x PHASE-BIN RECONSTRUCTION SWEEP
%  ============================================================
%
% For every target range:
%
%   1. generate the phase-zero angular target on that range;
%   2. rescale the cached rotating phase-library coordinates exactly as 1/D;
%   3. recompute rotating and frozen visibilities on the SAME natural-QPO
%      UV samples;
%   4. reconstruct every campaign / phase-bin case;
%   5. bring both reconstructions to the SAME fixed angular analysis beam;
%   6. compare against the SAME phase-zero truth at that beam.
%
% The frozen control always uses the identical sample mask as the rotating
% target.  Differences therefore isolate the effect of target rotation.

fprintf('\n============================================================\n');
fprintf('TARGET-RESOLUTION / RANGE RECONSTRUCTION SWEEP\n');
fprintf('============================================================\n');

for iRange = 1:nRangeCases

    thisRange_km = ...
        rangeStudy_km(iRange);

    angularScaleFactor = ...
        rangeItokawa_km/thisRange_km;

    fprintf('\n------------------------------------------------------------\n');
    fprintf('RANGE CASE %d / %d\n',iRange,nRangeCases);
    fprintf('  range                        = %.6e km\n',thisRange_km);
    fprintf('  angular scale / base         = %.6f\n',angularScaleFactor);
    fprintf('------------------------------------------------------------\n');


    %% Phase-zero target at this range

    simulatedRange = ...
        make_simulated_target_from_raster( ...
            imgFrozenRot, ...
            thisRange_km, ...
            shapeItokawa.type);

    imageRangeDifference = ...
        max(abs( ...
            simulatedRange.I(:) - ...
            simulatedFrozenRot.I(:)));

    if imageRangeDifference > 1e-12
        error(['Changing range unexpectedly changed the physical ', ...
            'fractional-flux raster. Max difference = %.3e.'], ...
            imageRangeDifference);
    end


    %% Principal projected target spans at phase zero
    %
    % Use the nonzero image support, rotate those sky-plane coordinates into
    % their principal axes, and measure the full support spans.  Because the
    % common analysis beam is circular, the resulting ratios
    %
    %   Nres,major = major angular span / analysis FWHM
    %   Nres,minor = minor angular span / analysis FWHM
    %
    % provide a transparent measure of how resolved the target is.

    [jSupport,iSupport] = ...
        find(simulatedRange.I > 0);

    if numel(iSupport) < 3
        error('Phase-zero target support contains fewer than three pixels.');
    end

    lSupport_mas = ...
        simulatedRange.l_mas(iSupport);

    mSupport_mas = ...
        simulatedRange.m_mas(jSupport);

    lSupport_mas = ...
        lSupport_mas(:);

    mSupport_mas = ...
        mSupport_mas(:);

    supportCoordinates_mas = [ ...
        lSupport_mas, ...
        mSupport_mas];

    supportCentroid_mas = ...
        mean(supportCoordinates_mas,1);

    supportCentered_mas = ...
        supportCoordinates_mas-supportCentroid_mas;

    supportCovariance = ...
        (supportCentered_mas.'*supportCentered_mas) / ...
        max(size(supportCentered_mas,1)-1,1);

    [supportEigVec,supportEigVal] = ...
        eig(0.5*(supportCovariance+supportCovariance.'));

    [~,supportOrder] = ...
        sort(real(diag(supportEigVal)),'descend');

    supportPrincipalCoordinates_mas = ...
        supportCentered_mas * ...
        supportEigVec(:,supportOrder);

    supportSpan_mas = ...
        max(supportPrincipalCoordinates_mas,[],1) - ...
        min(supportPrincipalCoordinates_mas,[],1);

    projectedMajorSpan_mas(iRange) = ...
        max(supportSpan_mas);

    projectedMinorSpan_mas(iRange) = ...
        min(supportSpan_mas);

    resolvedElementsMajor(iRange) = ...
        projectedMajorSpan_mas(iRange) / ...
        analysisBeamFwhm_mas;

    resolvedElementsMinor(iRange) = ...
        projectedMinorSpan_mas(iRange) / ...
        analysisBeamFwhm_mas;

    rangePixelScaleL_mas(iRange) = ...
        mean(diff(simulatedRange.l_mas));

    rangePixelScaleM_mas(iRange) = ...
        mean(diff(simulatedRange.m_mas));


    %% Fixed analysis beam and phase-zero truth at this range

    analysisBeamKernelRange = ...
        gaussian_kernel_from_covariance_mas( ...
            SigmaAnalysis_mas2, ...
            rangePixelScaleL_mas(iRange), ...
            rangePixelScaleM_mas(iRange), ...
            5.0);

    truthAnalysisRange = ...
        real(conv2( ...
            simulatedRange.I, ...
            analysisBeamKernelRange, ...
            'same'));

    peakTruthAnalysisRange = ...
        max(truthAnalysisRange(:));

    if peakTruthAnalysisRange <= 0
        error('Analysis-resolution truth has no positive peak.');
    end

    supportAnalysisRange = ...
        truthAnalysisRange >= ...
        correlationSupportLevel*peakTruthAnalysisRange;

    if nnz(supportAnalysisRange) < 2
        error('Analysis-resolution support mask contains too few pixels.');
    end

    rangeTruthAnalysis{iRange} = ...
        truthAnalysisRange;


    %% Frozen phase-zero master observation at this range

    obsFrozenRangeMaster = ...
        make_simulated_interferometric_observations( ...
            simulatedRange, ...
            uvMaster, ...
            accessMaster, ...
            baselineFields, ...
            baselineNames);

    if ~isequal( ...
            obsFrozenRangeMaster.accessibleIndex(:), ...
            accessibleMasterIndex(:))

        error(['Range-sweep frozen observation changed the master ', ...
            'accessible-epoch bookkeeping.']);

    end

    if numel(obsFrozenRangeMaster.V) ~= numel(sampleTime_s)
        error('Range-sweep frozen master observation has an incorrect sample count.');
    end

    uvDifferenceRange = ...
        max(abs([ ...
            obsFrozenRangeMaster.u(:)-obsFrozenMaster.u(:); ...
            obsFrozenRangeMaster.v(:)-obsFrozenMaster.v(:)]));

    if uvDifferenceRange > sampleOrderTolerance
        error(['Changing target range unexpectedly changed the UV ', ...
            'sample ordering or coordinates.']);
    end

    VFrozenRangeVector = ...
        obsFrozenRangeMaster.V(:);


    %% Rotating-target master visibilities at this range
    %
    % The physical raster and fractional fluxes do not change with range.
    % Orthographic angular coordinates scale exactly as 1/D.  Therefore the
    % already validated sparse phase-library coordinates can be rescaled
    % without re-rasterizing the Gaskell model.

    VRotAccessibleRange = ...
        complex(nan(nAccessibleMasterEpochs,nBaselinesMaster));

    ticRangeVisibility = tic;

    for iPhase = 1:nPhaseOffsets

        accessibleRowsThisPhase = ...
            find(accessibleMasterPhaseIndex == iPhase);

        if isempty(accessibleRowsThisPhase)
            continue
        end

        lPixRange_rad = ...
            angularScaleFactor * ...
            phaseLPix_rad{iPhase};

        mPixRange_rad = ...
            angularScaleFactor * ...
            phaseMPix_rad{iPhase};

        fluxPixRange = ...
            phaseFluxPix{iPhase};

        for b = 1:nBaselinesMaster

            uThis = ...
                uAccessibleMaster(accessibleRowsThisPhase,b);

            vThis = ...
                vAccessibleMaster(accessibleRowsThisPhase,b);

            nRowsThisPhase = ...
                numel(accessibleRowsThisPhase);

            for iBlockStart = 1:visibilityBlockSize:nRowsThisPhase

                iBlockEnd = ...
                    min( ...
                        iBlockStart+visibilityBlockSize-1, ...
                        nRowsThisPhase);

                localRows = ...
                    iBlockStart:iBlockEnd;

                physicalRows = ...
                    accessibleRowsThisPhase(localRows);

                fourierPhaseRange = ...
                    uThis(localRows)*lPixRange_rad.' + ...
                    vThis(localRows)*mPixRange_rad.';

                VRotAccessibleRange(physicalRows,b) = ...
                    exp(-2*pi*1i*fourierPhaseRange) * ...
                    fluxPixRange;

            end

        end

    end

    rangeVisibilityRuntime_s = ...
        toc(ticRangeVisibility);

    if any(~isfinite(real(VRotAccessibleRange)),'all') || ...
            any(~isfinite(imag(VRotAccessibleRange)),'all')

        error('At least one range-sweep rotating visibility was not evaluated.');

    end


    %% Flatten the rotating master data in the production sample order

    switch stroboscopicSampleOrder

        case 'baseline-major'

            VRotRangeVector = ...
                VRotAccessibleRange(:);

        case 'epoch-major'

            VRotRangeVector = ...
                reshape(VRotAccessibleRange.',[],1);

        otherwise

            error('Unexpected stroboscopicSampleOrder.');

    end


    %% Exact-phase range audit

    zeroPhaseDifferenceRange = ...
        VRotRangeVector(zeroPhaseSampleMask) - ...
        VFrozenRangeVector(zeroPhaseSampleMask);

    rangeExactPhaseVisibilityError(iRange) = ...
        max(abs(zeroPhaseDifferenceRange));

    if rangeExactPhaseVisibilityError(iRange) > ...
            zeroPhaseVisibilityTolerance

        error(['Exact-phase rotating/frozen visibility mismatch at ', ...
            'range %.6e km: %.3e.'], ...
            thisRange_km, ...
            rangeExactPhaseVisibilityError(iRange));
    end


    %% Report target resolution before reconstruction

    fprintf('  projected principal span     = %.6f x %.6f mas\n', ...
        projectedMajorSpan_mas(iRange), ...
        projectedMinorSpan_mas(iRange));

    fprintf('  resolved elements            = %.3f x %.3f\n', ...
        resolvedElementsMajor(iRange), ...
        resolvedElementsMinor(iRange));

    fprintf('  angular pixel scale          = %.6e x %.6e mas\n', ...
        rangePixelScaleL_mas(iRange), ...
        rangePixelScaleM_mas(iRange));

    fprintf('  exact-phase max |Vrot-V0|    = %.3e\n', ...
        rangeExactPhaseVisibilityError(iRange));

    fprintf('  rotating-visibility runtime  = %.3f s\n', ...
        rangeVisibilityRuntime_s);


    %% Campaign x bin reconstruction loop

    for iCampaign = pilotCampaignIndices

        stopAfterFirstFailure = false;

        campaignStart_s = ...
            tQPOStart_s;

        campaignEnd_s = ...
            campaignStart_s + ...
            campaignDuration_s(iCampaign);

        for iBin = pilotBinIndices

            if stopAfterFirstFailure
                break
            end

            fullBinWidth_deg = ...
                phaseBinWidth_deg(iBin);

            halfBinWidth_deg = ...
                0.5*fullBinWidth_deg;

            inCampaignSample = ...
                sampleTime_s >= campaignStart_s - 1e-9 & ...
                sampleTime_s <= campaignEnd_s + 1e-9;

            if fullBinWidth_deg == 0

                inPhaseBinSample = ...
                    abs(samplePhaseOffset_deg) <= ...
                    phaseGridTolerance_deg;

            else

                inPhaseBinSample = ...
                    abs(samplePhaseOffset_deg) <= ...
                    halfBinWidth_deg + ...
                    phaseGridTolerance_deg;

            end

            keepSample = ...
                inCampaignSample & ...
                inPhaseBinSample;

            selectedMasterEpochs = ...
                unique(sampleMasterEpoch(keepSample),'stable');

            nSelectedEpochs = ...
                numel(selectedMasterEpochs);

            if nSelectedEpochs ~= ...
                    caseAccessibleEpochs(iCampaign,iBin)

                error(['Range sweep selected a different number of ', ...
                    'accessible epochs than Section 22 for campaign ', ...
                    '%.3f TPrimary and W = %.1f deg.'], ...
                    campaignFractions(iCampaign), ...
                    fullBinWidth_deg);
            end

            if nnz(keepSample) ~= ...
                    casePhysicalSamples(iCampaign,iBin)

                error(['Range sweep selected a different number of ', ...
                    'physical samples than Section 22.']);
            end


            %% Explicit UNASSESSED states

            if ~caseSuccess(iCampaign,iBin)
                continue
            end


            try

                %% Frozen-control observation for this exact sample mask

                obsFrozenCase = ...
                    obsFrozenRangeMaster;

                obsFrozenCase.u = ...
                    obsFrozenRangeMaster.u(keepSample);

                obsFrozenCase.v = ...
                    obsFrozenRangeMaster.v(keepSample);

                if isfield(obsFrozenRangeMaster,'w') && ...
                        numel(obsFrozenRangeMaster.w) == ...
                        numel(VFrozenRangeVector)

                    obsFrozenCase.w = ...
                        obsFrozenRangeMaster.w(keepSample);

                end

                obsFrozenCase.V = ...
                    VFrozenRangeVector(keepSample);

                obsFrozenCase.rho = ...
                    hypot(obsFrozenCase.u,obsFrozenCase.v);

                obsFrozenCase.rho_Glambda = ...
                    obsFrozenCase.rho/1e9;

                obsFrozenCase.accessibleIndex = ...
                    selectedMasterEpochs(:);

                if isfield(obsFrozenCase,'amplitude')
                    obsFrozenCase.amplitude = abs(obsFrozenCase.V);
                end

                if isfield(obsFrozenCase,'phase')
                    obsFrozenCase.phase = angle(obsFrozenCase.V);
                end

                if isfield(obsFrozenCase,'phase_rad')
                    obsFrozenCase.phase_rad = angle(obsFrozenCase.V);
                end

                if isfield(obsFrozenCase,'phase_deg')
                    obsFrozenCase.phase_deg = ...
                        rad2deg(angle(obsFrozenCase.V));
                end


                %% Rotating-target observation on IDENTICAL samples

                obsRotatingCase = ...
                    obsFrozenCase;

                obsRotatingCase.V = ...
                    VRotRangeVector(keepSample);

                if isfield(obsRotatingCase,'amplitude')
                    obsRotatingCase.amplitude = abs(obsRotatingCase.V);
                end

                if isfield(obsRotatingCase,'phase')
                    obsRotatingCase.phase = angle(obsRotatingCase.V);
                end

                if isfield(obsRotatingCase,'phase_rad')
                    obsRotatingCase.phase_rad = angle(obsRotatingCase.V);
                end

                if isfield(obsRotatingCase,'phase_deg')
                    obsRotatingCase.phase_deg = ...
                        rad2deg(angle(obsRotatingCase.V));
                end


                %% Dirty images

                dirtyFrozenRangeCase = ...
                    make_dirty_image_from_observations( ...
                        obsFrozenCase, ...
                        simulatedRange.l_rad, ...
                        simulatedRange.m_rad, ...
                        addHermitian, ...
                        weightMode);

                dirtyRotatingRangeCase = ...
                    make_dirty_image_from_observations( ...
                        obsRotatingCase, ...
                        simulatedRange.l_rad, ...
                        simulatedRange.m_rad, ...
                        addHermitian, ...
                        weightMode);


                %% Same-sampling beam audit

                dirtyBeamDifferenceRange = ...
                    max(abs( ...
                        dirtyFrozenRangeCase.beam(:) - ...
                        dirtyRotatingRangeCase.beam(:)));

                if dirtyBeamDifferenceRange > beamEqualityTolerance
                    error(['Frozen and rotating dirty beams differ under ', ...
                        'identical sampling. Max difference = %.3e.'], ...
                        dirtyBeamDifferenceRange);
                end


                %% CLEAN frozen control

                cleanParamsFrozenRange = struct();

                cleanParamsFrozenRange.loopGain = ...
                    cleanLoopGain;

                cleanParamsFrozenRange.nIterMax = ...
                    cleanNIterMax;

                cleanParamsFrozenRange.cleanBeamFitLevel = ...
                    cleanBeamFitLevel;

                cleanParamsFrozenRange.cleanBeamTruncateSigma = ...
                    cleanBeamTruncateSigma;

                dirtyBeamPeakFrozenRange = ...
                    max(real(dirtyFrozenRangeCase.beam(:)));

                if dirtyBeamPeakFrozenRange <= 0
                    error('Frozen dirty beam has a nonpositive peak.');
                end

                dirtyImageNormalizedFrozenRange = ...
                    real(dirtyFrozenRangeCase.image) / ...
                    dirtyBeamPeakFrozenRange;

                cleanParamsFrozenRange.stopThreshold = ...
                    cleanStopFraction * ...
                    max(abs(dirtyImageNormalizedFrozenRange(:)));

                cleanFrozenRangeCase = ...
                    make_clean_image_from_dirty( ...
                        dirtyFrozenRangeCase, ...
                        cleanParamsFrozenRange);


                %% CLEAN rotating target

                cleanParamsRotatingRange = ...
                    cleanParamsFrozenRange;

                dirtyBeamPeakRotatingRange = ...
                    max(real(dirtyRotatingRangeCase.beam(:)));

                if dirtyBeamPeakRotatingRange <= 0
                    error('Rotating dirty beam has a nonpositive peak.');
                end

                dirtyImageNormalizedRotatingRange = ...
                    real(dirtyRotatingRangeCase.image) / ...
                    dirtyBeamPeakRotatingRange;

                cleanParamsRotatingRange.stopThreshold = ...
                    cleanStopFraction * ...
                    max(abs(dirtyImageNormalizedRotatingRange(:)));

                cleanRotatingRangeCase = ...
                    make_clean_image_from_dirty( ...
                        dirtyRotatingRangeCase, ...
                        cleanParamsRotatingRange);


                %% Native CLEAN-beam diagnostics
                %
                % Keep the fitted native beam only as a diagnostic.  Its fit
                % changes with the angular image-grid spacing even though the
                % physical UV sampling is unchanged.

                rangeNativeBeamMajor_mas(iRange,iCampaign,iBin) = ...
                    cleanFrozenRangeCase.params.cleanBeamFwhmMajor_mas;

                rangeNativeBeamMinor_mas(iRange,iCampaign,iBin) = ...
                    cleanFrozenRangeCase.params.cleanBeamFwhmMinor_mas;

                rangeNativeBeamPA_deg(iRange,iCampaign,iBin) = ...
                    cleanFrozenRangeCase.params.cleanBeamPaDeg;

                rangeNativeBeamMajorShift_mas(iRange,iCampaign,iBin) = ...
                    rangeNativeBeamMajor_mas(iRange,iCampaign,iBin) - ...
                    caseBeamMajor_mas(iCampaign,iBin);

                rangeNativeBeamMinorShift_mas(iRange,iCampaign,iBin) = ...
                    rangeNativeBeamMinor_mas(iRange,iCampaign,iBin) - ...
                    caseBeamMinor_mas(iCampaign,iBin);

                rangeCaseSuccess(iRange,iCampaign,iBin) = true;


                %% Restore the CLEAN component models DIRECTLY with the
                %  one explicit analysis beam.
                %
                % clean.model is the Högbom delta-component model.  Restoring
                % that same model with analysisBeamKernelRange gives every
                % target range and every phase-bin case exactly the same
                % angular comparison PSF by construction.
                %
                % This is preferable to first restoring with the range-grid
                % fitted CLEAN beam and then trying to deconvolve/reconvolve
                % covariance differences.  The latter made physically
                % identical UV resolution appear "UNASSESSED" solely because
                % the Gaussian beam fit moved slightly on a different angular
                % pixel grid.
                %
                % As in the earlier morphology study, the final residual is
                % intentionally excluded from the quantitative image metric;
                % the component model has a precisely controlled restoring
                % beam.

                imageFrozenAnalysisRange = ...
                    real(conv2( ...
                        real(cleanFrozenRangeCase.model), ...
                        analysisBeamKernelRange, ...
                        'same'));

                imageRotatingAnalysisRange = ...
                    real(conv2( ...
                        real(cleanRotatingRangeCase.model), ...
                        analysisBeamKernelRange, ...
                        'same'));

                rangeCaseResolutionComparable( ...
                    iRange,iCampaign,iBin) = true;


                %% Independent global morphology amplitude fits

                frozenNorm2Range = ...
                    sum(imageFrozenAnalysisRange(:).^2);

                rotatingNorm2Range = ...
                    sum(imageRotatingAnalysisRange(:).^2);

                if frozenNorm2Range <= 0 || ...
                        rotatingNorm2Range <= 0

                    error('At least one analysis-resolution image has zero norm.');
                end

                alphaFrozenRange = ...
                    sum( ...
                        truthAnalysisRange(:) .* ...
                        imageFrozenAnalysisRange(:)) / ...
                    frozenNorm2Range;

                alphaRotatingRange = ...
                    sum( ...
                        truthAnalysisRange(:) .* ...
                        imageRotatingAnalysisRange(:)) / ...
                    rotatingNorm2Range;

                imageFrozenAnalysisRange = ...
                    alphaFrozenRange * ...
                    imageFrozenAnalysisRange;

                imageRotatingAnalysisRange = ...
                    alphaRotatingRange * ...
                    imageRotatingAnalysisRange;


                %% Support-limited NRMSE

                errorFrozenRange = ...
                    imageFrozenAnalysisRange - ...
                    truthAnalysisRange;

                errorRotatingRange = ...
                    imageRotatingAnalysisRange - ...
                    truthAnalysisRange;

                rangeNRMSEFrozen(iRange,iCampaign,iBin) = ...
                    sqrt(mean( ...
                        errorFrozenRange(supportAnalysisRange).^2)) / ...
                    peakTruthAnalysisRange;

                rangeNRMSERotating(iRange,iCampaign,iBin) = ...
                    sqrt(mean( ...
                        errorRotatingRange(supportAnalysisRange).^2)) / ...
                    peakTruthAnalysisRange;


                %% Support-limited NCC

                truthVectorRange = ...
                    truthAnalysisRange(supportAnalysisRange);

                frozenVectorRange = ...
                    imageFrozenAnalysisRange(supportAnalysisRange);

                rotatingVectorRange = ...
                    imageRotatingAnalysisRange(supportAnalysisRange);

                truthVectorRange = ...
                    truthVectorRange-mean(truthVectorRange);

                frozenVectorRange = ...
                    frozenVectorRange-mean(frozenVectorRange);

                rotatingVectorRange = ...
                    rotatingVectorRange-mean(rotatingVectorRange);

                denominatorFrozenRange = ...
                    norm(truthVectorRange) * ...
                    norm(frozenVectorRange);

                denominatorRotatingRange = ...
                    norm(truthVectorRange) * ...
                    norm(rotatingVectorRange);

                if denominatorFrozenRange > 0

                    rangeNCCFrozen(iRange,iCampaign,iBin) = ...
                        dot( ...
                            truthVectorRange, ...
                            frozenVectorRange) / ...
                        denominatorFrozenRange;

                end

                if denominatorRotatingRange > 0

                    rangeNCCRotating(iRange,iCampaign,iBin) = ...
                        dot( ...
                            truthVectorRange, ...
                            rotatingVectorRange) / ...
                        denominatorRotatingRange;

                end


                %% 50-percent contour IoU

                truthContourRange = ...
                    truthAnalysisRange >= ...
                    contourLevel*peakTruthAnalysisRange;

                peakFrozenRange = ...
                    max(imageFrozenAnalysisRange(:));

                peakRotatingRange = ...
                    max(imageRotatingAnalysisRange(:));

                if peakFrozenRange > 0

                    frozenContourRange = ...
                        imageFrozenAnalysisRange >= ...
                        contourLevel*peakFrozenRange;

                    unionFrozenRange = ...
                        nnz(truthContourRange | frozenContourRange);

                    if unionFrozenRange > 0

                        rangeIoUFrozen(iRange,iCampaign,iBin) = ...
                            nnz( ...
                                truthContourRange & ...
                                frozenContourRange) / ...
                            unionFrozenRange;

                    end

                end

                if peakRotatingRange > 0

                    rotatingContourRange = ...
                        imageRotatingAnalysisRange >= ...
                        contourLevel*peakRotatingRange;

                    unionRotatingRange = ...
                        nnz(truthContourRange | rotatingContourRange);

                    if unionRotatingRange > 0

                        rangeIoURotating(iRange,iCampaign,iBin) = ...
                            nnz( ...
                                truthContourRange & ...
                                rotatingContourRange) / ...
                            unionRotatingRange;

                    end

                end


                %% Motion-specific reconstruction penalties

                rangeDeltaNRMSEBlur(iRange,iCampaign,iBin) = ...
                    rangeNRMSERotating(iRange,iCampaign,iBin) - ...
                    rangeNRMSEFrozen(iRange,iCampaign,iBin);

                rangeDeltaNCCBlur(iRange,iCampaign,iBin) = ...
                    rangeNCCFrozen(iRange,iCampaign,iBin) - ...
                    rangeNCCRotating(iRange,iCampaign,iBin);

                rangeDeltaIoUBlur(iRange,iCampaign,iBin) = ...
                    rangeIoUFrozen(iRange,iCampaign,iBin) - ...
                    rangeIoURotating(iRange,iCampaign,iBin);


                %% PASS / FAIL classification for an assessed case

                rangeCaseSuccess(iRange,iCampaign,iBin) = true;
                rangeCaseAssessed(iRange,iCampaign,iBin) = true;

                passesSelectedCriteria = ...
                    isfinite(rangeNRMSERotating(iRange,iCampaign,iBin)) && ...
                    isfinite(rangeNCCRotating(iRange,iCampaign,iBin)) && ...
                    isfinite(rangeDeltaNRMSEBlur(iRange,iCampaign,iBin)) && ...
                    isfinite(rangeDeltaNCCBlur(iRange,iCampaign,iBin)) && ...
                    rangeNRMSERotating(iRange,iCampaign,iBin) <= ...
                        fidelitySupportNRMSEMax && ...
                    rangeNCCRotating(iRange,iCampaign,iBin) >= ...
                        fidelityNCCMin && ...
                    rangeDeltaNRMSEBlur(iRange,iCampaign,iBin) <= ...
                        maxBlurNRMSEPenalty && ...
                    rangeDeltaNCCBlur(iRange,iCampaign,iBin) <= ...
                        maxBlurNCCPenalty;

                if passesSelectedCriteria

                    rangeCaseAssessmentCode(iRange,iCampaign,iBin) = ...
                        uint8(2);

                else

                    rangeCaseAssessmentCode(iRange,iCampaign,iBin) = ...
                        uint8(1);

                    % The pilot is designed only to locate the first
                    % contiguous-from-zero failure.  Once it occurs, wider
                    % bins are not needed to bracket that threshold.
                    stopAfterFirstFailure = true;

                end


                %% Store only representative longest-campaign images

                if iCampaign == iReferenceCampaign

                    representativeLocation = ...
                        find( ...
                            representativeBinIndex == iBin, ...
                            1, ...
                            'first');

                    if ~isempty(representativeLocation)

                        representativeFrozenAnalysis{ ...
                            iRange,representativeLocation} = ...
                            imageFrozenAnalysisRange;

                        representativeRotatingAnalysis{ ...
                            iRange,representativeLocation} = ...
                            imageRotatingAnalysisRange;

                    end

                end


            catch ME

                rangeCaseSuccess(iRange,iCampaign,iBin) = false;
                rangeCaseAssessed(iRange,iCampaign,iBin) = false;
                rangeCaseAssessmentCode(iRange,iCampaign,iBin) = uint8(0);

                warning( ...
                    ['Range-sweep reconstruction failed at D = %.6e km, ', ...
                     'campaign %.3f TPrimary, W = %.1f deg:\n%s'], ...
                    thisRange_km, ...
                    campaignFractions(iCampaign), ...
                    fullBinWidth_deg, ...
                    ME.message);

            end

        end

    end

end

fprintf('\n============================================================\n');
fprintf('TARGET-RESOLUTION / RANGE SWEEP COMPLETE\n');
fprintf('============================================================\n');

validBeamShift = ...
    isfinite(rangeNativeBeamMajorShift_mas) & ...
    isfinite(rangeNativeBeamMinorShift_mas);

if any(validBeamShift,'all')

    maxMajorBeamShift_mas = ...
        max(abs(rangeNativeBeamMajorShift_mas(validBeamShift)));

    maxMinorBeamShift_mas = ...
        max(abs(rangeNativeBeamMinorShift_mas(validBeamShift)));

    fprintf('\nRange-grid CLEAN-beam fit diagnostic\n');
    fprintf('  max |major FWHM shift|       = %.6e mas\n', ...
        maxMajorBeamShift_mas);
    fprintf('  max |minor FWHM shift|       = %.6e mas\n', ...
        maxMinorBeamShift_mas);
    fprintf(['  These are sampled-grid fit shifts only. They no longer ', ...
        'control assessment because the CLEAN delta model is restored ', ...
        'directly with the fixed analysis beam.\n']);

end


%% ============================================================
%  25. TRI-STATE ASSESSMENT + BLUR-LIMIT EXTRACTION + FIGURES
%  ============================================================
%
% Exact-phase status:
%
%   PASS       : W = 0 is assessed and satisfies the selected criteria.
%   FAIL       : W = 0 is assessed but violates at least one criterion.
%   UNASSESSED : W = 0 could not be reconstructed/evaluated numerically.
%
% Finite-bin limit status:
%
%   BRACKETED
%       exact phase passes and a wider tested bin fails.  The true threshold
%       lies between the last passing width and first failed width.
%
%   FULL_ROTATION_PASSES
%       every tested bin through W = 360 deg passes.  No blur-limited
%       threshold was found within one full asteroid rotation.
%
%   EXACT_PHASE_FAIL
%       the stroboscopic strategy itself does not satisfy the selected
%       criteria at this target resolution / campaign duration.
%
%   UNASSESSED
%       exact phase could not be reconstructed/evaluated numerically.
%
%   INCOMPLETE_AFTER_PASS
%       exact phase passes, but an unassessed case occurs before a failed
%       case, so a clean contiguous threshold cannot be bracketed.

exactPhaseAssessment = ...
    strings(nRangeCases,nCampaignCases);

binLimitStatus = ...
    strings(nRangeCases,nCampaignCases);

lastPassingBinWidth_deg = ...
    nan(nRangeCases,nCampaignCases);

lastPassingBinTime_min = ...
    nan(nRangeCases,nCampaignCases);

firstFailedBinWidth_deg = ...
    nan(nRangeCases,nCampaignCases);

firstFailedBinTime_min = ...
    nan(nRangeCases,nCampaignCases);

firstUnassessedBinWidth_deg = ...
    nan(nRangeCases,nCampaignCases);

for iRange = 1:nRangeCases

    for iCampaign = pilotCampaignIndices

        exactCode = ...
            rangeCaseAssessmentCode(iRange,iCampaign,iZeroBin);

        if exactCode == 0

            exactPhaseAssessment(iRange,iCampaign) = ...
                "UNASSESSED";

            binLimitStatus(iRange,iCampaign) = ...
                "UNASSESSED";

            continue

        elseif exactCode == 1

            exactPhaseAssessment(iRange,iCampaign) = ...
                "FAIL";

            binLimitStatus(iRange,iCampaign) = ...
                "EXACT_PHASE_FAIL";

            continue

        elseif exactCode == 2

            exactPhaseAssessment(iRange,iCampaign) = ...
                "PASS";

        else

            error('Unexpected rangeCaseAssessmentCode.');

        end


        %% Exact phase passes: walk outward contiguously in bin width

        lastPassingBinWidth_deg(iRange,iCampaign) = ...
            phaseBinWidth_deg(iZeroBin);

        statusResolved = false;

        for iPilotBin = 2:nPilotBins

            iBin = ...
                pilotBinIndices(iPilotBin);

            thisCode = ...
                rangeCaseAssessmentCode(iRange,iCampaign,iBin);

            if thisCode == 2

                lastPassingBinWidth_deg(iRange,iCampaign) = ...
                    phaseBinWidth_deg(iBin);

            elseif thisCode == 1

                firstFailedBinWidth_deg(iRange,iCampaign) = ...
                    phaseBinWidth_deg(iBin);

                binLimitStatus(iRange,iCampaign) = ...
                    "BRACKETED";

                statusResolved = true;
                break

            elseif thisCode == 0

                firstUnassessedBinWidth_deg(iRange,iCampaign) = ...
                    phaseBinWidth_deg(iBin);

                binLimitStatus(iRange,iCampaign) = ...
                    "INCOMPLETE_AFTER_PASS";

                statusResolved = true;
                break

            else

                error('Unexpected rangeCaseAssessmentCode.');

            end

        end

        if ~statusResolved

            binLimitStatus(iRange,iCampaign) = ...
                "FULL_ROTATION_PASSES";

        end


        %% Convert full phase width to full observing-window duration

        lastPassingBinTime_min(iRange,iCampaign) = ...
            lastPassingBinWidth_deg(iRange,iCampaign)/360 * ...
            TItokawa_s/60;

        if isfinite(firstFailedBinWidth_deg(iRange,iCampaign))

            firstFailedBinTime_min(iRange,iCampaign) = ...
                firstFailedBinWidth_deg(iRange,iCampaign)/360 * ...
                TItokawa_s/60;

        end

    end

end


%% Human-readable per-case assessment labels

rangeCaseAssessment = ...
    strings(nRangeCases,nCampaignCases,nBinCases);

rangeCaseAssessment(rangeCaseAssessmentCode == 0) = ...
    "UNASSESSED";

rangeCaseAssessment(rangeCaseAssessmentCode == 1) = ...
    "FAIL";

rangeCaseAssessment(rangeCaseAssessmentCode == 2) = ...
    "PASS";


%% Compact target-resolution table

resolutionStudyTable = ...
    table( ...
        rangeStudy_km(:), ...
        rangeScaleStudy(:), ...
        projectedMajorSpan_mas(:), ...
        projectedMinorSpan_mas(:), ...
        resolvedElementsMajor(:), ...
        resolvedElementsMinor(:), ...
        rangePixelScaleL_mas(:), ...
        rangeExactPhaseVisibilityError(:), ...
        'VariableNames',{ ...
            'Range_km', ...
            'RangeScaleFromBase', ...
            'ProjectedMajorSpan_mas', ...
            'ProjectedMinorSpan_mas', ...
            'ResolvedElementsMajor', ...
            'ResolvedElementsMinor', ...
            'PixelScaleL_mas', ...
            'ExactPhaseVisibilityError'});


%% Campaign-level threshold / bracket table

nSummaryRows = ...
    nRangeCases*nPilotCampaigns;

summaryRange_km = ...
    zeros(nSummaryRows,1);

summaryResolvedMajor = ...
    zeros(nSummaryRows,1);

summaryResolvedMinor = ...
    zeros(nSummaryRows,1);

summaryCampaign_TPrimary = ...
    zeros(nSummaryRows,1);

summaryCampaign_days = ...
    zeros(nSummaryRows,1);

summaryExactPhase = ...
    strings(nSummaryRows,1);

summaryLimitStatus = ...
    strings(nSummaryRows,1);

summaryLastPass_deg = ...
    nan(nSummaryRows,1);

summaryLastPass_min = ...
    nan(nSummaryRows,1);

summaryFirstFail_deg = ...
    nan(nSummaryRows,1);

summaryFirstFail_min = ...
    nan(nSummaryRows,1);

summaryFirstUnassessed_deg = ...
    nan(nSummaryRows,1);

summaryRow = 0;

for iRange = 1:nRangeCases

    for iCampaign = pilotCampaignIndices

        summaryRow = summaryRow+1;

        summaryRange_km(summaryRow) = ...
            rangeStudy_km(iRange);

        summaryResolvedMajor(summaryRow) = ...
            resolvedElementsMajor(iRange);

        summaryResolvedMinor(summaryRow) = ...
            resolvedElementsMinor(iRange);

        summaryCampaign_TPrimary(summaryRow) = ...
            campaignFractions(iCampaign);

        summaryCampaign_days(summaryRow) = ...
            campaignDuration_days(iCampaign);

        summaryExactPhase(summaryRow) = ...
            exactPhaseAssessment(iRange,iCampaign);

        summaryLimitStatus(summaryRow) = ...
            binLimitStatus(iRange,iCampaign);

        summaryLastPass_deg(summaryRow) = ...
            lastPassingBinWidth_deg(iRange,iCampaign);

        summaryLastPass_min(summaryRow) = ...
            lastPassingBinTime_min(iRange,iCampaign);

        summaryFirstFail_deg(summaryRow) = ...
            firstFailedBinWidth_deg(iRange,iCampaign);

        summaryFirstFail_min(summaryRow) = ...
            firstFailedBinTime_min(iRange,iCampaign);

        summaryFirstUnassessed_deg(summaryRow) = ...
            firstUnassessedBinWidth_deg(iRange,iCampaign);

    end

end

resolutionCampaignSummary = ...
    table( ...
        summaryRange_km, ...
        summaryResolvedMajor, ...
        summaryResolvedMinor, ...
        summaryCampaign_TPrimary, ...
        summaryCampaign_days, ...
        summaryExactPhase, ...
        summaryLimitStatus, ...
        summaryLastPass_deg, ...
        summaryLastPass_min, ...
        summaryFirstFail_deg, ...
        summaryFirstFail_min, ...
        summaryFirstUnassessed_deg, ...
        'VariableNames',{ ...
            'Range_km', ...
            'ResolvedElementsMajor', ...
            'ResolvedElementsMinor', ...
            'Campaign_TPrimary', ...
            'Campaign_days', ...
            'ExactPhaseAssessment', ...
            'BinLimitStatus', ...
            'LastPassingFullBinWidth_deg', ...
            'LastPassingFullBinTime_min', ...
            'FirstFailedFullBinWidth_deg', ...
            'FirstFailedFullBinTime_min', ...
            'FirstUnassessedFullBinWidth_deg'});


%% Complete long-form range x campaign x bin table

nLongRows = ...
    nRangeCases*nPilotCampaigns*nPilotBins;

longRange_km = zeros(nLongRows,1);
longResolvedMajor = zeros(nLongRows,1);
longResolvedMinor = zeros(nLongRows,1);
longCampaign = zeros(nLongRows,1);
longCampaignDays = zeros(nLongRows,1);
longBinWidth = zeros(nLongRows,1);
longBinTime = zeros(nLongRows,1);
longAssessment = strings(nLongRows,1);
longNRMSEFrozen = nan(nLongRows,1);
longNRMSERotating = nan(nLongRows,1);
longNCCFrozen = nan(nLongRows,1);
longNCCRotating = nan(nLongRows,1);
longIoUFrozen = nan(nLongRows,1);
longIoURotating = nan(nLongRows,1);
longDeltaNRMSE = nan(nLongRows,1);
longDeltaNCC = nan(nLongRows,1);
longDeltaIoU = nan(nLongRows,1);

longRow = 0;

for iRange = 1:nRangeCases

    for iCampaign = pilotCampaignIndices

        for iBin = pilotBinIndices

            longRow = longRow+1;

            longRange_km(longRow) = ...
                rangeStudy_km(iRange);

            longResolvedMajor(longRow) = ...
                resolvedElementsMajor(iRange);

            longResolvedMinor(longRow) = ...
                resolvedElementsMinor(iRange);

            longCampaign(longRow) = ...
                campaignFractions(iCampaign);

            longCampaignDays(longRow) = ...
                campaignDuration_days(iCampaign);

            longBinWidth(longRow) = ...
                phaseBinWidth_deg(iBin);

            longBinTime(longRow) = ...
                phaseBinWidth_deg(iBin)/360 * ...
                TItokawa_s/60;

            longAssessment(longRow) = ...
                rangeCaseAssessment(iRange,iCampaign,iBin);

            longNRMSEFrozen(longRow) = ...
                rangeNRMSEFrozen(iRange,iCampaign,iBin);

            longNRMSERotating(longRow) = ...
                rangeNRMSERotating(iRange,iCampaign,iBin);

            longNCCFrozen(longRow) = ...
                rangeNCCFrozen(iRange,iCampaign,iBin);

            longNCCRotating(longRow) = ...
                rangeNCCRotating(iRange,iCampaign,iBin);

            longIoUFrozen(longRow) = ...
                rangeIoUFrozen(iRange,iCampaign,iBin);

            longIoURotating(longRow) = ...
                rangeIoURotating(iRange,iCampaign,iBin);

            longDeltaNRMSE(longRow) = ...
                rangeDeltaNRMSEBlur(iRange,iCampaign,iBin);

            longDeltaNCC(longRow) = ...
                rangeDeltaNCCBlur(iRange,iCampaign,iBin);

            longDeltaIoU(longRow) = ...
                rangeDeltaIoUBlur(iRange,iCampaign,iBin);

        end

    end

end

resolutionSweepCaseTable = ...
    table( ...
        longRange_km, ...
        longResolvedMajor, ...
        longResolvedMinor, ...
        longCampaign, ...
        longCampaignDays, ...
        longBinWidth, ...
        longBinTime, ...
        longAssessment, ...
        longNRMSEFrozen, ...
        longNRMSERotating, ...
        longNCCFrozen, ...
        longNCCRotating, ...
        longIoUFrozen, ...
        longIoURotating, ...
        longDeltaNRMSE, ...
        longDeltaNCC, ...
        longDeltaIoU, ...
        'VariableNames',{ ...
            'Range_km', ...
            'ResolvedElementsMajor', ...
            'ResolvedElementsMinor', ...
            'Campaign_TPrimary', ...
            'Campaign_days', ...
            'FullBinWidth_deg', ...
            'FullBinTime_min', ...
            'Assessment', ...
            'FrozenSupportNRMSE', ...
            'RotatingSupportNRMSE', ...
            'FrozenNCC', ...
            'RotatingNCC', ...
            'FrozenIoU', ...
            'RotatingIoU', ...
            'DeltaNRMSE_Blur', ...
            'DeltaNCC_Blur', ...
            'DeltaIoU_Blur'});


%% Numerical summary

fprintf('\n============================================================\n');
fprintf('RESOLUTION-DEPENDENT STROBOSCOPIC RESULT\n');
fprintf('============================================================\n');
fprintf('Fixed analysis FWHM          = %.6f mas\n',analysisBeamFwhm_mas);
fprintf('Support NRMSE threshold      = %.6f\n',fidelitySupportNRMSEMax);
fprintf('NCC threshold                = %.6f\n',fidelityNCCMin);
fprintf('Max additional NRMSE penalty = %.6f\n',maxBlurNRMSEPenalty);
fprintf('Max additional NCC penalty   = %.6f\n',maxBlurNCCPenalty);
fprintf('------------------------------------------------------------\n');

for iRange = 1:nRangeCases

    fprintf('\nRange = %.6e km\n',rangeStudy_km(iRange));
    fprintf('  projected span              = %.4f x %.4f mas\n', ...
        projectedMajorSpan_mas(iRange), ...
        projectedMinorSpan_mas(iRange));
    fprintf('  resolved elements           = %.3f x %.3f\n', ...
        resolvedElementsMajor(iRange), ...
        resolvedElementsMinor(iRange));

    for iCampaign = pilotCampaignIndices

        fprintf('  campaign %.3f TPrimary (%.3f d): exact = %s, limit = %s', ...
            campaignFractions(iCampaign), ...
            campaignDuration_days(iCampaign), ...
            char(exactPhaseAssessment(iRange,iCampaign)), ...
            char(binLimitStatus(iRange,iCampaign)));

        if binLimitStatus(iRange,iCampaign) == "BRACKETED"

            fprintf(', W in (%.1f, %.1f) deg', ...
                lastPassingBinWidth_deg(iRange,iCampaign), ...
                firstFailedBinWidth_deg(iRange,iCampaign));

            fprintf(', time in (%.2f, %.2f) min', ...
                lastPassingBinTime_min(iRange,iCampaign), ...
                firstFailedBinTime_min(iRange,iCampaign));

        elseif binLimitStatus(iRange,iCampaign) == ...
                "FULL_ROTATION_PASSES"

            fprintf(', all tested widths through 360 deg pass');

        elseif binLimitStatus(iRange,iCampaign) == ...
                "INCOMPLETE_AFTER_PASS"

            fprintf(', first unassessed width = %.1f deg', ...
                firstUnassessedBinWidth_deg(iRange,iCampaign));

        end

        fprintf('\n');

    end

end

fprintf('============================================================\n');

fprintf('\nTARGET RESOLUTION TABLE\n');
disp(resolutionStudyTable);

fprintf('\nRANGE / CAMPAIGN SUMMARY\n');
disp(resolutionCampaignSummary);


%% Figure 1: number of resolved target elements versus range

figResolvedElements = figure( ...
    'Name','Resolved Itokawa elements versus range', ...
    'Color','w', ...
    'Position',[120 120 900 650]);

axResolvedElements = axes(figResolvedElements);

semilogx( ...
    axResolvedElements, ...
    rangeStudy_km, ...
    resolvedElementsMajor, ...
    'o-', ...
    'LineWidth',1.6, ...
    'MarkerSize',7, ...
    'DisplayName','Principal major span');

hold(axResolvedElements,'on');

semilogx( ...
    axResolvedElements, ...
    rangeStudy_km, ...
    resolvedElementsMinor, ...
    's-', ...
    'LineWidth',1.6, ...
    'MarkerSize',7, ...
    'DisplayName','Principal minor span');

grid(axResolvedElements,'on');
box(axResolvedElements,'on');

set(axResolvedElements,'XDir','reverse');

xlabel(axResolvedElements,'Target range [km]');
ylabel(axResolvedElements,'Resolved elements across target');

legend(axResolvedElements,'Location','best');

title( ...
    axResolvedElements, ...
    sprintf('Target resolution at fixed %.4f mas analysis FWHM', ...
        analysisBeamFwhm_mas), ...
    'Interpreter','none');


%% Figure 2: exact-phase proof of concept as target resolution increases

figExactPhaseResolution = figure( ...
    'Name','Exact-phase stroboscopic fidelity versus target resolution', ...
    'Color','w', ...
    'Position',[100 100 1200 800]);

 tlExactPhaseResolution = tiledlayout( ...
    figExactPhaseResolution, ...
    2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

campaignColors = ...
    lines(nPilotCampaigns);


%% Exact-phase NRMSE

ax = nexttile(tlExactPhaseResolution);
hold(ax,'on');
grid(ax,'on');
box(ax,'on');

for iPilotCampaign = 1:nPilotCampaigns

    iCampaign = ...
        pilotCampaignIndices(iPilotCampaign);

    yData = ...
        squeeze(rangeNRMSERotating(:,iCampaign,iZeroBin));

    assessed = ...
        squeeze(rangeCaseAssessed(:,iCampaign,iZeroBin));

    plot( ...
        ax, ...
        resolvedElementsMajor(assessed), ...
        yData(assessed), ...
        'o-', ...
        'Color',campaignColors(iPilotCampaign,:), ...
        'LineWidth',1.5, ...
        'MarkerSize',6, ...
        'DisplayName',sprintf('%.2f $T_p$',campaignFractions(iCampaign)));

end

yline(ax,fidelitySupportNRMSEMax,'k--','Study threshold', ...
    'HandleVisibility','off');

xlabel(ax,'Resolved elements across principal major span');
ylabel(ax,'Exact-phase support NRMSE');

title(ax,'Does exact-phase stroboscopic imaging remain faithful?','Interpreter','none');

legend(ax,'Location','best','Interpreter','latex');


%% Exact-phase NCC

ax = nexttile(tlExactPhaseResolution);
hold(ax,'on');
grid(ax,'on');
box(ax,'on');

for iPilotCampaign = 1:nPilotCampaigns

    iCampaign = ...
        pilotCampaignIndices(iPilotCampaign);

    yData = ...
        squeeze(rangeNCCRotating(:,iCampaign,iZeroBin));

    assessed = ...
        squeeze(rangeCaseAssessed(:,iCampaign,iZeroBin));

    plot( ...
        ax, ...
        resolvedElementsMajor(assessed), ...
        yData(assessed), ...
        'o-', ...
        'Color',campaignColors(iPilotCampaign,:), ...
        'LineWidth',1.5, ...
        'MarkerSize',6, ...
        'DisplayName',sprintf('%.2f $T_p$',campaignFractions(iCampaign)));

end

yline(ax,fidelityNCCMin,'k--','Study threshold', ...
    'HandleVisibility','off');

xlabel(ax,'Resolved elements across principal major span');
ylabel(ax,'Exact-phase NCC');

title(ax,'Exact repeated phase: morphology correlation','Interpreter','none');

legend(ax,'Location','best','Interpreter','latex');


%% Figure 3: motion penalty versus phase-bin width for the longest campaign

figMotionPenaltyResolution = figure( ...
    'Name','Fast pilot: motion penalty versus bin width and target resolution', ...
    'Color','w', ...
    'Position',[100 100 1250 800]);

 tlMotionPenaltyResolution = tiledlayout( ...
    figMotionPenaltyResolution, ...
    2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

rangeColors = ...
    lines(nRangeCases);


%% NRMSE motion penalty

ax = nexttile(tlMotionPenaltyResolution);
hold(ax,'on');
grid(ax,'on');
box(ax,'on');

for iRange = 1:nRangeCases

    yData = ...
        squeeze(rangeDeltaNRMSEBlur(iRange,iReferenceCampaign,:));

    assessed = ...
        squeeze(rangeCaseAssessed(iRange,iReferenceCampaign,:));

    plot( ...
        ax, ...
        phaseBinWidth_deg(assessed), ...
        yData(assessed), ...
        'o-', ...
        'Color',rangeColors(iRange,:), ...
        'LineWidth',1.5, ...
        'MarkerSize',6, ...
        'DisplayName',sprintf('$N_{res}=%.2f$',resolvedElementsMajor(iRange)));

end

yline(ax,maxBlurNRMSEPenalty,'k--','Study threshold', ...
    'HandleVisibility','off');

xlabel(ax,'Full phase-bin width [deg]');
ylabel(ax,'NRMSE$_{rot}$ - NRMSE$_{frozen}$');

title(ax,'Additional morphology error caused by rotation inside the bin');

legend(ax,'Location','best','Interpreter','latex');


%% NCC motion penalty

ax = nexttile(tlMotionPenaltyResolution);
hold(ax,'on');
grid(ax,'on');
box(ax,'on');

for iRange = 1:nRangeCases

    yData = ...
        squeeze(rangeDeltaNCCBlur(iRange,iReferenceCampaign,:));

    assessed = ...
        squeeze(rangeCaseAssessed(iRange,iReferenceCampaign,:));

    plot( ...
        ax, ...
        phaseBinWidth_deg(assessed), ...
        yData(assessed), ...
        'o-', ...
        'Color',rangeColors(iRange,:), ...
        'LineWidth',1.5, ...
        'MarkerSize',6, ...
        'DisplayName',sprintf('$N_{res}=%.2f$',resolvedElementsMajor(iRange)));

end

yline(ax,maxBlurNCCPenalty,'k--','Study threshold', ...
    'HandleVisibility','off');

xlabel(ax,'Full phase-bin width [deg]');
ylabel(ax,'NCC$_{frozen}$ - NCC$_{rot}$');

title(ax,'Additional morphology-correlation loss caused by binning');

legend(ax,'Location','best','Interpreter','latex');


%% Figure 4: blur-limit bracket versus resolved target size
%
% A solid circle is the last passing tested width.
% If the next wider tested bin fails, a vertical segment reaches the first
% failed time, explicitly showing the current numerical bracket.
% An upward triangle at one full rotation means all tested widths pass.

figBinLimitResolution = figure( ...
    'Name','Stroboscopic bin limit versus target resolution', ...
    'Color','w', ...
    'Position',[100 100 1050 720]);

axBinLimitResolution = axes(figBinLimitResolution);

hold(axBinLimitResolution,'on');
grid(axBinLimitResolution,'on');
box(axBinLimitResolution,'on');

for iPilotCampaign = 1:nPilotCampaigns

    iCampaign = ...
        pilotCampaignIndices(iPilotCampaign);

    for iRange = 1:nRangeCases

        xValue = ...
            resolvedElementsMajor(iRange);

        status = ...
            binLimitStatus(iRange,iCampaign);

        if status == "BRACKETED"

            yLow = ...
                lastPassingBinTime_min(iRange,iCampaign);

            yHigh = ...
                firstFailedBinTime_min(iRange,iCampaign);

            plot( ...
                axBinLimitResolution, ...
                [xValue xValue], ...
                [yLow yHigh], ...
                '-', ...
                'Color',campaignColors(iPilotCampaign,:), ...
                'LineWidth',1.1, ...
                'HandleVisibility','off');

            plot( ...
                axBinLimitResolution, ...
                xValue, ...
                yLow, ...
                'o', ...
                'Color',campaignColors(iPilotCampaign,:), ...
                'MarkerFaceColor',campaignColors(iPilotCampaign,:), ...
                'MarkerSize',7, ...
                'HandleVisibility','off');

            plot( ...
                axBinLimitResolution, ...
                xValue, ...
                yHigh, ...
                'x', ...
                'Color',campaignColors(iPilotCampaign,:), ...
                'MarkerSize',8, ...
                'LineWidth',1.4, ...
                'HandleVisibility','off');

        elseif status == "FULL_ROTATION_PASSES"

            plot( ...
                axBinLimitResolution, ...
                xValue, ...
                lastPassingBinTime_min(iRange,iCampaign), ...
                '^', ...
                'Color',campaignColors(iPilotCampaign,:), ...
                'MarkerFaceColor','none', ...
                'MarkerSize',8, ...
                'LineWidth',1.4, ...
                'HandleVisibility','off');

        end

    end

    % Dummy handle for one compact legend entry per campaign.
    plot( ...
        axBinLimitResolution, ...
        nan, ...
        nan, ...
        'o-', ...
        'Color',campaignColors(iPilotCampaign,:), ...
        'LineWidth',1.4, ...
        'MarkerSize',6, ...
        'DisplayName',sprintf('%.2f $T_p$',campaignFractions(iCampaign)));

end

xlabel(axBinLimitResolution,'Resolved elements across principal major span');
ylabel(axBinLimitResolution,'Full phase-bin time [min]');

title( ...
    axBinLimitResolution, ...
    'Fast-pilot stroboscopic bin-limit bracket versus target resolution', ...
    'Interpreter','none');

legend( ...
    axBinLimitResolution, ...
    'Location','best', ...
    'Interpreter','latex');


%% Figure 5: most-resolved target, longest campaign, representative bins

[~,iMostResolvedRange] = ...
    max(resolvedElementsMajor);

availableRepresentative = ...
    false(1,nRepresentativeBins);

for iRep = 1:nRepresentativeBins

    availableRepresentative(iRep) = ...
        ~isempty( ...
            representativeRotatingAnalysis{ ...
                iMostResolvedRange,iRep});

end

if any(availableRepresentative)

    selectedRepresentative = ...
        find(availableRepresentative);

    nRepresentativePlot = ...
        numel(selectedRepresentative);

    figMostResolved = figure( ...
        'Name','Most-resolved Itokawa finite-bin comparison', ...
        'Color','w', ...
        'Position',[60 60 1450 310*nRepresentativePlot]);

    tlMostResolved = tiledlayout( ...
        figMostResolved, ...
        nRepresentativePlot, ...
        3, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    truthMostResolved = ...
        rangeTruthAnalysis{iMostResolvedRange};

    truthShow = ...
        truthMostResolved/max(truthMostResolved(:));

    simulatedMostResolved = ...
        make_simulated_target_from_raster( ...
            imgFrozenRot, ...
            rangeStudy_km(iMostResolvedRange), ...
            shapeItokawa.type);

    for iPlot = 1:nRepresentativePlot

        iRep = ...
            selectedRepresentative(iPlot);

        iBin = ...
            representativeBinIndex(iRep);

        frozenImage = ...
            representativeFrozenAnalysis{ ...
                iMostResolvedRange,iRep};

        rotatingImage = ...
            representativeRotatingAnalysis{ ...
                iMostResolvedRange,iRep};

        frozenShow = ...
            frozenImage/max(frozenImage(:));

        rotatingShow = ...
            rotatingImage/max(rotatingImage(:));


        %% Truth

        ax = nexttile(tlMostResolved);

        imagesc( ...
            ax, ...
            simulatedMostResolved.l_mas, ...
            simulatedMostResolved.m_mas, ...
            truthShow);

        axis(ax,'image');
        set(ax,'YDir','normal');
        colormap(ax,hot(256));
        clim(ax,[0 1]);

        xlabel(ax,'$l$ [mas]');
        ylabel(ax,'$m$ [mas]');

        title( ...
            ax, ...
            sprintf('Truth; W = %.0f deg',phaseBinWidth_deg(iBin)), ...
            'Interpreter','none');


        %% Frozen same-sampling control

        ax = nexttile(tlMostResolved);

        imagesc( ...
            ax, ...
            simulatedMostResolved.l_mas, ...
            simulatedMostResolved.m_mas, ...
            frozenShow);

        axis(ax,'image');
        set(ax,'YDir','normal');
        colormap(ax,hot(256));
        clim(ax,[0 1]);

        xlabel(ax,'$l$ [mas]');
        ylabel(ax,'$m$ [mas]');

        title( ...
            ax, ...
            sprintf('Frozen: NRMSE %.3f, NCC %.3f', ...
                rangeNRMSEFrozen( ...
                    iMostResolvedRange, ...
                    iReferenceCampaign, ...
                    iBin), ...
                rangeNCCFrozen( ...
                    iMostResolvedRange, ...
                    iReferenceCampaign, ...
                    iBin)), ...
            'Interpreter','none');


        %% Rotating target

        ax = nexttile(tlMostResolved);

        imagesc( ...
            ax, ...
            simulatedMostResolved.l_mas, ...
            simulatedMostResolved.m_mas, ...
            rotatingShow);

        axis(ax,'image');
        set(ax,'YDir','normal');
        colormap(ax,hot(256));
        clim(ax,[0 1]);

        xlabel(ax,'$l$ [mas]');
        ylabel(ax,'$m$ [mas]');

        title( ...
            ax, ...
            sprintf('Rotating: NRMSE %.3f, NCC %.3f', ...
                rangeNRMSERotating( ...
                    iMostResolvedRange, ...
                    iReferenceCampaign, ...
                    iBin), ...
                rangeNCCRotating( ...
                    iMostResolvedRange, ...
                    iReferenceCampaign, ...
                    iBin)), ...
            'Interpreter','none');

    end

    title( ...
        tlMostResolved, ...
        sprintf([ ...
            'Most-resolved range: D = %.3e km, ', ...
            'Nres = %.2f x %.2f, campaign = %.3f days'], ...
            rangeStudy_km(iMostResolvedRange), ...
            resolvedElementsMajor(iMostResolvedRange), ...
            resolvedElementsMinor(iMostResolvedRange), ...
            campaignDuration_days(iReferenceCampaign)), ...
        'Interpreter','none');

end


%% Final interpretation flags

fprintf('\nInterpretation status by target resolution\n');
fprintf('------------------------------------------\n');

for iRange = 1:nRangeCases

    fprintf('Nres major/minor = %.3f / %.3f at D = %.3e km\n', ...
        resolvedElementsMajor(iRange), ...
        resolvedElementsMinor(iRange), ...
        rangeStudy_km(iRange));

    for iCampaign = pilotCampaignIndices

        fprintf('  %.2f Tp: exact %s; bin limit %s\n', ...
            campaignFractions(iCampaign), ...
            char(exactPhaseAssessment(iRange,iCampaign)), ...
            char(binLimitStatus(iRange,iCampaign)));

    end

end

fprintf('\n');
fprintf(['A BRACKETED case is the useful blur-limit result: the true ', ...
    'maximum faithful full bin width lies between the last passing and ', ...
    'first failed tested widths.\n']);

fprintf(['FULL_ROTATION_PASSES means no reconstruction-defined blur limit ', ...
    'was found within one complete asteroid rotation at that resolved ', ...
    'target size.\n']);

fprintf(['UNASSESSED is not a failure. It means the case could not be ', ...
    'reconstructed/evaluated numerically.\n']);

fprintf(['This file is the fast pilot: only the longest campaign and the ', ...
    'reduced phase-bin grid were reconstructed in the range sweep.\n']);
