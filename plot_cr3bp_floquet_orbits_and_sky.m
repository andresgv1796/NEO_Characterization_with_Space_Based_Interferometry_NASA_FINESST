function [fig,overlays] = plot_cr3bp_floquet_orbits_and_sky( ...
    systemName,poFile,parentRows,epochUTC,opts)
%PLOT_CR3BP_FLOQUET_ORBITS_AND_SKY
%
% With the same core inputs used by make_cr3bp_floquet_normal_overlays,
% produce a linked two-panel diagnostic:
%
%   LEFT  : selected parent periodic orbit(s) in synodic CR3BP 3-D space
%   RIGHT : corresponding ideal inertial Floquet-normal pattern(s) on the
%           celestial unit sphere
%
% Each parent uses exactly the same color in both panels and in the overlay
% structures returned for plot_JPL_asteroid_sky.
%
% EXAMPLE
% -------
%   epoch = datetime('now','TimeZone','UTC');
%
%   [fig,ov] = plot_cr3bp_floquet_orbits_and_sky( ...
%       "EarthMoon", ...
%       "L1_northern_halo_orbits_JPL_IC_Earth_Moon.csv", ...
%       [82 125 197], ...
%       epoch, ...
%       CenterMode="vertical", ...
%       SkyFrame="ecliptic");
%
%   plot_JPL_asteroid_sky( ...
%       "geocentric",epoch, ...
%       CoordinateFrame="ecliptic", ...
%       NormalOverlay=ov);
%
% Sun-Earth is identical in use:
%
%   [fig,ov] = plot_cr3bp_floquet_orbits_and_sky( ...
%       "SunEarth",sunEarthFamilyFile,[i1 i2 i3],epoch);
%

arguments
    systemName (1,1) string = "EarthMoon"
    poFile (1,1) string = "L1_northern_halo_orbits_JPL_IC_Earth_Moon.csv"
    parentRows (1,:) double {mustBeInteger,mustBePositive} = 116
    epochUTC (1,1) datetime = datetime('now','TimeZone','UTC')

    opts.CenterMode string = "vertical"
    opts.NPrimaryPeriods (1,1) double {mustBeInteger,mustBePositive} = 1
    opts.SamplesPerPeriod (1,1) double {mustBeInteger,mustBePositive} = 600
    opts.FloquetPhaseSamples (1,1) double {mustBeInteger,mustBePositive} = 61
    opts.OrbitSamples (1,1) double {mustBeInteger,mustBePositive} = 800
    opts.NormalReferenceSynodic (1,3) double = [1 0 0]

    opts.EmbeddingMode (1,1) string ...
        {mustBeMember(opts.EmbeddingMode, ...
        ["ideal_epoch_anchor","identity","custom"])} ...
        = "ideal_epoch_anchor"

    opts.PrimaryHorizonsID (1,1) string = ""
    opts.SecondaryHorizonsID (1,1) string = ""
    opts.InitialSynodicToEcliptic double = []

    opts.MassRatioSource (1,1) string ...
        {mustBeMember(opts.MassRatioSource,["catalog","system","manual"])} ...
        = "catalog"
    opts.MassRatio (1,1) double = NaN
    opts.CheckOrbitClosure (1,1) logical = true
    opts.OrbitClosureTolerance (1,1) double {mustBePositive} = 1e-8

    opts.ShowAntipode (1,1) logical = true
    opts.LineWidth (1,1) double {mustBePositive} = 2.4
    opts.Colors double = []
    opts.ColorMode (1,1) string ...
        {mustBeMember(opts.ColorMode,["lines","gradient"])} = "lines"
    opts.ColorMapName (1,1) string = "turbo"
    opts.ColorParameter (1,1) string ...
        {mustBeMember(opts.ColorParameter,["sequence","parentRow"])} = "sequence"
    opts.Labels string = strings(0,1)
    opts.Verbose (1,1) logical = true

    % Celestial-sphere display frame only; the stored overlay vectors remain
    % J2000 ecliptic so they are directly compatible with the asteroid map.
    opts.SkyFrame (1,1) string ...
        {mustBeMember(opts.SkyFrame,["ecliptic","equatorial","galactic"])} ...
        = "ecliptic"

    opts.ShowSphereGrid (1,1) logical = true
    opts.ShowStartMarkers (1,1) logical = true

    % Reuse precomputed overlays to avoid rebuilding the same Floquet tracks
    % (and re-querying Horizons) when the asteroid map already has them.
    opts.Overlays = struct([])

    % Body display. White-background plotting is now fixed by design.
    % ShowPrimaryBody=[] selects a sensible system-dependent default:
    % Earth-Moon -> true, Sun-Earth -> false (the Sun otherwise dominates
    % the scale visually). Set explicitly to true/false to override.
    opts.ShowPrimaryBody = []
    opts.ShowSecondaryBody (1,1) logical = true
    opts.BodyResolution (1,1) double {mustBeInteger,mustBePositive} = 80
    opts.BodyGST_deg (1,1) double = 0
    opts.PrimaryBodyScaleFactor (1,1) double = NaN
    opts.SecondaryBodyScaleFactor (1,1) double = NaN

    opts.FigurePosition (1,4) double = [80 80 1550 720]
end

%% Build all parent geometries once, unless the caller already did it.
if isempty(opts.Overlays)
    overlays = make_cr3bp_floquet_normal_overlays( ...
    systemName,poFile,parentRows,epochUTC, ...
    CenterMode=opts.CenterMode, ...
    NPrimaryPeriods=opts.NPrimaryPeriods, ...
    SamplesPerPeriod=opts.SamplesPerPeriod, ...
    FloquetPhaseSamples=opts.FloquetPhaseSamples, ...
    OrbitSamples=opts.OrbitSamples, ...
    NormalReferenceSynodic=opts.NormalReferenceSynodic, ...
    EmbeddingMode=opts.EmbeddingMode, ...
    PrimaryHorizonsID=opts.PrimaryHorizonsID, ...
    SecondaryHorizonsID=opts.SecondaryHorizonsID, ...
    InitialSynodicToEcliptic=opts.InitialSynodicToEcliptic, ...
    MassRatioSource=opts.MassRatioSource, ...
    MassRatio=opts.MassRatio, ...
    CheckOrbitClosure=opts.CheckOrbitClosure, ...
    OrbitClosureTolerance=opts.OrbitClosureTolerance, ...
    ShowAntipode=opts.ShowAntipode, ...
    LineWidth=opts.LineWidth, ...
    Colors=opts.Colors, ...
    ColorMode=opts.ColorMode, ...
    ColorMapName=opts.ColorMapName, ...
    ColorParameter=opts.ColorParameter, ...
    Labels=opts.Labels, ...
    Verbose=opts.Verbose);
else
    overlays = opts.Overlays;
    if ~isstruct(overlays) || isempty(overlays)
        error('Overlays must be a nonempty struct array.');
    end
end

sys = overlays(1).system;

%% Figure
fig = figure( ...
    'Theme','light', ...
    'Color','w', ...
    'Name','CR3BP parent orbits and Floquet normal sky patterns', ...
    'Position',opts.FigurePosition);

tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

%% ------------------------------------------------------------------------
% LEFT: periodic orbits in synodic coordinates
% -------------------------------------------------------------------------
axOrbit = nexttile(tl,1);
axOrbit.Tag = 'CR3BP_OrbitAxes';
hold(axOrbit,'on');
grid(axOrbit,'on');
axis(axOrbit,'equal');
box(axOrbit,'on');
axOrbit.TickLabelInterpreter = 'latex';

hOrbit = gobjects(numel(overlays),1);

for k = 1:numel(overlays)
    r = overlays(k).orbitSynodic;
    c = overlays(k).color;

    hOrbit(k) = plot3( ...
        axOrbit,r(:,1),r(:,2),r(:,3), ...
        'Color',c, ...
        'LineWidth',opts.LineWidth, ...
        'Tag',sprintf('CR3BP_Orbit_%d',k), ...
        'DisplayName',char(latex_text(overlays(k).label)));

    if opts.ShowStartMarkers && ~isempty(r)
        plot3(axOrbit,r(1,1),r(1,2),r(1,3),'o', ...
            'MarkerSize',6, ...
            'MarkerFaceColor',c, ...
            'MarkerEdgeColor','k', ...
            'Tag',sprintf('CR3BP_OrbitStart_%d',k), ...
            'HandleVisibility','off');
    end
end

[pName,sName] = body_names_from_overlays(overlays);
showPrimaryBody = resolve_show_primary_body(systemName,opts.ShowPrimaryBody);
render_cr3bp_bodies( ...
    axOrbit,systemName,sys,pName,sName,showPrimaryBody,opts.ShowSecondaryBody,opts);

xlabel(axOrbit,'$x\;[\mathrm{LU}]$','Interpreter','latex');
ylabel(axOrbit,'$y\;[\mathrm{LU}]$','Interpreter','latex');
zlabel(axOrbit,'$z\;[\mathrm{LU}]$','Interpreter','latex');
title(axOrbit,latex_text(system_display_name(sys,systemName)) + " parent periodic orbits", ...
    'Interpreter','latex');
view(axOrbit,38,25);
if numel(overlays) <= 5
    lgdOrbit = legend(axOrbit,hOrbit,'Location','best','Interpreter','latex');
else
    lgdOrbit = gobjects(0);
end

%% ------------------------------------------------------------------------
% RIGHT: celestial unit sphere normal patterns
% -------------------------------------------------------------------------
axSky = nexttile(tl,2);
axSky.Tag = 'CR3BP_SkyAxes';
hold(axSky,'on');
axis(axSky,'equal');
box(axSky,'on');
axSky.TickLabelInterpreter = 'latex';

[xs,ys,zs] = sphere(60);
surf(axSky,xs,ys,zs, ...
    'FaceColor',[0.84 0.84 0.84], ...
    'FaceAlpha',0.07, ...
    'EdgeColor',[0.70 0.70 0.70], ...
    'EdgeAlpha',0.13, ...
    'HandleVisibility','off');

if opts.ShowSphereGrid
    plot_sphere_reference_grid(axSky);
end

hSky = gobjects(numel(overlays),1);

for k = 1:numel(overlays)
    n = transform_ecliptic_vectors( ...
        overlays(k).vectorsEclipticJ2000,opts.SkyFrame);
    c = overlays(k).color;

    hSky(k) = plot3( ...
        axSky,n(:,1),n(:,2),n(:,3), ...
        'Color',c, ...
        'LineWidth',opts.LineWidth, ...
        'Tag',sprintf('CR3BP_Sky_%d',k), ...
        'DisplayName',char(latex_text(overlays(k).label)));

    if opts.ShowAntipode
        plot3(axSky,-n(:,1),-n(:,2),-n(:,3), ...
            '--','Color',c, ...
            'LineWidth',0.8*opts.LineWidth, ...
            'Tag',sprintf('CR3BP_SkyAntipode_%d',k), ...
            'HandleVisibility','off');
    end

    if opts.ShowStartMarkers && ~isempty(n)
        plot3(axSky,n(1,1),n(1,2),n(1,3),'o', ...
            'MarkerSize',6, ...
            'MarkerFaceColor',c, ...
            'MarkerEdgeColor','k', ...
            'Tag',sprintf('CR3BP_SkyStart_%d',k), ...
            'HandleVisibility','off');
    end
end

xlabel(axSky,sphere_axis_label(opts.SkyFrame,'x'),'Interpreter','latex');
ylabel(axSky,sphere_axis_label(opts.SkyFrame,'y'),'Interpreter','latex');
zlabel(axSky,sphere_axis_label(opts.SkyFrame,'z'),'Interpreter','latex');
xlim(axSky,[-1.05 1.05]);
ylim(axSky,[-1.05 1.05]);
zlim(axSky,[-1.05 1.05]);
view(axSky,38,24);

title(axSky,"Linear Floquet normal patterns: " + latex_text(opts.SkyFrame) + " frame", ...
    'Interpreter','latex');
if numel(overlays) <= 5
    lgdSky = legend(axSky,hSky,'Location','best','Interpreter','latex');
else
    lgdSky = gobjects(0);
end

%% Overall title
[~,familyName,~] = fileparts(char(poFile));
hSuper = title(tl,latex_text(systemName) + " $\mid$ " + latex_text(familyName) + ...
    " $\mid$ " + "$N_{\mathrm{orbits}}=" + string(numel(overlays)) + "$", ...
    'Interpreter','latex','FontWeight','bold');

apply_white_plot_theme(fig,[axOrbit axSky],[lgdOrbit lgdSky],hSuper);

end


%% =========================================================================
% CR3BP body rendering through the general plot_body_rotating helper
% =========================================================================

function render_cr3bp_bodies(ax,systemName,sys,pName,sName,showPrimaryBody,showSecondaryBody,opts)
mu = sys.mu;
rPrimary = [-mu 0 0];
rSecondary = [1-mu 0 0];
[primaryScale,secondaryScale] = resolve_body_scale_factors(systemName,opts);

usedPrimary = ~showPrimaryBody;
usedSecondary = ~showSecondaryBody;

if showPrimaryBody
    try
        h = plot_body_rotating(ax,pName,rPrimary, ...
            RotationDeg=opts.BodyGST_deg, ...
            ScaleFactor=primaryScale, ...
            LengthUnit_km=sys.Lstar_km, ...
            Resolution=opts.BodyResolution);
        hide_from_legend(h);
        usedPrimary = true;
    catch ME
        warning('General rendering of %s failed: %s',char(pName),ME.message);
    end
end

if showSecondaryBody
    try
        h = plot_body_rotating(ax,sName,rSecondary, ...
            RotationDeg=opts.BodyGST_deg, ...
            ScaleFactor=secondaryScale, ...
            LengthUnit_km=sys.Lstar_km, ...
            Resolution=opts.BodyResolution);
        hide_from_legend(h);
        usedSecondary = true;
    catch ME
        warning('General rendering of %s failed: %s',char(sName),ME.message);
    end
end

if showPrimaryBody && ~usedPrimary
    scatter3(ax,rPrimary(1),rPrimary(2),rPrimary(3),90,[0.20 0.20 0.20], ...
        'filled','MarkerEdgeColor','k','HandleVisibility','off');
end
if showSecondaryBody && ~usedSecondary
    scatter3(ax,rSecondary(1),rSecondary(2),rSecondary(3),60,[0.70 0.70 0.70], ...
        'filled','MarkerEdgeColor','k','HandleVisibility','off');
end
end


function tf = resolve_show_primary_body(systemName,userValue)
if isempty(userValue)
    key = lower(regexprep(char(systemName),'[^a-zA-Z0-9]',''));
    tf = ~ismember(key,{'sunearth','sunearthsystem'});
else
    if ~(islogical(userValue) && isscalar(userValue)) && ...
            ~(isnumeric(userValue) && isscalar(userValue) && ismember(userValue,[0 1]))
        error('ShowPrimaryBody must be empty, true, or false.');
    end
    tf = logical(userValue);
end
end


function [primaryScale,secondaryScale] = resolve_body_scale_factors(systemName,opts)
if isfinite(opts.PrimaryBodyScaleFactor)
    primaryScale = opts.PrimaryBodyScaleFactor;
else
    primaryScale = NaN;
end
if isfinite(opts.SecondaryBodyScaleFactor)
    secondaryScale = opts.SecondaryBodyScaleFactor;
else
    secondaryScale = NaN;
end

key = lower(regexprep(char(systemName),'[^a-zA-Z0-9]',''));
switch key
    case {'sunearth','sunearthsystem'}
        if ~isfinite(primaryScale), primaryScale = 6; end
        if ~isfinite(secondaryScale), secondaryScale = 400; end
    otherwise
        if ~isfinite(primaryScale), primaryScale = 1; end
        if ~isfinite(secondaryScale), secondaryScale = 1; end
end
end


function hide_from_legend(h)
if isempty(h) || ~isgraphics(h), return, end
try
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
catch
    set(h,'HandleVisibility','off');
end
end


function apply_white_plot_theme(fig,axesList,legendList,hSuper)
set(fig,'Color','w');
for ax = axesList
    if ~isgraphics(ax), continue, end
    set(ax,'Color','w','XColor','k','YColor','k','ZColor','k');
    ax.Title.Color = 'k';
    ax.XLabel.Color = 'k';
    ax.YLabel.Color = 'k';
    ax.ZLabel.Color = 'k';
end
for lgd = legendList
    if ~isgraphics(lgd), continue, end
    set(lgd,'TextColor','k','Color','w','EdgeColor',[0.75 0.75 0.75]);
end
if isgraphics(hSuper)
    hSuper.Color = 'k';
end
end


%% =========================================================================
% Celestial-frame transforms
% =========================================================================

function out = transform_ecliptic_vectors(vEcl,frame)
vEcl = normalize_rows_local(double(vEcl));

switch frame
    case "ecliptic"
        out = vEcl;

    case "equatorial"
        out = ecliptic_to_equatorial(vEcl);

    case "galactic"
        vEq = ecliptic_to_equatorial(vEcl);
        R_eq2gal = [ ...
           -0.0548755604162154  -0.8734370902348850  -0.4838350155487132
            0.4941094278755837  -0.4448296299600112   0.7469822444972189
           -0.8676661490190047  -0.1980763734312015   0.4559837761750669];
        out = (R_eq2gal*vEq.').';
end

out = normalize_rows_local(out);
end


function vEq = ecliptic_to_equatorial(vEcl)
epsJ2000 = deg2rad(84381.448/3600);
c = cos(epsJ2000);
s = sin(epsJ2000);
R = [1 0 0; 0 c -s; 0 s c];
vEq = (R*vEcl.').';
end


function A = normalize_rows_local(A)
n = vecnorm(A,2,2);
good = isfinite(n) & n > 0;
A(good,:) = A(good,:) ./ n(good);
A(~good,:) = NaN;
end


%% =========================================================================
% Sphere styling
% =========================================================================

function plot_sphere_reference_grid(ax)
phi = linspace(0,2*pi,361);

% Equator of the selected display frame.
plot3(ax,cos(phi),sin(phi),zeros(size(phi)), ...
    ':','Color',[0.45 0.45 0.45],'LineWidth',0.9, ...
    'HandleVisibility','off');

% Meridians every 45 deg.
lat = linspace(-pi/2,pi/2,181);
for lonDeg = 0:45:315
    lon = deg2rad(lonDeg);
    x = cos(lat)*cos(lon);
    y = cos(lat)*sin(lon);
    z = sin(lat);
    plot3(ax,x,y,z,':','Color',[0.76 0.76 0.76], ...
        'LineWidth',0.55,'HandleVisibility','off');
end

% Latitude circles at +/-30 and +/-60 deg.
for latDeg = [-60 -30 30 60]
    la = deg2rad(latDeg);
    x = cos(la)*cos(phi);
    y = cos(la)*sin(phi);
    z = sin(la)*ones(size(phi));
    plot3(ax,x,y,z,':','Color',[0.76 0.76 0.76], ...
        'LineWidth',0.55,'HandleVisibility','off');
end
end


function label = sphere_axis_label(frame,axisName)
switch frame
    case "ecliptic"
        prefix = '\hat{e}_{';
        switch axisName
            case 'x', suffix = '\lambda,0}';
            case 'y', suffix = '\lambda,90}';
            otherwise, suffix = '\beta,+}';
        end
        label = ['$' prefix suffix '$'];

    case "equatorial"
        switch axisName
            case 'x', label = '$\hat{x}_{\rm eq}$';
            case 'y', label = '$\hat{y}_{\rm eq}$';
            otherwise, label = '$\hat{z}_{\rm eq}$';
        end

    case "galactic"
        switch axisName
            case 'x', label = '$\hat{x}_{\rm gal}$';
            case 'y', label = '$\hat{y}_{\rm gal}$';
            otherwise, label = '$\hat{z}_{\rm gal}$';
        end
end
end


%% =========================================================================
% LaTeX-safe plain text
% =========================================================================

function out = latex_text(in)
out = string(in);
out = replace(out,"_","\_");
out = replace(out,"%","\%");
out = replace(out,"&","\&");
out = replace(out,"#","\#");
end


%% =========================================================================
% Labels
% =========================================================================

function [pName,sName] = body_names_from_overlays(overlays)
pName = "Primary";
sName = "Secondary";
if isfield(overlays(1),'embedding') && isstruct(overlays(1).embedding)
    E = overlays(1).embedding;
    if isfield(E,'primaryName') && strlength(string(E.primaryName)) > 0
        pName = string(E.primaryName);
    end
    if isfield(E,'secondaryName') && strlength(string(E.secondaryName)) > 0
        sName = string(E.secondaryName);
    end
end
end


function name = system_display_name(sys,fallback)
if isstruct(sys) && isfield(sys,'name') && strlength(string(sys.name)) > 0
    name = string(sys.name);
else
    name = string(fallback);
end
end
