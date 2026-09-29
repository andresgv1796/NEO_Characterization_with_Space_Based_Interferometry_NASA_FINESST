function overlays = make_cr3bp_floquet_normal_overlays( ...
    systemName,poFile,parentRows,epochUTC,opts)
%MAKE_CR3BP_FLOQUET_NORMAL_OVERLAYS
%
% Build ideal-CR3BP inertial celestial-sky tracks for the LINEAR transported
% Floquet formation normals of one or more periodic-orbit parents.
%
% The function is deliberately family-agnostic and system-agnostic.  The
% selected periodic-orbit CSV is passed directly to choose_periodic_orbit,
% while systemName is passed directly to cr3bp_system_parameters.
%
% BASIC USE
% ---------
% Earth-Moon, several parents from one halo family:
%
%   epoch = datetime('now','TimeZone','UTC');
%
%   ov = make_cr3bp_floquet_normal_overlays( ...
%       "EarthMoon", ...
%       "L1_northern_halo_orbits_JPL_IC_Earth_Moon.csv", ...
%       [82 125 197], ...
%       epoch, ...
%       CenterMode="vertical");
%
%   plot_JPL_asteroid_sky( ...
%       "geocentric",epoch, ...
%       CoordinateFrame="ecliptic", ...
%       NormalOverlay=ov);
%
% Sun-Earth works identically; only the system and family file change:
%
%   ov = make_cr3bp_floquet_normal_overlays( ...
%       "SunEarth",sunEarthFamilyFile,[i1 i2 i3],epoch);
%
% CENTER MODE
% -----------
% CenterMode may be one scalar string applied to every parent, or a string
% array with one entry per ParentRows element:
%
%   CenterMode="vertical"
%   CenterMode=["vertical" "planar" "vertical"]
%
% IDEAL SKY EMBEDDING
% -------------------
% The CR3BP transformation itself is
%
%   n_I(tau) = C0 * R3(tau) * n_S(tau),
%
% where tau is nondimensional CR3BP time and R3 is the active +z rotation.
%
% EmbeddingMode="ideal_epoch_anchor" (default)
%   C0 is constructed at epochUTC from the instantaneous primary->secondary
%   direction and orbital angular-momentum direction returned by Horizons.
%   After that one-time anchoring, the CR3BP frame rotates ideally at unit
%   nondimensional rate.  This is intentionally NOT an ephemeris propagation.
%
% EmbeddingMode="identity"
%   Treat the model inertial basis as J2000 ecliptic directly: C0 = I.
%
% EmbeddingMode="custom"
%   Supply InitialSynodicToEcliptic, a 3x3 orthonormal matrix whose columns
%   are the CR3BP synodic x,y,z axes at tau=0 expressed in J2000 ecliptic.
%
% For known systems the Horizons pair is inferred automatically:
%   EarthMoon -> Earth (399), Moon (301)
%   SunEarth  -> Sun (10), Earth (399)
%
% For another system, supply PrimaryHorizonsID and SecondaryHorizonsID.
%
% OUTPUT
% ------
% overlays is a struct array, one element per parent.  Fields used directly
% by plot_JPL_asteroid_sky include:
%
%   vectorsEclipticJ2000
%   label
%   color
%   lineWidth
%   lineStyle
%   showAntipode
%   antipodeLineStyle
%   showEndpoints
%   showLabel
%
% The same structures also retain the parent periodic-orbit trajectory,
% Floquet geometry, synodic normal history, timing, and metadata so they can
% be reused by plot_cr3bp_floquet_orbits_and_sky without recomputation.
%
% Required FINESST functions on path:
%   cr3bp_system_parameters
%   choose_periodic_orbit
%   monodromy_CR3BP
%   gmos_center_eigpair_from_eigs
%   transport_floquet_center_geometry_CR3BP
%   CR3BP_mu
%

arguments
    systemName (1,1) string = "EarthMoon"
    poFile (1,1) string = "L1_northern_halo_orbits_JPL_IC_Earth_Moon.csv"
    parentRows (1,:) double {mustBeInteger,mustBePositive} = 116
    epochUTC (1,1) datetime = datetime('now','TimeZone','UTC')

    opts.CenterMode string = "vertical"

    opts.NPrimaryPeriods (1,1) double ...
        {mustBeInteger,mustBePositive} = 1

    opts.SamplesPerPeriod (1,1) double ...
        {mustBeInteger,mustBePositive} = 600

    opts.FloquetPhaseSamples (1,1) double ...
        {mustBeInteger,mustBePositive} = 61

    opts.OrbitSamples (1,1) double ...
        {mustBeInteger,mustBePositive} = 800

    opts.NormalReferenceSynodic (1,3) double = [1 0 0]

    opts.EmbeddingMode (1,1) string ...
        {mustBeMember(opts.EmbeddingMode, ...
        ["ideal_epoch_anchor","identity","custom"])} ...
        = "ideal_epoch_anchor"

    opts.PrimaryHorizonsID (1,1) string = ""
    opts.SecondaryHorizonsID (1,1) string = ""
    opts.InitialSynodicToEcliptic double = []

    % The JPL periodic-orbit initial conditions are tied to the mass ratio
    % stored in their family file. Use that value by default so the
    % catalog orbit remains periodic under propagation.
    opts.MassRatioSource (1,1) string ...
        {mustBeMember(opts.MassRatioSource,["catalog","system","manual"])} ...
        = "catalog"

    opts.MassRatio (1,1) double = NaN

    opts.CheckOrbitClosure (1,1) logical = true
    opts.OrbitClosureTolerance (1,1) double {mustBePositive} = 1e-8

    opts.ShowAntipode (1,1) logical = true
    opts.LineWidth (1,1) double {mustBePositive} = 2.4

    % [] -> automatic colors according to ColorMode. Otherwise nParents x 3 RGB.
    opts.Colors double = []
    opts.ColorMode (1,1) string ...
        {mustBeMember(opts.ColorMode,["lines","gradient"])} = "lines"
    opts.ColorMapName (1,1) string = "turbo"
    opts.ColorParameter (1,1) string ...
        {mustBeMember(opts.ColorParameter,["sequence","parentRow"])} = "sequence"

    % Empty -> automatic labels. Scalar -> same prefix for all. nParents
    % entries -> one label per parent.
    opts.Labels string = strings(0,1)

    opts.Verbose (1,1) logical = true
end

%% Normalize inputs
if isempty(epochUTC.TimeZone)
    epochUTC.TimeZone = 'UTC';
else
    epochUTC.TimeZone = 'UTC';
end

parentRows = parentRows(:).';
nParents = numel(parentRows);
centerModes = expand_string_option(opts.CenterMode,nParents,'CenterMode');
labels = expand_optional_labels(opts.Labels,nParents);
[colors,colorValues,colorLabel] = resolve_colors( ...
    opts.Colors,nParents,parentRows,opts.ColorMode, ...
    opts.ColorMapName,opts.ColorParameter);

for k = 1:nParents
    if ~ismember(centerModes(k),["vertical","planar"])
        error('CenterMode entries must be "vertical" or "planar".');
    end
end

requiredFunctions = [ ...
    "cr3bp_system_parameters"
    "choose_periodic_orbit"
    "monodromy_CR3BP"
    "gmos_center_eigpair_from_eigs"
    "transport_floquet_center_geometry_CR3BP"
    "CR3BP_mu"
    ];

for k = 1:numel(requiredFunctions)
    if exist(requiredFunctions(k),'file') ~= 2
        error(['Required FINESST function "%s" is not on the MATLAB path. ' ...
               'Add the FINESST Code directory before calling this helper.'], ...
               requiredFunctions(k));
    end
end

if ~isfile(poFile)
    error('Periodic-orbit family file not found: %s',poFile);
end

sysNominal = cr3bp_system_parameters(char(systemName));

% IMPORTANT: JPL family initial conditions are generated for the mass ratio
% stored in the catalog itself. A small mismatch in mu is enough to make an
% otherwise periodic orbit appear to drift or fail to close. Synchronize
% the dynamics with the selected family before any monodromy or propagation.
[sys,massRatioMeta] = synchronize_system_mass_ratio( ...
    sysNominal,poFile,parentRows,opts.MassRatioSource,opts.MassRatio,opts.Verbose);

%% One common ideal inertial embedding for every parent in this call
[C0,embeddingMeta] = resolve_initial_embedding( ...
    systemName,epochUTC,opts.EmbeddingMode, ...
    opts.PrimaryHorizonsID,opts.SecondaryHorizonsID, ...
    opts.InitialSynodicToEcliptic,opts.Verbose);

%% Console header
vfprintf_local(opts.Verbose,'\n============================================================\n');
vfprintf_local(opts.Verbose,' GENERIC CR3BP FLOQUET NORMAL SKY TRACKS\n');
vfprintf_local(opts.Verbose,'============================================================\n');
vfprintf_local(opts.Verbose,'System        : %s\n',systemName);
vfprintf_local(opts.Verbose,'Family file   : %s\n',poFile);
vfprintf_local(opts.Verbose,'mu source     : %s\n',massRatioMeta.source);
vfprintf_local(opts.Verbose,'mu used       : %.16e\n',sys.mu);
if isfinite(massRatioMeta.systemMu)
    vfprintf_local(opts.Verbose,'mu(system)    : %.16e\n',massRatioMeta.systemMu);
end
if isfinite(massRatioMeta.catalogMu)
    vfprintf_local(opts.Verbose,'mu(catalog)   : %.16e\n',massRatioMeta.catalogMu);
end
vfprintf_local(opts.Verbose,'Parent rows   : ');
if opts.Verbose
    fprintf('%d ',parentRows);
    fprintf('\n');
end
vfprintf_local(opts.Verbose,'Periods       : %d\n',opts.NPrimaryPeriods);
vfprintf_local(opts.Verbose,'Embedding     : %s\n',opts.EmbeddingMode);
vfprintf_local(opts.Verbose,'Epoch UTC     : %s\n', ...
    char(string(epochUTC,'yyyy-MM-dd HH:mm:ss')));

%% Preallocate output struct
emptyOverlay = struct( ...
    'vectorsEclipticJ2000',[], ...
    'tau',[], ...
    'elapsedDays',[], ...
    'epochUTC',epochUTC, ...
    'systemName',systemName, ...
    'system',sys, ...
    'massRatio',sys.mu, ...
    'massRatioMeta',massRatioMeta, ...
    'familyFile',poFile, ...
    'familyName',string(fileparts_name(poFile)), ...
    'parentRow',NaN, ...
    'parentOrbit',struct(), ...
    'parentPeriod_TU',NaN, ...
    'parentPeriod_days',NaN, ...
    'centerMode',"", ...
    'centerPairIndex',NaN, ...
    'centerPlanarFraction',NaN, ...
    'centerVerticalFraction',NaN, ...
    'centerCharacter',"", ...
    'centerRequestedModeDominant',false, ...
    'floquetEigenvalue',NaN, ...
    'synodicNormals',[], ...
    'normalFieldPath',"", ...
    'rawNormalSamples',NaN, ...
    'resampledNormalSamples',NaN, ...
    'normalTimePath',"", ...
    'floquetGeom',struct(), ...
    'initialSynodicToEcliptic',C0, ...
    'embedding',embeddingMeta, ...
    'orbitTau',[], ...
    'orbitSynodic',[], ...
    'orbitClosureError',NaN, ...
    'label',"", ...
    'color',[0 0 0], ...
    'colorMode',opts.ColorMode, ...
    'colorMapName',opts.ColorMapName, ...
    'colorParameter',opts.ColorParameter, ...
    'colorValue',NaN, ...
    'colorLabel',"", ...
    'lineWidth',opts.LineWidth, ...
    'lineStyle','-', ...
    'showAntipode',opts.ShowAntipode, ...
    'antipodeLineStyle','--', ...
    'showEndpoints',true, ...
    'showLabel',true);

overlays = repmat(emptyOverlay,1,nParents);
nValid = 0;
nSkipped = 0;

odeOpts = odeset('RelTol',1e-12,'AbsTol',1e-12);

%% Parent loop
for iParent = 1:nParents

    parentRow = parentRows(iParent);
    centerMode = centerModes(iParent);

    vfprintf_local(opts.Verbose,'\n------------------------------------------------------------\n');
    vfprintf_local(opts.Verbose,'Parent %d / %d : row %d\n', ...
        iParent,nParents,parentRow);
    vfprintf_local(opts.Verbose,'Center mode          : %s\n',centerMode);

    po = choose_periodic_orbit(char(poFile),parentRow,sys);

    %% Monodromy + selected center pair
    M_eigs = monodromy_CR3BP([0,po.T],po.x0(:),sys.mu);
    V = M_eigs(:,1:6);
    D = M_eigs(:,7:12);

    centerOpts = struct();
    centerOpts.tolUnit = 1e-6;
    centerOpts.tolImag = 1e-10;
    centerOpts.warnMismatch = false;

    try
        cen = gmos_center_eigpair_from_eigs( ...
            V,D,char(centerMode),centerOpts);
    catch ME
        warning(['Skipping parent row %d: could not extract the requested ' ...
                 'center subspace (%s). Reason: %s'], ...
                 parentRow,char(centerMode),ME.message);
        nSkipped = nSkipped + 1;
        continue
    end

    if ~isstruct(cen) || ~isfield(cen,'w') || isempty(cen.w) || ...
            any(~isfinite(cen.w(:)))
        warning(['Skipping parent row %d: no usable %s center subspace ' ...
                 'was returned by gmos_center_eigpair_from_eigs.'], ...
                 parentRow,char(centerMode));
        nSkipped = nSkipped + 1;
        continue
    end

    %% Linear Floquet geometry over one parent period
    nFloquetSamples = opts.SamplesPerPeriod + 1;

    floquetGeom = transport_floquet_center_geometry_CR3BP( ...
        [0,po.T],po.x0(:),sys.mu,cen.w, ...
        nFloquetSamples,odeOpts,opts.FloquetPhaseSamples);

    [nRaw,pathUsed,tRaw,timePath] = extract_synodic_normal_track( ...
        floquetGeom,nFloquetSamples,po.T);

    nRaw = normalize_rows(nRaw);
    nRaw = enforce_axial_sign_continuity( ...
        nRaw,opts.NormalReferenceSynodic);

    rawNormalSamples = size(nRaw,1);

    % transport_floquet_center_geometry_CR3BP may return its normal track on
    % the adaptive ODE mesh rather than exactly on SamplesPerPeriod+1 points.
    % Resample the ordered, sign-continuous unit-normal track onto the
    % requested uniform CR3BP time grid so every parent has the same number
    % of samples and the ideal R3(tau) embedding uses the correct time tags.
    [tTemplate,nTemplate,resampleTimeSource] = resample_normal_track( ...
        nRaw,tRaw,po.T,nFloquetSamples,timePath);

    %% Repeat periodic axial template
    [tau,nSyn] = repeat_axial_template( ...
        tTemplate,nTemplate,po.T,opts.NPrimaryPeriods);

    %% Ideal synodic -> inertial J2000-ecliptic embedding
    nEcl = ideal_synodic_vectors_to_ecliptic(nSyn,tau,C0);

    %% Parent periodic orbit for the companion 3-D visualization
    tOrbit = linspace(0,po.T,opts.OrbitSamples).';
    [~,XOrbit] = ode89_or_fallback( ...
        @(t,x) CR3BP_mu(t,x,sys.mu), ...
        tOrbit,po.x0(:),odeOpts);

    orbitClosureError = norm(XOrbit(end,1:6).' - po.x0(:));

    if opts.CheckOrbitClosure && orbitClosureError > opts.OrbitClosureTolerance
        warning( ...
            ['Parent row %d closes with ||x(T)-x(0)||_2 = %.3e, larger ' ...
             'than the requested tolerance %.3e. Check that the family file ' ...
             'and mass ratio correspond to the same CR3BP system.'], ...
            parentRow,orbitClosureError,opts.OrbitClosureTolerance);
    end

    %% Package
    ov = emptyOverlay;
    ov.vectorsEclipticJ2000 = nEcl;
    ov.tau = tau;

    % Prefer the dimensional period supplied with the selected JPL orbit.
    % This makes the displayed time scale consistent with the same catalog
    % that supplied x0, T, and mu. Fall back to the system scale otherwise.
    timeUnitDays = NaN;
    if isfield(po,'T_days') && isfinite(po.T_days) && po.T > 0
        timeUnitDays = po.T_days/po.T;
    elseif isfield(sys,'Tstar_days') && isfinite(sys.Tstar_days)
        timeUnitDays = sys.Tstar_days;
    end

    if isfinite(timeUnitDays)
        ov.elapsedDays = tau*timeUnitDays;
    else
        ov.elapsedDays = NaN(size(tau));
    end

    ov.parentRow = parentRow;
    ov.parentOrbit = po;
    ov.parentPeriod_TU = po.T;
    if isfield(po,'T_days') && isfinite(po.T_days)
        ov.parentPeriod_days = po.T_days;
    elseif isfinite(timeUnitDays)
        ov.parentPeriod_days = po.T*timeUnitDays;
    end
    ov.centerMode = centerMode;
    if isfield(cen,'idx')
        ov.centerPairIndex = cen.idx;
    end
    if isfield(cen,'proj') && isstruct(cen.proj)
        if isfield(cen.proj,'planarFraction')
            ov.centerPlanarFraction = cen.proj.planarFraction;
        end
        if isfield(cen.proj,'verticalFraction')
            ov.centerVerticalFraction = cen.proj.verticalFraction;
        end
        if isfield(cen.proj,'character')
            ov.centerCharacter = string(cen.proj.character);
        end
        if isfield(cen.proj,'requestedModeDominant')
            ov.centerRequestedModeDominant = logical(cen.proj.requestedModeDominant);
        end
    end
    if isfield(cen,'lambda')
        ov.floquetEigenvalue = cen.lambda;
    end
    ov.synodicNormals = nSyn;
    ov.normalFieldPath = string(pathUsed);
    ov.rawNormalSamples = rawNormalSamples;
    ov.resampledNormalSamples = size(nTemplate,1);
    ov.normalTimePath = string(resampleTimeSource);
    ov.floquetGeom = floquetGeom;
    ov.orbitTau = tOrbit;
    ov.orbitSynodic = XOrbit(:,1:3);
    ov.orbitClosureError = orbitClosureError;
    ov.color = colors(iParent,:);
    ov.colorMode = opts.ColorMode;
    ov.colorMapName = opts.ColorMapName;
    ov.colorParameter = opts.ColorParameter;
    ov.colorValue = colorValues(iParent);
    ov.colorLabel = colorLabel;

    if strlength(labels(iParent)) > 0
        ov.label = labels(iParent);
    else
        ov.label = sprintf('%s row %d (%s)', ...
            short_family_name(poFile),parentRow,centerMode);
    end

    nValid = nValid + 1;
    overlays(nValid) = ov;

    vfprintf_local(opts.Verbose,'TPrimary             : %.12f TU\n',po.T);
    if isfinite(ov.parentPeriod_days)
        vfprintf_local(opts.Verbose,'TPrimary             : %.6f days\n', ...
            ov.parentPeriod_days);
    end
    if isfinite(ov.centerVerticalFraction)
        vfprintf_local(opts.Verbose,'Center character     : %s (f_v = %.3f)\n', ...
            char(ov.centerCharacter),ov.centerVerticalFraction);
    end
    vfprintf_local(opts.Verbose,'Normal field         : %s\n',pathUsed);
    vfprintf_local(opts.Verbose,'Raw normal samples   : %d\n',rawNormalSamples);
    vfprintf_local(opts.Verbose,'Normal time source   : %s\n',resampleTimeSource);
    vfprintf_local(opts.Verbose,'Resampled / period   : %d\n',size(nTemplate,1));
    vfprintf_local(opts.Verbose,'Returned samples     : %d\n',size(nEcl,1));
    vfprintf_local(opts.Verbose,'PO closure error     : %.3e\n',orbitClosureError);

end

overlays = overlays(1:nValid);

if nValid == 0
    error(['None of the requested parent rows produced a usable center ' ...
           'subspace for the selected CenterMode.']);
end

vfprintf_local(opts.Verbose,'\n============================================================\n');
vfprintf_local(opts.Verbose,'Returned %d color-linked overlay(s).\n',nValid);
if nSkipped > 0
    vfprintf_local(opts.Verbose,'Skipped %d parent row(s) with no usable center subspace.\n',nSkipped);
end
vfprintf_local(opts.Verbose,'============================================================\n\n');

end


%% =========================================================================
% Catalog mass-ratio synchronization
% =========================================================================

function [sys,meta] = synchronize_system_mass_ratio( ...
    sysIn,poFile,parentRows,source,manualMu,verbose)

sys = sysIn;

meta = struct();
meta.source = string(source);
meta.systemMu = NaN;
meta.catalogMu = NaN;
meta.catalogSpread = NaN;
meta.catalogColumn = "";

if isstruct(sysIn) && isfield(sysIn,'mu') && isfinite(sysIn.mu)
    meta.systemMu = double(sysIn.mu);
end

switch source

    case "system"
        if ~isfinite(meta.systemMu)
            error('MassRatioSource="system" requires sys.mu to be finite.');
        end
        muUse = meta.systemMu;

    case "manual"
        if ~isfinite(manualMu) || manualMu <= 0 || manualMu >= 0.5
            error(['MassRatioSource="manual" requires MassRatio to be a ' ...
                   'finite scalar in (0,0.5).']);
        end
        muUse = double(manualMu);

    case "catalog"
        [muCatalog,spread,columnName] = ...
            read_catalog_mass_ratio(poFile,parentRows);

        meta.catalogMu = muCatalog;
        meta.catalogSpread = spread;
        meta.catalogColumn = columnName;
        muUse = muCatalog;

    otherwise
        error('Unknown MassRatioSource "%s".',source);
end

if ~isfinite(muUse) || muUse <= 0 || muUse >= 0.5
    error('Resolved CR3BP mass ratio is invalid: %.16e.',muUse);
end

sys.mu = muUse;

if isfinite(meta.systemMu)
    delta = abs(muUse-meta.systemMu);
    rel = delta/max(abs(muUse),eps);

    if rel > 1e-12
        vfprintf_local(verbose, ...
            ['Mass-ratio synchronization: overriding sys.mu %.16e with ' ...
             '%.16e from %s (relative difference %.3e).\n'], ...
            meta.systemMu,muUse,source,rel);
    end
end

end


function [muCatalog,spread,columnName] = ...
    read_catalog_mass_ratio(poFile,parentRows)

try
    T = readtable(char(poFile),'VariableNamingRule','preserve');
catch ME
    error('Could not read periodic-orbit catalog "%s": %s',poFile,ME.message);
end

if any(parentRows < 1) || any(parentRows > height(T))
    error(['At least one requested parent row lies outside catalog "%s" ' ...
           '(height = %d).'],poFile,height(T));
end

names = string(T.Properties.VariableNames);
keys = lower(regexprep(names,'[^a-zA-Z0-9]',''));

idx = find( ...
    contains(keys,'massratio') | ...
    keys == "mu" | ...
    contains(keys,'massparameter'), ...
    1,'first');

% JPL periodic-orbit CSV files used in this project have the mass ratio as
% column 12. Keep this only as a conservative fallback if the header is not
% descriptive.
if isempty(idx) && width(T) >= 12
    candidate = table_column_to_double(T,parentRows,12);
    if all(isfinite(candidate)) && all(candidate > 0) && all(candidate < 0.5)
        idx = 12;
    end
end

if isempty(idx)
    error([ ...
        'Could not identify the JPL catalog mass-ratio column in "%s". ' ...
        'Use MassRatioSource="manual", MassRatio=<mu> if this catalog uses ' ...
        'a nonstandard format.'],poFile);
end

vals = table_column_to_double(T,parentRows,idx);

if any(~isfinite(vals)) || any(vals <= 0) || any(vals >= 0.5)
    error('Catalog mass-ratio values are missing or invalid in column "%s".', ...
        names(idx));
end

muCatalog = median(vals);
spread = max(vals)-min(vals);
columnName = names(idx);

tol = max(1e-14,1e-10*abs(muCatalog));

if spread > tol
    error([ ...
        'Requested parent rows do not share one CR3BP mass ratio. ' ...
        'Catalog spread = %.3e around mu = %.16e.'], ...
        spread,muCatalog);
end

end


function vals = table_column_to_double(T,rows,idx)

raw = T{rows,idx};

if isnumeric(raw) || islogical(raw)
    vals = double(raw(:));
else
    vals = str2double(string(raw(:)));
end

end


%% =========================================================================
% Input expansion / colors / labels
% =========================================================================

function out = expand_string_option(in,n,name)
in = string(in(:));
if isscalar(in)
    out = repmat(in,n,1);
elseif numel(in) == n
    out = in;
else
    error('%s must be scalar or have one entry per parent row.',name);
end
end


function labels = expand_optional_labels(in,n)
in = string(in(:));
if isempty(in)
    labels = repmat("",n,1);
elseif isscalar(in)
    labels = repmat(in,n,1);
elseif numel(in) == n
    labels = in;
else
    error('Labels must be empty, scalar, or have one entry per parent row.');
end
end


function [colors,values,label] = resolve_colors( ...
    colorsIn,n,parentRows,colorMode,colorMapName,colorParameter)

if ~isempty(colorsIn)
    if ~isnumeric(colorsIn) || size(colorsIn,2) ~= 3 || size(colorsIn,1) ~= n
        error('Colors must be empty or an nParents x 3 numeric RGB array.');
    end
    colors = max(0,min(1,double(colorsIn)));
    values = (1:n).';
    label = "custom color index";
    return
end

switch colorParameter
    case "sequence"
        values = (1:n).';
        label = "family sequence";
    case "parentRow"
        values = parentRows(:);
        label = "parent row";
end

switch colorMode
    case "lines"
        colors = lines(n);

    case "gradient"
        % Color continuously with the selected scalar parameter, rather than
        % simply taking equally spaced colors by array index. This keeps the
        % gradient meaningful when parentRows are irregularly spaced.
        cmap = local_named_colormap(colorMapName,256);
        vmin = min(values);
        vmax = max(values);
        if n == 1 || vmax <= vmin
            q = 0.5*ones(n,1);
        else
            q = (values-vmin)/(vmax-vmin);
        end
        idx = 1 + round(q*(size(cmap,1)-1));
        idx = max(1,min(size(cmap,1),idx));
        colors = cmap(idx,:);
end
end


function cmap = local_named_colormap(name,m)
name = lower(char(name));
switch name
    case 'parula'
        cmap = parula(m);
    case 'turbo'
        cmap = turbo(m);
    case 'hsv'
        cmap = hsv(m);
    case 'jet'
        cmap = jet(m);
    case 'hot'
        cmap = hot(m);
    case 'cool'
        cmap = cool(m);
    case 'spring'
        cmap = spring(m);
    case 'summer'
        cmap = summer(m);
    case 'autumn'
        cmap = autumn(m);
    case 'winter'
        cmap = winter(m);
    case 'gray'
        cmap = gray(m);
    case 'bone'
        cmap = bone(m);
    case 'copper'
        cmap = copper(m);
    case 'pink'
        cmap = pink(m);
    case 'lines'
        cmap = lines(m);
    otherwise
        warning('Unknown ColorMapName "%s". Falling back to turbo.',string(name));
        cmap = turbo(m);
end
end


function [C0,meta] = resolve_initial_embedding( ...
    systemName,epochUTC,embeddingMode,primaryID,secondaryID,Ccustom,verbose)

meta = struct();
meta.mode = embeddingMode;
meta.primaryHorizonsID = "";
meta.secondaryHorizonsID = "";
meta.primaryName = "Primary";
meta.secondaryName = "Secondary";
meta.description = "";

switch embeddingMode

    case "identity"
        C0 = eye(3);
        meta.description = ...
            "Model inertial axes identified directly with J2000 ecliptic axes.";

    case "custom"
        if ~isequal(size(Ccustom),[3 3])
            error(['EmbeddingMode="custom" requires ' ...
                   'InitialSynodicToEcliptic as a 3x3 matrix.']);
        end
        C0 = validate_rotation_matrix(Ccustom,'InitialSynodicToEcliptic');
        meta.description = ...
            "User-supplied tau=0 synodic-to-J2000-ecliptic rotation.";

    case "ideal_epoch_anchor"
        [primaryID,secondaryID,pName,sName] = resolve_horizons_pair( ...
            systemName,primaryID,secondaryID);

        jdUTC = posixtime(epochUTC)/86400 + 2440587.5;

        vfprintf_local(verbose, ...
            'Querying Horizons for %s -> %s epoch anchor...\n',pName,sName);

        [rPS,vPS] = get_relative_state_Horizons( ...
            char(secondaryID),char(primaryID),jdUTC);

        xHat = rPS/norm(rPS);
        zHat = cross(rPS,vPS);

        if norm(zHat) < 1e-14
            error('Primary-secondary state produced a degenerate orbital normal.');
        end

        zHat = zHat/norm(zHat);
        yHat = cross(zHat,xHat);
        yHat = yHat/norm(yHat);
        xHat = cross(yHat,zHat);
        xHat = xHat/norm(xHat);

        C0 = [xHat(:),yHat(:),zHat(:)];
        C0 = validate_rotation_matrix(C0,'Horizons epoch anchor');

        meta.primaryHorizonsID = secondary_or_empty(primaryID);
        meta.secondaryHorizonsID = secondary_or_empty(secondaryID);
        meta.primaryName = pName;
        meta.secondaryName = sName;
        meta.description = sprintf( ...
            ['Ideal CR3BP anchored at epoch to the instantaneous %s->%s ' ...
             'direction and angular momentum; thereafter C=C0*R3(tau).'], ...
             pName,sName);
end

meta.C0 = C0;
end


function s = secondary_or_empty(x)
s = string(x);
end


function [primaryID,secondaryID,pName,sName] = resolve_horizons_pair( ...
    systemName,primaryID,secondaryID)

if strlength(primaryID) > 0 || strlength(secondaryID) > 0
    if strlength(primaryID) == 0 || strlength(secondaryID) == 0
        error(['PrimaryHorizonsID and SecondaryHorizonsID must either both ' ...
               'be supplied or both be left empty.']);
    end
    pName = "Primary";
    sName = "Secondary";
    return
end

key = lower(regexprep(char(systemName),'[^a-zA-Z0-9]',''));

switch key
    case {'earthmoon','earthmoonsystem'}
        primaryID = "399";
        secondaryID = "301";
        pName = "Earth";
        sName = "Moon";

    case {'sunearth','sunearthsystem'}
        primaryID = "10";
        secondaryID = "399";
        pName = "Sun";
        sName = "Earth";

    otherwise
        error([ ...
            'No automatic Horizons body-pair mapping is defined for system "%s". ' ...
            'Supply PrimaryHorizonsID and SecondaryHorizonsID, or use ' ...
            'EmbeddingMode="identity"/"custom".'],systemName);
end
end


function C = validate_rotation_matrix(C,name)
C = double(C);
errOrth = norm(C.'*C-eye(3),'fro');
detC = det(C);
if errOrth > 1e-8 || abs(detC-1) > 1e-8
    error('%s must be a proper orthonormal rotation matrix.',name);
end
end


%% =========================================================================
% Floquet normal extraction
% =========================================================================

function [n,pathUsed,tRaw,timePath] = extract_synodic_normal_track(S,nExpected,T)

tRaw = [];
timePath = "";

% Preferred explicit field used by the current FINESST helper.
if isstruct(S) && isfield(S,'nSynodic') && isnumeric(S.nSynodic)
    A = orient_Nx3(S.nSynodic);
    if size(A,1) >= 2
        n = A;
        pathUsed = 'floquetGeom.nSynodic';
        [tRaw,timePath] = find_matching_time_vector(S,size(A,1),T);
        return
    end
end

candidates = collect_numeric_vector_candidates(S,"floquetGeom");

if isempty(candidates)
    error(['Could not find an N x 3 numeric vector track inside the output ' ...
           'of transport_floquet_center_geometry_CR3BP.']);
end

bestScore = -Inf;
bestIdx = NaN;

for k = 1:numel(candidates)
    A = candidates(k).data;
    path = lower(candidates(k).path);

    if size(A,2) ~= 3
        continue
    end

    norms = vecnorm(A,2,2);
    finiteNorms = norms(isfinite(norms));
    if isempty(finiteNorms)
        continue
    end

    score = 0;
    if contains(path,'synod'), score = score + 40; end
    if contains(path,'normal'), score = score + 35; end
    if contains(path,'nhat') || endsWith(path,'.n'), score = score + 20; end
    if contains(path,'inert'), score = score - 40; end
    if size(A,1) == nExpected, score = score + 25; end

    medNorm = median(finiteNorms);
    if abs(medNorm-1) < 0.05
        score = score + 15;
    elseif medNorm > 0
        score = score + 2;
    end

    if score > bestScore
        bestScore = score;
        bestIdx = k;
    end
end

if ~isfinite(bestIdx) || bestScore < 20
    msg = "Available N x 3 candidates were:";
    for k = 1:numel(candidates)
        if size(candidates(k).data,2) == 3
            msg = msg + newline + "  " + candidates(k).path + ...
                sprintf(' [%d x 3]',size(candidates(k).data,1));
        end
    end
    error(['Could not confidently identify the synodic Floquet-normal ' ...
           'track.' newline '%s'],msg);
end

n = candidates(bestIdx).data;
pathUsed = candidates(bestIdx).path;
[tRaw,timePath] = find_matching_time_vector(S,size(n,1),T);
end


function [tUniform,nUniform,timeSource] = resample_normal_track( ...
    nRaw,tRaw,T,nExpected,timePath)

nRaw = normalize_rows(nRaw);

if isempty(tRaw)
    tWork = linspace(0,T,size(nRaw,1)).';
    timeSource = "uniform-index fallback";
else
    tWork = double(tRaw(:));
    good = isfinite(tWork) & all(isfinite(nRaw),2);
    tWork = tWork(good);
    nRaw = nRaw(good,:);

    [tWork,ia] = unique(tWork,'stable');
    nRaw = nRaw(ia,:);

    if numel(tWork) < 2 || any(diff(tWork) <= 0)
        tWork = linspace(0,T,size(nRaw,1)).';
        timeSource = "uniform-index fallback";
    else
        % If the helper returned an adaptive mesh that spans one parent
        % period, preserve its relative timing but force exact [0,T]
        % endpoints so repeated templates meet cleanly.
        span = tWork(end)-tWork(1);
        if span <= 0
            tWork = linspace(0,T,size(nRaw,1)).';
            timeSource = "uniform-index fallback";
        else
            if abs(tWork(1)) <= 1e-8*max(1,T) && ...
               abs(tWork(end)-T) <= 1e-6*max(1,T)
                % already in CR3BP time
            else
                tWork = (tWork-tWork(1))*(T/span);
            end
            if strlength(string(timePath)) > 0
                timeSource = string(timePath);
            else
                timeSource = "detected adaptive time vector";
            end
        end
    end
end

tUniform = linspace(0,T,nExpected).';

if size(nRaw,1) == nExpected && ...
   max(abs(tWork-tUniform)) <= 1e-12*max(1,T)
    nUniform = nRaw;
else
    nUniform = zeros(nExpected,3);
    for j = 1:3
        nUniform(:,j) = interp1(tWork,nRaw(:,j),tUniform,'pchip');
    end
    nUniform = normalize_rows(nUniform);
    nUniform = enforce_axial_sign_continuity(nUniform,[]);
end
end


function [tBest,pathBest] = find_matching_time_vector(S,nRows,T)

tBest = [];
pathBest = "";
candidates = collect_numeric_scalar_vector_candidates(S,"floquetGeom");

bestScore = -Inf;
for k = 1:numel(candidates)
    v = double(candidates(k).data(:));
    if numel(v) ~= nRows || any(~isfinite(v)) || nRows < 2
        continue
    end
    dv = diff(v);
    if any(dv < 0) || all(dv == 0)
        continue
    end

    path = lower(string(candidates(k).path));
    score = 0;
    if contains(path,"time"), score = score + 50; end
    if contains(path,"tau"), score = score + 45; end
    if endsWith(path,".t") || contains(path,".t_") || contains(path,".tgrid")
        score = score + 40;
    end
    if contains(path,"state") || contains(path,"normal") || contains(path,"sigma")
        score = score - 15;
    end

    span = v(end)-v(1);
    if span > 0
        score = score + 10;
        if abs(v(1)) <= 1e-6*max(1,T), score = score + 5; end
        if abs(v(end)-T) <= 1e-4*max(1,T), score = score + 10; end
    end

    if score > bestScore
        bestScore = score;
        tBest = v;
        pathBest = string(candidates(k).path);
    end
end

if bestScore < 20
    tBest = [];
    pathBest = "";
end
end


function out = collect_numeric_scalar_vector_candidates(S,path)
out = struct('path',{},'data',{});

if isnumeric(S) && isvector(S) && numel(S) >= 2
    out(1).path = char(path);
    out(1).data = S(:);
    return
end

if isstruct(S)
    if numel(S) ~= 1, return, end
    f = fieldnames(S);
    for i = 1:numel(f)
        childPath = string(path) + "." + string(f{i});
        c = collect_numeric_scalar_vector_candidates(S.(f{i}),childPath);
        out = [out,c]; %#ok<AGROW>
    end
elseif iscell(S) && numel(S) <= 10
    for i = 1:numel(S)
        childPath = string(path) + sprintf('{%d}',i);
        c = collect_numeric_scalar_vector_candidates(S{i},childPath);
        out = [out,c]; %#ok<AGROW>
    end
end
end


function A = orient_Nx3(A)
A = double(A);
if size(A,2) == 3
    return
elseif size(A,1) == 3
    A = A.';
else
    A = [];
end
end


function out = collect_numeric_vector_candidates(S,path)
out = struct('path',{},'data',{});

if isnumeric(S)
    A = orient_Nx3(S);
    if ~isempty(A) && size(A,1) >= 2
        out(1).path = char(path);
        out(1).data = A;
    end
    return
end

if isstruct(S)
    if numel(S) ~= 1, return, end
    f = fieldnames(S);
    for i = 1:numel(f)
        childPath = string(path) + "." + string(f{i});
        c = collect_numeric_vector_candidates(S.(f{i}),childPath);
        out = [out,c]; %#ok<AGROW>
    end
elseif iscell(S) && numel(S) <= 10
    for i = 1:numel(S)
        childPath = string(path) + sprintf('{%d}',i);
        c = collect_numeric_vector_candidates(S{i},childPath);
        out = [out,c]; %#ok<AGROW>
    end
end
end


%% =========================================================================
% Axial continuity and repetition
% =========================================================================

function n = enforce_axial_sign_continuity(n,reference)
n = normalize_rows(n);
for k = 2:size(n,1)
    if dot(n(k,:),n(k-1,:)) < 0
        n(k,:) = -n(k,:);
    end
end

if ~isempty(reference) && all(isfinite(reference)) && norm(reference) > 0
    reference = reference/norm(reference);
    if dot(n(1,:),reference) < 0
        n = -n;
    end
end
end


function [tau,nSyn] = repeat_axial_template(tTemplate,nTemplate,T,nPeriods)
seamParity = sign(dot(nTemplate(end,:),nTemplate(1,:)));
if seamParity == 0, seamParity = 1; end

nBlocks = cell(nPeriods,1);
tBlocks = cell(nPeriods,1);

for iPeriod = 0:(nPeriods-1)
    if iPeriod < nPeriods-1
        idx = 1:(size(nTemplate,1)-1);
    else
        idx = 1:size(nTemplate,1);
    end

    signFactor = seamParity^iPeriod;
    nBlocks{iPeriod+1} = signFactor*nTemplate(idx,:);
    tBlocks{iPeriod+1} = tTemplate(idx) + iPeriod*T;
end

nSyn = vertcat(nBlocks{:});
tau = vertcat(tBlocks{:});
end


function nEcl = ideal_synodic_vectors_to_ecliptic(nSyn,tau,C0)
nEcl = NaN(size(nSyn));
for k = 1:numel(tau)
    c = cos(tau(k));
    s = sin(tau(k));
    R3 = [c -s 0; s c 0; 0 0 1];
    nEcl(k,:) = (C0*R3*nSyn(k,:).').';
end
nEcl = normalize_rows(nEcl);
end


function A = normalize_rows(A)
n = vecnorm(A,2,2);
valid = isfinite(n) & n > 0;
A(valid,:) = A(valid,:) ./ n(valid);
A(~valid,:) = NaN;
end


%% =========================================================================
% Periodic-orbit propagation
% =========================================================================

function [t,X] = ode89_or_fallback(rhs,tspan,x0,odeOpts)
if exist('ode89','file') == 2
    [t,X] = ode89(rhs,tspan,x0,odeOpts);
else
    [t,X] = ode113(rhs,tspan,x0,odeOpts);
end
end


%% =========================================================================
% Horizons relative state at one epoch
% =========================================================================

function [rRel,vRel] = get_relative_state_Horizons(targetID,centerID,jdUTC)
baseURL = 'https://ssd.jpl.nasa.gov/api/horizons.api';
webOpts = weboptions('Timeout',120,'ContentType','json');

args = { ...
    'format','json', ...
    'COMMAND',sprintf("'%s'",targetID), ...
    'OBJ_DATA','''NO''', ...
    'MAKE_EPHEM','''YES''', ...
    'EPHEM_TYPE','''VECTORS''', ...
    'CENTER',sprintf("'500@%s'",centerID), ...
    'TLIST',sprintf('''%.9f''',jdUTC), ...
    'TLIST_TYPE','''JD''', ...
    'TIME_TYPE','''UT''', ...
    'REF_PLANE','''ECLIPTIC''', ...
    'REF_SYSTEM','''ICRF''', ...
    'OUT_UNITS','''AU-D''', ...
    'VEC_TABLE','''2''', ...
    'VEC_CORR','''NONE''', ...
    'VEC_LABELS','''YES''' ...
    };

payload = webread_retry_local(baseURL,args,webOpts);
txt = char(payload.result);

iSOE = strfind(txt,'$$SOE');
iEOE = strfind(txt,'$$EOE');
if isempty(iSOE) || isempty(iEOE)
    error('Could not locate Horizons $$SOE/$$EOE block.');
end
block = txt(iSOE(1):iEOE(1));

rRel = [ ...
    horizons_number_local(block,'X'), ...
    horizons_number_local(block,'Y'), ...
    horizons_number_local(block,'Z')];

vRel = [ ...
    horizons_number_local(block,'VX'), ...
    horizons_number_local(block,'VY'), ...
    horizons_number_local(block,'VZ')];
end


function payload = webread_retry_local(url,args,webOpts)
maxAttempts = 5;
for attempt = 1:maxAttempts
    try
        payload = webread(url,args{:},webOpts);
        return
    catch ME
        if attempt == maxAttempts, rethrow(ME), end
        pause(2^(attempt-1));
    end
end
end


function value = horizons_number_local(block,label)
numberPattern = ...
    '([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+\-]?\d+)?)';
pattern = sprintf( ...
    '(?:^|[\\s,])%s\\s*=\\s*%s', ...
    regexptranslate('escape',label),numberPattern);

token = regexp(block,pattern,'tokens','once');
if isempty(token)
    error('Could not parse Horizons quantity "%s".',label);
end
value = str2double(regexprep(token{1},'[Dd]','E'));
end


%% =========================================================================
% Small utilities
% =========================================================================

function s = short_family_name(poFile)
[~,name,~] = fileparts(char(poFile));
s = string(name);
s = regexprep(s,'_orbits_JPL_IC_.*$','');
s = replace(s,'_',' ');
end


function name = fileparts_name(file)
[~,name,~] = fileparts(char(file));
end


function vfprintf_local(tf,varargin)
if tf, fprintf(varargin{:}); end
end
