function [catalog, snapshot, fig] = plot_JPL_asteroid_sky(observerMode, epochUTC, opts)
%PLOT_JPL_ASTEROID_SKY
%
% Instantaneous all-sky map of NASA/JPL asteroids.
%
% Each asteroid is plotted at its instantaneous direction at epochUTC:
%
%   observerMode = "geocentric"
%       Direction as seen from the center of Earth.
%
%   observerMode = "heliocentric"
%       Direction as seen from the center of the Sun.
%
% Visual encoding:
%
%   POSITION     = instantaneous sky longitude/latitude in the selected
%                  equatorial or ecliptic frame
%   MARKER AREA  ∝ asteroid physical projected area = pi D^2 / 4
%   COLOR        = spin rate = 24 / P_rot [rev/day]
%
% The asteroid catalog is obtained from the NASA/JPL Small-Body Database
% Query API. Orbital elements are propagated locally using vectorized
% two-body Kepler propagation.
%
% For geocentric viewing, the current Earth heliocentric state is obtained
% from JPL Horizons and asteroid light time is iterated.  In heliocentric
% mode, the eight major planets can also be plotted from Horizons states.
%
% -------------------------------------------------------------------------
% EXAMPLES
% -------------------------------------------------------------------------
%
% Current geocentric sky:
%
%   plot_JPL_asteroid_sky
%
% Current heliocentric sky:
%
%   plot_JPL_asteroid_sky("heliocentric")
%
% Equatorial map with RA labelled in degrees:
%
%   plot_JPL_asteroid_sky( ...
%       "geocentric", ...
%       datetime('now','TimeZone','UTC'), ...
%       CoordinateFrame="equatorial", ...
%       RAUnits="degrees")
%
% Current heliocentric sky in J2000 ecliptic coordinates, including planets:
%
%   plot_JPL_asteroid_sky( ...
%       "heliocentric", ...
%       datetime('now','TimeZone','UTC'), ...
%       CoordinateFrame="ecliptic", ...
%       ShowPlanets=true)
%
% Specific epoch:
%
%   t = datetime(2027,1,1,0,0,0,'TimeZone','UTC');
%   plot_JPL_asteroid_sky("geocentric",t)
%
% Download the entire selected asteroid population rather than only
% objects with both diameter and rotation period:
%
%   plot_JPL_asteroid_sky( ...
%       "geocentric", ...
%       datetime('now','TimeZone','UTC'), ...
%       DownloadMode="all")
%
% Plot all JPL asteroids (all orbital classes):
%
%   plot_JPL_asteroid_sky( ...
%       "geocentric", ...
%       datetime('now','TimeZone','UTC'), ...
%       Population="all_asteroids")
%
% Quiet mode:
%
%   plot_JPL_asteroid_sky( ...
%       "geocentric", ...
%       datetime('now','TimeZone','UTC'), ...
%       Verbose=false)
%
% Force a new SBDB download:
%
%   plot_JPL_asteroid_sky( ...
%       "geocentric", ...
%       datetime('now','TimeZone','UTC'), ...
%       ForceCatalogRefresh=true)
%
% OUTPUTS
%
%   catalog
%       Cached/downloaded SBDB catalog.
%
%   snapshot
%       Objects actually plotted, including:
%
%           RA_deg
%           Dec_deg
%           eclipticLongitude_deg
%           eclipticLatitude_deg
%           range_AU
%           diameter_km
%           projectedArea_km2
%           rot_per_h
%           spinRate_revDay
%
%   fig
%       Figure handle.
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
%
% 1. SBDB osculating elements are heliocentric J2000-ecliptic elements
%    referenced to their individual epoch of osculation.
%
% 2. Local propagation is two-body Keplerian. It is appropriate for this
%    population-scale visualization, but it should NOT replace Horizons
%    for precision ephemerides of an individual science target.
%
% 3. "Geocentric" coordinates here are astrometric-like directions:
%       - Earth state from Horizons
%       - asteroid light-time iteration
%       - no annual stellar-aberration correction
%
% 4. A small graphical marker-size floor is applied so that small
%    asteroids remain visible. Above that floor, graphical marker area is
%    exactly proportional to physical asteroid projected area.
%
% -------------------------------------------------------------------------
%
% Population options:
%
%   Population="MBA"            main-belt asteroids only
%   Population="all_asteroids"  all JPL asteroid orbital classes
%
% Display filters (applied locally after loading the reusable catalog and
% before propagation):
%
%   DiameterRange_km=[0 Inf]          no diameter filtering (default)
%   SpinRateRange_revDay=[0 Inf]      no spin-rate filtering (default)
%
% Examples:
%
%   DiameterRange_km=[50 Inf]         display D >= 50 km
%   DiameterRange_km=[10 100]         display 10 <= D <= 100 km
%   SpinRateRange_revDay=[0 3]        display spin <= 3 rev/day
%   SpinRateRange_revDay=[2 10]       display 2--10 rev/day
%
% The two filters are combined with AND when both are specified.
%
% Console modes:
%
%   Verbose=true   detailed query, cache, and propagation diagnostics
%   Verbose=false  quiet mode; errors and warnings are still shown
%
% Coordinate frames:
%
%   CoordinateFrame="equatorial"  J2000 Right Ascension / Declination
%   CoordinateFrame="ecliptic"    J2000 ecliptic longitude / latitude
%   CoordinateFrame="galactic"    IAU Galactic longitude / latitude
%
% Equatorial longitude labels:
%
%   RAUnits="hours"    0--24 h right ascension labels
%   RAUnits="degrees"  0--360 deg right ascension labels
%
% Longitude display convention:
%
%   LongitudeDirection="increasing_right"
%       Longitude increases left-to-right. This matches the original
%       milkyway.m RA/Dec map and the reference plots used in this project.
%
%   LongitudeDirection="increasing_left"
%       Conventional "viewed from inside the celestial sphere" sky-map
%       orientation.
%
% Planet overlay:
%
%   ShowPlanets=true       plots Mercury through Neptune in heliocentric mode
%   ShowPlanetLabels=true  labels the plotted planets
%
% Milky Way overlay:
%
%   ShowMilkyWay=true with the default MilkyWaySource="legacy_footprint"
%   requires milkyway.m on the MATLAB path.  The historical irregular
%   RA/Dec footprint is transformed on the celestial sphere into the
%   selected plotting frame.
%
%   MilkyWaySource="galactic_band" instead draws an idealized constant-|b|
%   Galactic latitude band for geometric validation.
%

%% ========================================================================
%  INPUTS
% =========================================================================

arguments

    observerMode (1,1) string ...
        {mustBeMember(observerMode,["geocentric","heliocentric"])} ...
        = "geocentric"

    epochUTC (1,1) datetime = ...
        datetime('now','TimeZone','UTC')

    opts.DownloadMode (1,1) string ...
        {mustBeMember(opts.DownloadMode,["plottable","all"])} ...
        = "plottable"

    opts.Population (1,1) string ...
        {mustBeMember(opts.Population,["MBA","all_asteroids"])} ...
        = "MBA"

    % Local display filters. Inclusive bounds. Inf is allowed.
    opts.DiameterRange_km (1,2) double = [0 Inf]

    opts.SpinRateRange_revDay (1,2) double = [0 Inf]

    opts.CoordinateFrame (1,1) string ...
        {mustBeMember(opts.CoordinateFrame,["equatorial","ecliptic","galactic"])} ...
        = "equatorial"

    opts.RAUnits (1,1) string ...
        {mustBeMember(opts.RAUnits,["hours","degrees"])} ...
        = "hours"

    opts.LongitudeDirection (1,1) string ...
        {mustBeMember(opts.LongitudeDirection,["increasing_right","increasing_left"])} ...
        = "increasing_right"

    opts.Verbose (1,1) logical = true

    % When false, compute and return catalog/snapshot data but skip all
    % graphics construction. Used by the live visualizer to avoid building
    % a large temporary figure only to copy its contents.
    opts.RenderPlot (1,1) logical = true

    % Optional formation-normal sky-track overlay.  Supply either a scalar
    % struct or struct array produced by make_cr3bp_floquet_normal_overlays
% (the older make_halo_floquet_normal_overlay wrapper also remains valid).
    % Each overlay is expected to contain vectorsEclipticJ2000 (N x 3).
    opts.NormalOverlay = struct([])

    opts.ShowMilkyWay (1,1) logical = true

    opts.MilkyWaySource (1,1) string ...
        {mustBeMember(opts.MilkyWaySource,["legacy_footprint","legacy_polygon","galactic_band"])} ...
        = "legacy_footprint"

    opts.MilkyWayHalfWidth_deg (1,1) double ...
        {mustBePositive} = 12

    opts.ShowGalacticEquator (1,1) logical = true

    opts.ShowPlanets (1,1) logical = true

    opts.ShowPlanetLabels (1,1) logical = true

    opts.MilkyWayFill (1,1) logical = true

    opts.MilkyWayFillAlpha (1,1) double = 0.08

    opts.MilkyWayLineWidth (1,1) double ...
        {mustBePositive} = 0.9

    opts.MilkyWayGridStep_deg (1,1) double ...
        {mustBePositive} = 0.5

    opts.ForceCatalogRefresh (1,1) logical = false

    opts.MaxCacheAge_days (1,1) double ...
        {mustBePositive} = 7

    opts.CacheDirectory (1,1) string = ...
        string(fullfile(pwd,'JPL_SBDB_cache'))

    opts.PageSize (1,1) double ...
        {mustBeInteger,mustBePositive} = 10000

    opts.WriteCSV (1,1) logical = false

    opts.CenterLongitude_deg (1,1) double = 180

    opts.LongitudeGrid_deg (1,1) double ...
        {mustBePositive} = 30

    opts.LatitudeGrid_deg (1,1) double ...
        {mustBePositive} = 30

    opts.MaxMarkerArea_pt2 (1,1) double ...
        {mustBePositive} = 900

    opts.MinMarkerArea_pt2 (1,1) double ...
        {mustBeNonnegative} = 1.5

    opts.MarkerAlpha (1,1) double = 0.70

    opts.NLabels (1,1) double ...
        {mustBeInteger,mustBeNonnegative} = 12

    opts.ColorScale (1,1) string ...
        {mustBeMember(opts.ColorScale,["linear","log"])} ...
        = "linear"

    opts.ColorPercentiles (1,2) double = [0 100]

    opts.LightTimeIterations (1,1) double ...
        {mustBeInteger,mustBeNonnegative} = 3

end


if opts.MarkerAlpha < 0 || opts.MarkerAlpha > 1
    error('MarkerAlpha must lie between 0 and 1.');
end

if opts.MilkyWayFillAlpha < 0 || opts.MilkyWayFillAlpha > 1
    error('MilkyWayFillAlpha must lie between 0 and 1.');
end

if opts.MilkyWayHalfWidth_deg <= 0 || opts.MilkyWayHalfWidth_deg >= 90
    error('MilkyWayHalfWidth_deg must satisfy 0 < value < 90 deg.');
end

if any(isnan(opts.DiameterRange_km)) || ...
   opts.DiameterRange_km(1) < 0 || ...
   opts.DiameterRange_km(2) < opts.DiameterRange_km(1)

    error(['DiameterRange_km must be [Dmin Dmax] with ' ...
           '0 <= Dmin <= Dmax. Inf is allowed.']);

end

if any(isnan(opts.SpinRateRange_revDay)) || ...
   opts.SpinRateRange_revDay(1) < 0 || ...
   opts.SpinRateRange_revDay(2) < opts.SpinRateRange_revDay(1)

    error(['SpinRateRange_revDay must be [fmin fmax] with ' ...
           '0 <= fmin <= fmax. Inf is allowed.']);

end

if any(opts.ColorPercentiles < 0) || ...
   any(opts.ColorPercentiles > 100) || ...
   opts.ColorPercentiles(2) <= opts.ColorPercentiles(1)

    error('ColorPercentiles must satisfy 0 <= p1 < p2 <= 100.');

end


%% ========================================================================
%  POPULATION DEFINITION
% =========================================================================

switch opts.Population

    case "MBA"

        populationTag = "MBA";
        populationTitle = "Main-Belt Asteroids";

        populationArgs = { ...
            'sb-class','MBA' ...
            };

    case "all_asteroids"

        populationTag = "ALL_ASTEROIDS";
        populationTitle = "All Asteroids";

        % JPL SBDB:
        %   sb-kind = a   -> asteroids only
        %
        % No orbital-class restriction is imposed.
        populationArgs = { ...
            'sb-kind','a' ...
            };

end


%% ========================================================================
%  EPOCH
% =========================================================================

% An unzoned datetime is interpreted as UTC.
if isempty(epochUTC.TimeZone)
    epochUTC.TimeZone = 'UTC';
else
    epochUTC.TimeZone = 'UTC';
end

jdUTC = ...
    posixtime(epochUTC)/86400 + 2440587.5;

% SBDB orbital epochs use TDB.
%
% For the present-era application of this function, this conversion is
% much more than sufficient for an all-sky population map.
jdTDB = utc_to_jdtdb_approx(epochUTC);

epochPrint = epochUTC;
epochPrint.Format = 'yyyy-MM-dd HH:mm:ss.SSS';


vfprintf(opts.Verbose,'\n');
vfprintf(opts.Verbose,'============================================================\n');
vfprintf(opts.Verbose,' JPL ASTEROID SKY MAP\n');
vfprintf(opts.Verbose,'============================================================\n');
vfprintf(opts.Verbose,'Observer:       %s\n',upper(observerMode));
vfprintf(opts.Verbose,'Population:     %s\n',populationTitle);
vfprintf(opts.Verbose,'Coord. frame:   %s\n',opts.CoordinateFrame);
vfprintf(opts.Verbose,'RA labels:      %s\n',opts.RAUnits);
vfprintf(opts.Verbose,'Lon. direction: %s\n',opts.LongitudeDirection);
if opts.ShowMilkyWay
    vfprintf(opts.Verbose,'Milky Way:      %s',opts.MilkyWaySource);
    if opts.MilkyWaySource == "galactic_band"
        vfprintf(opts.Verbose,'  (|b| <= %.1f deg)',opts.MilkyWayHalfWidth_deg);
    end
    vfprintf(opts.Verbose,'\n');
end
vfprintf(opts.Verbose,'Epoch UTC:      %s\n',char(epochPrint));
vfprintf(opts.Verbose,'JD UTC:         %.9f\n',jdUTC);
vfprintf(opts.Verbose,'JD TDB approx:  %.9f\n',jdTDB);
vfprintf(opts.Verbose,'Download mode:  %s\n',opts.DownloadMode);
vfprintf(opts.Verbose,'Verbose:        %s\n',string(opts.Verbose));
vfprintf(opts.Verbose,'============================================================\n\n');


%% ========================================================================
%  CACHE DIRECTORY
% =========================================================================

cacheDirectory = char(opts.CacheDirectory);

if ~exist(cacheDirectory,'dir')
    mkdir(cacheDirectory);
end


cacheFile = fullfile( ...
    cacheDirectory, ...
    sprintf( ...
        'JPL_SBDB_%s_%s.mat', ...
        char(populationTag), ...
        char(opts.DownloadMode)));

csvFile = fullfile( ...
    cacheDirectory, ...
    sprintf( ...
        'JPL_SBDB_%s_%s.csv', ...
        char(populationTag), ...
        char(opts.DownloadMode)));


%% ========================================================================
%  DETERMINE WHETHER CATALOG MUST BE REFRESHED
% =========================================================================

refreshCatalog = ...
    opts.ForceCatalogRefresh || ...
    ~isfile(cacheFile);


if isfile(cacheFile) && ~opts.ForceCatalogRefresh

    d = dir(cacheFile);

    cacheAge_days = ...
        now - d.datenum; %#ok<TNOW1>

    if cacheAge_days > opts.MaxCacheAge_days
        refreshCatalog = true;
    end

end


%% ========================================================================
%  LOAD / DOWNLOAD SBDB
% =========================================================================

if refreshCatalog

    vfprintf(opts.Verbose,'Downloading NASA/JPL SBDB catalog...\n\n');

    catalog = download_JPL_catalog(opts,populationArgs,populationTag);

    catalogDownloadUTC = ...
        datetime('now','TimeZone','UTC');

    save( ...
        cacheFile, ...
        'catalog', ...
        'catalogDownloadUTC', ...
        '-v7.3');

    vfprintf(opts.Verbose,'\nSaved MAT cache:\n  %s\n',cacheFile);

    if opts.WriteCSV

        writetable(catalog,csvFile);

        vfprintf(opts.Verbose,'Saved CSV copy:\n  %s\n',csvFile);

    end

else

    vfprintf(opts.Verbose,'Using cached SBDB catalog:\n  %s\n\n',cacheFile);

    S = load(cacheFile,'catalog');

    catalog = S.catalog;

end


%% ========================================================================
%  IDENTIFY USABLE OBJECTS
% =========================================================================

validOrbit = ...
    isfinite(catalog.epoch) & ...
    isfinite(catalog.a) & ...
    catalog.a > 0 & ...
    isfinite(catalog.e) & ...
    catalog.e >= 0 & ...
    catalog.e < 1 & ...
    isfinite(catalog.i) & ...
    isfinite(catalog.om) & ...
    isfinite(catalog.w) & ...
    isfinite(catalog.ma);

validDiameter = ...
    isfinite(catalog.diameter_km) & ...
    catalog.diameter_km > 0;

validRotation = ...
    isfinite(catalog.rot_per_h) & ...
    catalog.rot_per_h > 0;


use = ...
    validOrbit & ...
    validDiameter & ...
    validRotation;


snapshot = catalog(use,:);

NusableBeforeFilter = height(snapshot);

% -------------------------------------------------------------------------
% LOCAL DISPLAY FILTERS
%
% Apply these before ephemeris propagation so restrictive filters also
% reduce the computational work. The cached/downloaded catalog itself is
% left untouched and can be reused with different filters on later calls.
% -------------------------------------------------------------------------

candidateSpinRate_revDay = ...
    24 ./ snapshot.rot_per_h;

diameterFilter = ...
    snapshot.diameter_km >= opts.DiameterRange_km(1) & ...
    snapshot.diameter_km <= opts.DiameterRange_km(2);

spinFilter = ...
    candidateSpinRate_revDay >= opts.SpinRateRange_revDay(1) & ...
    candidateSpinRate_revDay <= opts.SpinRateRange_revDay(2);

displayFilter = ...
    diameterFilter & spinFilter;

snapshot = ...
    snapshot(displayFilter,:);


vfprintf(opts.Verbose,'\n');
vfprintf(opts.Verbose,'------------------------------------------------------------\n');
vfprintf(opts.Verbose,' CATALOG CONTENT\n');
vfprintf(opts.Verbose,'------------------------------------------------------------\n');
vfprintf(opts.Verbose,'Downloaded/cached objects:           %8d\n',height(catalog));
vfprintf(opts.Verbose,'Valid orbital elements:              %8d\n',nnz(validOrbit));
vfprintf(opts.Verbose,'Measured diameter:                   %8d\n',nnz(validDiameter));
vfprintf(opts.Verbose,'Measured rotation period:            %8d\n',nnz(validRotation));
vfprintf(opts.Verbose,'Usable before display filters:       %8d\n',NusableBeforeFilter);
vfprintf(opts.Verbose,'Displayed after size/spin filters:    %8d\n',height(snapshot));

if opts.DiameterRange_km(1) > 0 || isfinite(opts.DiameterRange_km(2))
    vfprintf(opts.Verbose, ...
        'Diameter filter [km]:              [%g, %g]\n', ...
        opts.DiameterRange_km(1),opts.DiameterRange_km(2));
end

if opts.SpinRateRange_revDay(1) > 0 || isfinite(opts.SpinRateRange_revDay(2))
    vfprintf(opts.Verbose, ...
        'Spin-rate filter [rev/day]:        [%g, %g]\n', ...
        opts.SpinRateRange_revDay(1),opts.SpinRateRange_revDay(2));
end

vfprintf(opts.Verbose,'------------------------------------------------------------\n\n');


if isempty(snapshot)
    error(['No asteroids satisfy the requested population, physical-data, ' ...
           'diameter, and spin-rate filters.']);
end


%% ========================================================================
%  PROPAGATION SPAN DIAGNOSTIC
% =========================================================================

propagationAge_days = ...
    jdTDB - snapshot.epoch;

vfprintf(opts.Verbose,'Propagation from SBDB osculating epochs:\n');
vfprintf(opts.Verbose,'  Median |dt| : %.2f days\n', ...
    median(abs(propagationAge_days),'omitnan'));

vfprintf(opts.Verbose,'  Maximum |dt|: %.2f days\n\n', ...
    max(abs(propagationAge_days),[],'omitnan'));


%% ========================================================================
%  HORIZONS REFERENCE STATES
% =========================================================================

rEarth_ecl_AU = [];
planetTable = table();

switch observerMode

    case "geocentric"

        vfprintf(opts.Verbose, ...
            'Querying JPL Horizons for Earth heliocentric state...\n');

        [rEarth_ecl_AU,~] = ...
            get_body_heliocentric_state_Horizons("399",jdUTC);

        vfprintf(opts.Verbose, ...
            'Earth heliocentric r [AU] = [%+.8f  %+.8f  %+.8f]\n\n', ...
            rEarth_ecl_AU);


    case "heliocentric"

        if opts.ShowPlanets

            vfprintf(opts.Verbose, ...
                'Querying JPL Horizons for major-planet heliocentric states...\n');

            planetTable = ...
                get_major_planets_heliocentric_states_Horizons( ...
                    jdUTC, ...
                    opts.Verbose);

        end

end


%% ========================================================================
%  OBSERVER GEOMETRY
% =========================================================================

N = height(snapshot);

c_AU_day = ...
    299792.458 * 86400 / 149597870.700;


switch observerMode

    % =====================================================================
    case "geocentric"
    % =====================================================================

        vfprintf(opts.Verbose,'Computing GEOCENTRIC directions...\n');

        % Start with the reception epoch.
        tEmit_TDB = ...
            jdTDB .* ones(N,1);

        % ---------------------------------------------------------------
        % Iterative one-way light-time correction:
        %
        % t_emit = t_obs - |r_ast(t_emit) - r_E(t_obs)| / c
        % ---------------------------------------------------------------

        for iter = 1:opts.LightTimeIterations

            rAst_ecl_AU = ...
                propagate_asteroids_kepler( ...
                    snapshot, ...
                    tEmit_TDB);

            rho_ecl_AU = ...
                rAst_ecl_AU - rEarth_ecl_AU;

            range_AU = ...
                vecnorm(rho_ecl_AU,2,2);

            lightTime_days = ...
                range_AU ./ c_AU_day;

            tEmit_TDB = ...
                jdTDB - lightTime_days;

        end

        % Final target evaluation at the converged emission epoch.
        rAst_ecl_AU = ...
            propagate_asteroids_kepler( ...
                snapshot, ...
                tEmit_TDB);

        rho_ecl_AU = ...
            rAst_ecl_AU - rEarth_ecl_AU;

        range_AU = ...
            vecnorm(rho_ecl_AU,2,2);

        lightTime_min = ...
            range_AU ./ c_AU_day .* 1440;


    % =====================================================================
    case "heliocentric"
    % =====================================================================

        vfprintf(opts.Verbose,'Computing HELIOCENTRIC directions...\n');

        tNow_TDB = ...
            jdTDB .* ones(N,1);

        rAst_ecl_AU = ...
            propagate_asteroids_kepler( ...
                snapshot, ...
                tNow_TDB);

        % Observer is the Sun.
        rho_ecl_AU = ...
            rAst_ecl_AU;

        range_AU = ...
            vecnorm(rho_ecl_AU,2,2);

        lightTime_min = ...
            NaN(N,1);

end


%% ========================================================================
%  BOTH SKY COORDINATE SYSTEMS
%
%  SBDB/Horizons vectors are already J2000 ecliptic.  We keep both the
%  J2000 ecliptic longitude/latitude and J2000 equatorial RA/Dec so the
%  returned snapshot remains useful regardless of the plotted frame.
% =========================================================================

[eclLon_deg,eclLat_deg] = ...
    cartesian_to_lonlat(rho_ecl_AU);

rho_eq_AU = ...
    ecliptic_to_equatorial_J2000(rho_ecl_AU);

[RA_deg,Dec_deg] = ...
    cartesian_to_lonlat(rho_eq_AU);


snapshot.RA_deg = ...
    RA_deg;

snapshot.Dec_deg = ...
    Dec_deg;

snapshot.eclipticLongitude_deg = ...
    eclLon_deg;

snapshot.eclipticLatitude_deg = ...
    eclLat_deg;


% Galactic longitude/latitude are derived from the same J2000/ICRS
% equatorial direction vectors.  This is useful both as an output and as a
% direct validation frame for the Milky Way overlay.
rho_gal = ...
    equatorial_to_galactic_J2000(rho_eq_AU);

[galLon_deg,galLat_deg] = ...
    cartesian_to_lonlat(rho_gal);

snapshot.galacticLongitude_deg = ...
    galLon_deg;

snapshot.galacticLatitude_deg = ...
    galLat_deg;


snapshot.range_AU = ...
    range_AU;

snapshot.lightTime_min = ...
    lightTime_min;

snapshot.spinRate_revDay = ...
    24 ./ snapshot.rot_per_h;

snapshot.projectedArea_km2 = ...
    (pi/4) .* snapshot.diameter_km.^2;

snapshot.observationUTC = ...
    repmat(epochUTC,height(snapshot),1);


%% ========================================================================
%  SELECT PLOTTING FRAME
% =========================================================================

switch opts.CoordinateFrame

    case "equatorial"

        snapshot.plotLongitude_deg = ...
            snapshot.RA_deg;

        snapshot.plotLatitude_deg = ...
            snapshot.Dec_deg;

        coordinateFrameTitle = ...
            '\mathrm{Equatorial}\;(\mathrm{J2000}\;\alpha/\delta)';

    case "ecliptic"

        snapshot.plotLongitude_deg = ...
            snapshot.eclipticLongitude_deg;

        snapshot.plotLatitude_deg = ...
            snapshot.eclipticLatitude_deg;

        coordinateFrameTitle = ...
            '\mathrm{Ecliptic}\;(\mathrm{J2000}\;\lambda/\beta)';


    case "galactic"

        snapshot.plotLongitude_deg = ...
            snapshot.galacticLongitude_deg;

        snapshot.plotLatitude_deg = ...
            snapshot.galacticLatitude_deg;

        coordinateFrameTitle = ...
            '\mathrm{Galactic}\;(\mathrm{IAU}\;l/b)';

end


%% ========================================================================
%  REMOVE NUMERIC FAILURES
% =========================================================================

good = ...
    isfinite(snapshot.plotLongitude_deg) & ...
    isfinite(snapshot.plotLatitude_deg) & ...
    isfinite(snapshot.range_AU) & ...
    isfinite(snapshot.spinRate_revDay) & ...
    snapshot.spinRate_revDay > 0;


snapshot = ...
    snapshot(good,:);


vfprintf(opts.Verbose,'Objects surviving propagation: %d\n\n',height(snapshot));


%% ========================================================================
%  SORT BY DIAMETER
%
%  Largest first so smaller objects tend to remain visible on top.
% =========================================================================

[~,idxSort] = ...
    sort(snapshot.diameter_km,'descend');

snapshot = ...
    snapshot(idxSort,:);


%% ========================================================================
%  MOLLWEIDE COORDINATES
% =========================================================================

[xSky,ySky] = ...
    mollweide_lonlat( ...
        snapshot.plotLongitude_deg, ...
        snapshot.plotLatitude_deg, ...
        opts.CenterLongitude_deg, ...
        opts.LongitudeDirection);


%% ========================================================================
%  MARKER AREA
%
%  Physical cross-sectional area:
%
%       A_phys = pi D^2 / 4
%
%  MATLAB scatter() SizeData is graphical marker area.
%
%  Therefore:
%
%       SizeData ∝ A_phys
%
%  except for the explicitly defined visibility floor.
% =========================================================================

Aphys = ...
    snapshot.projectedArea_km2;

scaleArea = ...
    opts.MaxMarkerArea_pt2 / max(Aphys);

markerArea_pt2 = ...
    scaleArea .* Aphys;

markerArea_pt2 = ...
    max(markerArea_pt2,opts.MinMarkerArea_pt2);


%% ========================================================================
%  COLOR VARIABLE
% =========================================================================

spinRate = ...
    snapshot.spinRate_revDay;

% Cache plotting-ready values in the returned snapshot so callers such as
% the unified live visualizer can render directly without constructing an
% intermediate MATLAB figure.
snapshot.mollweideX = xSky;
snapshot.mollweideY = ySky;
snapshot.markerArea_pt2 = markerArea_pt2;

if ~opts.RenderPlot
    fig = gobjects(0);
    return
end


%% ========================================================================
%  FIGURE
% =========================================================================

fig = figure( ...
    'Theme','light', ...
    'Color','w', ...
    'Name',sprintf('JPL %s Sky',populationTitle), ...
    'Position',[80 80 1500 850]);


ax = axes(fig);
ax.Tag = 'JPL_AsteroidSkyAxes';

hold(ax,'on');
axis(ax,'equal');
axis(ax,'off');

set(fig,'Color','w');
set(ax,'Color','w','XColor','k','YColor','k');


%% ========================================================================
%  SKY GRID
% =========================================================================

draw_mollweide_grid( ...
    ax, ...
    opts.CenterLongitude_deg, ...
    opts.LongitudeGrid_deg, ...
    opts.LatitudeGrid_deg, ...
    opts.CoordinateFrame, ...
    opts.RAUnits, ...
    opts.LongitudeDirection);


%% ========================================================================
%  MILKY WAY OVERLAY
% =========================================================================

if opts.ShowMilkyWay

    draw_milkyway_mollweide( ...
        ax, ...
        opts.CenterLongitude_deg, ...
        opts.CoordinateFrame, ...
        opts.LongitudeDirection, ...
        opts.MilkyWaySource, ...
        opts.MilkyWayHalfWidth_deg, ...
        opts.ShowGalacticEquator, ...
        opts.MilkyWayFill, ...
        opts.MilkyWayFillAlpha, ...
        opts.MilkyWayLineWidth, ...
        opts.MilkyWayGridStep_deg);

end


%% ========================================================================
%  ECLIPTIC
% =========================================================================

draw_ecliptic( ...
    ax, ...
    opts.CenterLongitude_deg, ...
    opts.CoordinateFrame, ...
    opts.LongitudeDirection);


%% ========================================================================
%  ASTEROIDS
% =========================================================================

hs = scatter( ...
    ax, ...
    xSky, ...
    ySky, ...
    markerArea_pt2, ...
    spinRate, ...
    'filled', ...
    'Tag','JPL_AsteroidScatter', ...
    'MarkerFaceAlpha',opts.MarkerAlpha, ...
    'MarkerEdgeColor',[0.12 0.12 0.12], ...
    'MarkerEdgeAlpha',0.15);


colormap(ax,turbo(256));

ax.ColorScale = char(opts.ColorScale);


%% ========================================================================
%  COLOR LIMITS
% =========================================================================

cLo = ...
    percentile_simple( ...
        spinRate, ...
        opts.ColorPercentiles(1));

cHi = ...
    percentile_simple( ...
        spinRate, ...
        opts.ColorPercentiles(2));


if opts.ColorScale == "log"
    cLo = max(cLo,min(spinRate(spinRate > 0)));
end


if isfinite(cLo) && ...
   isfinite(cHi) && ...
   cHi > cLo

    clim(ax,[cLo cHi]);

end


cb = colorbar(ax);

cb.Label.String = ...
    '$24/P_{\mathrm{rot}}\;[\mathrm{rev\,day}^{-1}]$';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';
cb.Label.FontSize = 12;
cb.Color = 'k';
cb.Label.Color = 'k';


%% ========================================================================
%  SOLAR-SYSTEM REFERENCE OBJECTS
% =========================================================================

switch observerMode

    case "geocentric"

        % From Earth, the heliocentric Sun direction is -r_Earth.
        sun_ecl = ...
            -rEarth_ecl_AU;

        sun_plot = ...
            ecliptic_vectors_to_selected_frame( ...
                sun_ecl, ...
                opts.CoordinateFrame);

        [sunLon,sunLat] = ...
            cartesian_to_lonlat(sun_plot);

        [xSun,ySun] = ...
            mollweide_lonlat( ...
                sunLon, ...
                sunLat, ...
                opts.CenterLongitude_deg, ...
                opts.LongitudeDirection);

        plot( ...
            ax, ...
            xSun, ...
            ySun, ...
            'p', ...
            'MarkerSize',12, ...
            'MarkerFaceColor',[1.00 0.73 0.05], ...
            'MarkerEdgeColor',[0.15 0.15 0.15], ...
            'LineWidth',1.0, ...
            'HandleVisibility','off');

        text( ...
            ax, ...
            xSun + 0.04, ...
            ySun, ...
            '$\mathrm{Sun}$', ...
            'FontSize',9, ...
            'FontWeight','bold', ...
            'Color',[0.12 0.12 0.12], ...
            'Interpreter','latex', ...
            'Clipping','on');


    case "heliocentric"

        % The observer is at the Sun, so the Sun itself has no sky
        % direction. Plot all eight major planets instead.
        if opts.ShowPlanets && ~isempty(planetTable)

            plot_major_planets_mollweide( ...
                ax, ...
                planetTable, ...
                opts.CenterLongitude_deg, ...
                opts.CoordinateFrame, ...
                opts.LongitudeDirection, ...
                opts.ShowPlanetLabels);

        end

end


%% ========================================================================
%  LINEAR FORMATION-NORMAL OVERLAY
% =========================================================================

if ~isempty(opts.NormalOverlay)

    plot_normal_overlays_mollweide( ...
        ax, ...
        opts.NormalOverlay, ...
        opts.CenterLongitude_deg, ...
        opts.CoordinateFrame, ...
        opts.LongitudeDirection);

end


%% ========================================================================
%  LABEL LARGEST OBJECTS
% =========================================================================

nLabel = ...
    min(opts.NLabels,height(snapshot));


for k = 1:nLabel

    objectName = ...
        strtrim(snapshot.full_name(k));

    if ismissing(objectName) || ...
       strlength(objectName) == 0

        objectName = ...
            snapshot.pdes(k);

    end


    text( ...
        ax, ...
        xSky(k), ...
        ySky(k), ...
        "  " + objectName, ...
        'FontSize',8.3, ...
        'Color',[0.12 0.12 0.12], ...
        'Interpreter','none', ...
        'Clipping','on');

end


%% ========================================================================
%  DIAMETER / AREA LEGEND
% =========================================================================

Dexample = ...
    choose_diameter_legend(snapshot.diameter_km);


hSize = ...
    gobjects(numel(Dexample),1);


for k = 1:numel(Dexample)

    Aexample = ...
        pi/4 * Dexample(k)^2;

    Sexample = ...
        max( ...
            scaleArea*Aexample, ...
            opts.MinMarkerArea_pt2);

    hSize(k) = scatter( ...
        ax, ...
        NaN, ...
        NaN, ...
        Sexample, ...
        [0.48 0.48 0.48], ...
        'filled', ...
        'MarkerEdgeColor',[0.15 0.15 0.15]);

end


legendLabels = ...
    arrayfun(@(d) sprintf('$D = %g\\,\\mathrm{km}$',d), ...
        Dexample,'UniformOutput',false);


lgd = legend( ...
    hSize, ...
    legendLabels, ...
    'Location','southoutside', ...
    'Orientation','horizontal', ...
    'Interpreter','latex');

lgd.Box = 'off';
lgd.TextColor = 'k';
lgd.Color = 'w';


%% ========================================================================
%  TITLES
% =========================================================================

switch observerMode

    case "geocentric"

        viewTitle = ...
            'Geocentric instantaneous sky';

        viewDescription = ...
            'position = direction from Earth';


    case "heliocentric"

        viewTitle = ...
            'Heliocentric instantaneous sky';

        viewDescription = ...
            'position = direction from Sun';

end


epochShort = epochUTC;
epochShort.Format = 'yyyy-MM-dd HH:mm:ss';


filterLatex = "";

if opts.DiameterRange_km(1) > 0 || isfinite(opts.DiameterRange_km(2))
    filterLatex = filterLatex + ...
        "\quad|\quad D\in[" + ...
        string(opts.DiameterRange_km(1)) + "," + ...
        string(opts.DiameterRange_km(2)) + "]\,\mathrm{km}";
end

if opts.SpinRateRange_revDay(1) > 0 || isfinite(opts.SpinRateRange_revDay(2))
    filterLatex = filterLatex + ...
        "\quad|\quad f_{\mathrm{spin}}\in[" + ...
        string(opts.SpinRateRange_revDay(1)) + "," + ...
        string(opts.SpinRateRange_revDay(2)) + ...
        "]\,\mathrm{rev\,day}^{-1}";
end

switch observerMode
    case "geocentric"
        viewDescriptionLatex = "\mathrm{position=direction\;from\;Earth}";
    case "heliocentric"
        viewDescriptionLatex = "\mathrm{position=direction\;from\;Sun}";
end

line1 = "$" + latex_roman_text_local( ...
    "NASA/JPL " + string(populationTitle) + " -- " + string(viewTitle)) + "$";
line2 = "$" + latex_roman_text_local(string(epochShort) + " UTC") + "$";

overlaySummaryLatex = build_overlay_family_summary_latex(opts.NormalOverlay);

subtitleText = "$" + ...
    viewDescriptionLatex + ...
    "\quad|\quad \mathrm{frame=}\;" + string(coordinateFrameTitle) + ...
    "\quad|\quad A_{\mathrm{marker}}\propto A_{\mathrm{proj}}" + ...
    "\quad|\quad \mathrm{color=spin\;rate}" + ...
    "\quad|\quad N=" + string(height(snapshot)) + ...
    overlaySummaryLatex + ...
    filterLatex + "$";

hTitle = title( ...
    ax, ...
    {char(line1),char(line2)}, ...
    'FontSize',17, ...
    'FontWeight','bold', ...
    'Interpreter','latex');

hSubtitle = subtitle( ...
    ax, ...
    char(subtitleText), ...
    'FontSize',11, ...
    'Interpreter','latex');

hTitle.Color = 'k';
hSubtitle.Color = 'k';


%% ========================================================================
%  AXIS LIMITS
% =========================================================================

xMax = ...
    2*sqrt(2);

yMax = ...
    sqrt(2);


xlim(ax,1.07*[-xMax xMax]);
ylim(ax,1.16*[-yMax yMax]);


%% ========================================================================
%  OPTIONAL INTERACTIVE DATA TIPS
% =========================================================================

try

    hs.DataTipTemplate.DataTipRows = [ ...
        dataTipTextRow( ...
            'Object', ...
            snapshot.full_name)
        dataTipTextRow( ...
            'RA [deg]', ...
            snapshot.RA_deg)
        dataTipTextRow( ...
            'Dec [deg]', ...
            snapshot.Dec_deg)
        dataTipTextRow( ...
            'Ecl. lon [deg]', ...
            snapshot.eclipticLongitude_deg)
        dataTipTextRow( ...
            'Ecl. lat [deg]', ...
            snapshot.eclipticLatitude_deg)
        dataTipTextRow( ...
            'Diameter [km]', ...
            snapshot.diameter_km)
        dataTipTextRow( ...
            'Spin [rev/day]', ...
            snapshot.spinRate_revDay)
        dataTipTextRow( ...
            'Range [AU]', ...
            snapshot.range_AU)
        ];

catch
    % DataTipTemplate behavior varies slightly between MATLAB releases.
    % The plot itself does not depend on this feature.
end


%% ========================================================================
%  FINAL SUMMARY
% =========================================================================

vfprintf(opts.Verbose,'\n');
vfprintf(opts.Verbose,'============================================================\n');
vfprintf(opts.Verbose,' SKY MAP COMPLETE\n');
vfprintf(opts.Verbose,'============================================================\n');
vfprintf(opts.Verbose,'Observer:        %s\n',observerMode);
vfprintf(opts.Verbose,'Population:      %s\n',populationTitle);
vfprintf(opts.Verbose,'Coordinate frame: %s\n',opts.CoordinateFrame);
vfprintf(opts.Verbose,'RA label units:   %s\n',opts.RAUnits);
vfprintf(opts.Verbose,'Objects plotted: %d\n',height(snapshot));
vfprintf(opts.Verbose,'Epoch:           %s UTC\n',char(epochShort));
vfprintf(opts.Verbose,'============================================================\n\n');


end


%% =========================================================================
%
%  DOWNLOAD JPL SBDB CATALOG
%
% =========================================================================

function T = download_JPL_catalog(opts,populationArgs,populationTag)

baseURL = ...
    'https://ssd-api.jpl.nasa.gov/sbdb_query.api';


% -------------------------------------------------------------------------
% Fields retained locally
% -------------------------------------------------------------------------

fields = [ ...
    "spkid"
    "pdes"
    "full_name"
    "class"
    "epoch"
    "equinox"
    "e"
    "a"
    "i"
    "om"
    "w"
    "ma"
    "n"
    "H"
    "diameter"
    "diameter_sigma"
    "rot_per"
    "albedo"
    ];


fieldString = ...
    char(strjoin(fields,','));


% -------------------------------------------------------------------------
% Optional physical-data constraint
% -------------------------------------------------------------------------

switch opts.DownloadMode

    case "plottable"

        % Only objects for which both quantities required for the plot are
        % actually defined.
        filterJSON = ...
            '{"AND":["diameter|DF","rot_per|DF"]}';


    case "all"

        filterJSON = ...
            '';

end


webOpts = ...
    weboptions( ...
        'Timeout',180, ...
        'ContentType','json');


%% ========================================================================
%  LIVE OBJECT COUNT
% =========================================================================

countArgs = populationArgs;


if ~isempty(filterJSON)

    countArgs = [ ...
        countArgs, ...
        {'sb-cdata',filterJSON} ...
        ];

end


countPayload = ...
    webread_retry( ...
        baseURL, ...
        countArgs, ...
        webOpts);


Nexpected = ...
    str2double(string(countPayload.count));


vfprintf(opts.Verbose,'JPL query:\n');
vfprintf(opts.Verbose,'  population: %s\n',populationTag);
vfprintf(opts.Verbose,'  mode:       %s\n',opts.DownloadMode);

if opts.DownloadMode == "plottable"
    vfprintf(opts.Verbose,'  diameter:   defined\n');
    vfprintf(opts.Verbose,'  rot_per:    defined\n');
end

vfprintf(opts.Verbose,'\n');
vfprintf(opts.Verbose,'JPL reports %d matching objects.\n\n',Nexpected);


%% ========================================================================
%  PAGINATED DOWNLOAD
% =========================================================================

nPages = ...
    ceil(Nexpected/opts.PageSize);


chunks = ...
    cell(nPages,1);


for page = 1:nPages

    offset = ...
        (page-1)*opts.PageSize;


    vfprintf(opts.Verbose, ...
        'Page %3d / %3d : records %8d -- %8d\n', ...
        page, ...
        nPages, ...
        offset+1, ...
        min(offset+opts.PageSize,Nexpected));


    args = [ ...
        { ...
        'fields',fieldString, ...
        'sort','spkid', ...
        'limit',opts.PageSize, ...
        'limit-from',offset, ...
        'full-prec','true' ...
        }, ...
        populationArgs ...
        ];


    if ~isempty(filterJSON)

        args = [ ...
            args, ...
            {'sb-cdata',filterJSON} ...
            ];

    end


    payload = ...
        webread_retry( ...
            baseURL, ...
            args, ...
            webOpts);


    raw = ...
        normalize_JPL_data( ...
            payload.data, ...
            numel(fields));


    chunks{page} = ...
        convert_JPL_chunk( ...
            raw, ...
            fields);

end


%% ========================================================================
%  MERGE
% =========================================================================

T = ...
    vertcat(chunks{:});


% SBDB documentation notes that the underlying database can change between
% paginated requests. Remove a duplicated SPK-ID if an update happens during
% a long download.

if ~isempty(T)

    [~,idxUnique] = ...
        unique(T.spkid,'stable');

    if numel(idxUnique) ~= height(T)

        warning( ...
            ['Duplicate SPK-IDs detected during paginated download. ' ...
             'Removing duplicates.']);

        T = ...
            T(idxUnique,:);

    end

end


vfprintf(opts.Verbose,'\n');
vfprintf(opts.Verbose,'Downloaded %d unique objects.\n',height(T));


if height(T) ~= Nexpected

    warning( ...
        ['Initial SBDB count was %d but %d unique rows were downloaded. ' ...
         'The live database may have changed during paging.'], ...
        Nexpected, ...
        height(T));

end

end


%% =========================================================================
%
%  WEBREAD WITH RETRIES
%
% =========================================================================

function payload = webread_retry(url,args,webOpts)

maxAttempts = ...
    5;


for attempt = 1:maxAttempts

    try

        payload = ...
            webread( ...
                url, ...
                args{:}, ...
                webOpts);

        return


    catch ME

        if attempt == maxAttempts
            rethrow(ME);
        end


        delay_s = ...
            2^(attempt-1);


        warning( ...
            'JPL request failed. Retrying in %g s.', ...
            delay_s);


        pause(delay_s);

    end

end

end


%% =========================================================================
%
%  NORMALIZE JPL JSON DATA
%
% =========================================================================

function C = normalize_JPL_data(data,nFields)

if isempty(data)

    C = ...
        cell(0,nFields);

    return

end


% Most MATLAB releases produce an N x nFields cell matrix.
if iscell(data) && ...
   ismatrix(data) && ...
   size(data,2) == nFields

    C = data;

    return

end


% Some releases may decode it as a cell vector of cell rows.
if iscell(data) && isvector(data)

    N = ...
        numel(data);

    C = ...
        cell(N,nFields);


    for k = 1:N

        row = ...
            data{k};

        if ~iscell(row)
            error('Unexpected SBDB JSON row format.');
        end

        row = ...
            row(:).';

        nCopy = ...
            min(numel(row),nFields);

        C(k,1:nCopy) = ...
            row(1:nCopy);

    end

    return

end


error('Unable to interpret the JPL SBDB JSON response.');

end


%% =========================================================================
%
%  CONVERT JPL PAGE TO MATLAB TABLE
%
% =========================================================================

function T = convert_JPL_chunk(C,fields)

N = ...
    size(C,1);


T = ...
    table();


for j = 1:numel(fields)

    field = ...
        fields(j);

    column = ...
        C(:,j);


    switch field

        % ---------------------------------------------------------------
        % Numeric fields
        % ---------------------------------------------------------------
        case { ...
                "epoch", ...
                "e", ...
                "a", ...
                "i", ...
                "om", ...
                "w", ...
                "ma", ...
                "n", ...
                "H", ...
                "diameter", ...
                "diameter_sigma", ...
                "rot_per", ...
                "albedo"}

            value = ...
                NaN(N,1);


            for k = 1:N

                value(k) = ...
                    JPL_to_double(column{k});

            end


            switch field

                case "diameter"

                    T.diameter_km = ...
                        value;


                case "diameter_sigma"

                    T.diameter_sigma_km = ...
                        value;


                case "rot_per"

                    T.rot_per_h = ...
                        value;


                otherwise

                    T.(field) = ...
                        value;

            end


        % ---------------------------------------------------------------
        % String fields
        % ---------------------------------------------------------------
        otherwise

            value = ...
                strings(N,1);


            for k = 1:N

                value(k) = ...
                    JPL_to_string(column{k});

            end


            T.(field) = ...
                value;

    end

end

end


%% =========================================================================
%
%  JPL VALUE -> DOUBLE
%
% =========================================================================

function x = JPL_to_double(v)

if isempty(v)

    x = NaN;


elseif isnumeric(v)

    x = double(v);


elseif islogical(v)

    x = double(v);


elseif ischar(v) || isstring(v)

    x = str2double(string(v));


else

    x = NaN;

end

end


%% =========================================================================
%
%  JPL VALUE -> STRING
%
% =========================================================================

function s = JPL_to_string(v)

if isempty(v)

    s = missing;


elseif ischar(v) || isstring(v)

    s = string(v);


elseif isnumeric(v) || islogical(v)

    s = string(v);


else

    s = missing;

end

end


%% =========================================================================
%
%  GENERIC HELIOCENTRIC BODY STATE FROM JPL HORIZONS
%
%  commandID examples:
%       "199" Mercury
%       "299" Venus
%       "399" Earth
%       "499" Mars
%       "599" Jupiter
%       "699" Saturn
%       "799" Uranus
%       "899" Neptune
%
%  Output frame:
%       Sun-centered J2000 ecliptic, AU and AU/day.
%
% =========================================================================

function [rBody,vBody] = ...
    get_body_heliocentric_state_Horizons(commandID,jdUTC)

baseURL = ...
    'https://ssd.jpl.nasa.gov/api/horizons.api';


webOpts = ...
    weboptions( ...
        'Timeout',120, ...
        'ContentType','json');


commandQuoted = ...
    sprintf('''%s''',char(commandID));


args = { ...
    'format','json', ...
    'COMMAND',commandQuoted, ...
    'OBJ_DATA','''NO''', ...
    'MAKE_EPHEM','''YES''', ...
    'EPHEM_TYPE','''VECTORS''', ...
    'CENTER','''500@10''', ...
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


payload = ...
    webread_retry( ...
        baseURL, ...
        args, ...
        webOpts);


txt = ...
    char(payload.result);


iSOE = ...
    strfind(txt,'$$SOE');

iEOE = ...
    strfind(txt,'$$EOE');


if isempty(iSOE) || isempty(iEOE)

    error( ...
        'Could not locate Horizons $$SOE/$$EOE ephemeris block for command %s.', ...
        commandID);

end


block = ...
    txt(iSOE(1):iEOE(1));


X  = horizons_number(block,'X');
Y  = horizons_number(block,'Y');
Z  = horizons_number(block,'Z');

VX = horizons_number(block,'VX');
VY = horizons_number(block,'VY');
VZ = horizons_number(block,'VZ');


rBody = ...
    [X Y Z];

vBody = ...
    [VX VY VZ];


if any(~isfinite(rBody)) || any(~isfinite(vBody))

    error( ...
        ['Horizons returned a non-finite state for command %s.\n' ...
         'r = [%g %g %g] AU\n' ...
         'v = [%g %g %g] AU/day'], ...
        commandID, ...
        rBody(1), ...
        rBody(2), ...
        rBody(3), ...
        vBody(1), ...
        vBody(2), ...
        vBody(3));

end

end


%% =========================================================================
%
%  MAJOR PLANET HELIOCENTRIC STATES
%
% =========================================================================

function planetTable = ...
    get_major_planets_heliocentric_states_Horizons(jdUTC,verbose)

names = [ ...
    "Mercury"
    "Venus"
    "Earth"
    "Mars"
    "Jupiter"
    "Saturn"
    "Uranus"
    "Neptune"
    ];

commands = [ ...
    "199"
    "299"
    "399"
    "499"
    "599"
    "699"
    "799"
    "899"
    ];


N = numel(names);

R = NaN(N,3);
V = NaN(N,3);


for k = 1:N

    vfprintf( ...
        verbose, ...
        '  Horizons: %-8s ... ', ...
        names(k));

    [R(k,:),V(k,:)] = ...
        get_body_heliocentric_state_Horizons( ...
            commands(k), ...
            jdUTC);

    vfprintf(verbose,'done\n');

end


planetTable = table( ...
    names, ...
    commands, ...
    R(:,1), ...
    R(:,2), ...
    R(:,3), ...
    V(:,1), ...
    V(:,2), ...
    V(:,3), ...
    'VariableNames',{ ...
        'Name', ...
        'Command', ...
        'x_AU', ...
        'y_AU', ...
        'z_AU', ...
        'vx_AUday', ...
        'vy_AUday', ...
        'vz_AUday'});


vfprintf(verbose,'\n');

end


%% =========================================================================
%
%  PARSE ONE HORIZONS VECTOR VALUE
%
%  Horizons commonly writes several quantities on the same line, e.g.
%
%      X = ...  Y = ...  Z = ...
%
%  so the label must NOT be required to occur at the beginning of a line.
%
% =========================================================================

function value = horizons_number(block,label)

% Robust floating-point pattern.
%
% Handles:
%   1.234
%   -1.234
%   +1.234E-03
%   1.234D+02
%
numberPattern = ...
    '([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+\-]?\d+)?)';


% Require either whitespace/comma or beginning of string before the label.
%
% This is important because searching for "Y" must NOT accidentally match
% the Y inside "VY".
pattern = ...
    sprintf( ...
        '(?:^|[\\s,])%s\\s*=\\s*%s', ...
        regexptranslate('escape',label), ...
        numberPattern);


token = ...
    regexp( ...
        block, ...
        pattern, ...
        'tokens', ...
        'once');


if isempty(token)

    fprintf('\nHorizons block returned:\n');
    fprintf('------------------------------------------------------------\n');
    fprintf('%s\n',block);
    fprintf('------------------------------------------------------------\n\n');

    error( ...
        'Could not parse Horizons quantity "%s".', ...
        label);

end


txt = token{1};

% Horizons can in principle use Fortran-style D exponents.
txt = ...
    regexprep(txt,'[Dd]','E');


value = ...
    str2double(txt);

end
%% =========================================================================
%
%  VECTORIZED TWO-BODY KEPLER PROPAGATION
%
% =========================================================================

function r = propagate_asteroids_kepler(T,jd)

a = ...
    T.a;

e = ...
    T.e;


inc = ...
    deg2rad(T.i);

Omega = ...
    deg2rad(T.om);

omega = ...
    deg2rad(T.w);

M0 = ...
    deg2rad(T.ma);


%% ========================================================================
%  MEAN MOTION
% =========================================================================

n = ...
    deg2rad(T.n);


badN = ...
    ~isfinite(n) | ...
    n <= 0;


% Gaussian gravitational constant:
%
% sqrt(mu_sun) in AU^(3/2)/day
kGaussian = ...
    0.01720209895;


n(badN) = ...
    kGaussian ./ ...
    a(badN).^(3/2);


%% ========================================================================
%  MEAN ANOMALY
% =========================================================================

dt = ...
    jd - T.epoch;


M = ...
    M0 + n.*dt;


M = ...
    mod(M,2*pi);


%% ========================================================================
%  SOLVE KEPLER EQUATION
%
%      E - e sin E = M
% =========================================================================

E = ...
    M;


for iter = 1:15

    f = ...
        E - e.*sin(E) - M;

    fp = ...
        1 - e.*cos(E);

    dE = ...
        -f./fp;

    E = ...
        E + dE;


    if max(abs(dE),[],'omitnan') < 1e-13
        break
    end

end


%% ========================================================================
%  PERIFOCAL COORDINATES
% =========================================================================

xp = ...
    a .* ...
    (cos(E) - e);


yp = ...
    a .* ...
    sqrt(1-e.^2) .* ...
    sin(E);


%% ========================================================================
%  ROTATION INTO J2000 ECLIPTIC FRAME
% =========================================================================

cO = cos(Omega);
sO = sin(Omega);

co = cos(omega);
so = sin(omega);

ci = cos(inc);
si = sin(inc);


R11 = ...
    cO.*co - ...
    sO.*so.*ci;

R12 = ...
    -cO.*so - ...
    sO.*co.*ci;


R21 = ...
    sO.*co + ...
    cO.*so.*ci;

R22 = ...
    -sO.*so + ...
    cO.*co.*ci;


R31 = ...
    so.*si;

R32 = ...
    co.*si;


x = ...
    R11.*xp + ...
    R12.*yp;

y = ...
    R21.*xp + ...
    R22.*yp;

z = ...
    R31.*xp + ...
    R32.*yp;


r = ...
    [x y z];

end


%% =========================================================================
%
%  J2000 ECLIPTIC -> J2000 EQUATORIAL / ICRF
%
% =========================================================================

function rEq = ecliptic_to_equatorial_J2000(rEcl)

% IAU76/80 J2000 obliquity.
epsJ2000 = ...
    deg2rad(84381.448/3600);

x = rEcl(:,1);
y = rEcl(:,2);
z = rEcl(:,3);

xEq = x;

yEq = ...
    cos(epsJ2000).*y - ...
    sin(epsJ2000).*z;

zEq = ...
    sin(epsJ2000).*y + ...
    cos(epsJ2000).*z;

rEq = ...
    [xEq yEq zEq];

end


%% =========================================================================
%
%  J2000 EQUATORIAL / ICRF -> J2000 ECLIPTIC
%
% =========================================================================

function rEcl = equatorial_to_ecliptic_J2000(rEq)

epsJ2000 = ...
    deg2rad(84381.448/3600);

x = rEq(:,1);
y = rEq(:,2);
z = rEq(:,3);

xEcl = x;

yEcl = ...
    cos(epsJ2000).*y + ...
    sin(epsJ2000).*z;

zEcl = ...
   -sin(epsJ2000).*y + ...
    cos(epsJ2000).*z;

rEcl = ...
    [xEcl yEcl zEcl];

end


%% =========================================================================
%
%  CONVERT J2000 ECLIPTIC VECTORS TO THE SELECTED SKY FRAME
%
% =========================================================================

function rPlot = ...
    ecliptic_vectors_to_selected_frame(rEcl,coordinateFrame)

switch coordinateFrame

    case "equatorial"
        rPlot = ...
            ecliptic_to_equatorial_J2000(rEcl);

    case "ecliptic"
        rPlot = ...
            rEcl;

    case "galactic"
        rEq = ...
            ecliptic_to_equatorial_J2000(rEcl);

        rPlot = ...
            equatorial_to_galactic_J2000(rEq);

end

end


%% =========================================================================
%
%  CARTESIAN -> LONGITUDE / LATITUDE
%
% =========================================================================

function [lon_deg,lat_deg] = cartesian_to_lonlat(r)

x = r(:,1);
y = r(:,2);
z = r(:,3);

rho = ...
    vecnorm(r,2,2);

lon = ...
    atan2(y,x);

lon = ...
    mod(lon,2*pi);

lat = ...
    asin(z./rho);

lon_deg = ...
    rad2deg(lon);

lat_deg = ...
    rad2deg(lat);

end


%% =========================================================================
%
%  LONGITUDE / LATITUDE -> UNIT CARTESIAN VECTOR
%
% =========================================================================

function r = lonlat_to_cartesian(lon_deg,lat_deg)

lon = ...
    deg2rad(lon_deg);

lat = ...
    deg2rad(lat_deg);

r = [ ...
    cos(lat).*cos(lon), ...
    cos(lat).*sin(lon), ...
    sin(lat) ...
    ];

end


%% =========================================================================
%
%  MOLLWEIDE PROJECTION
%
%  Longitude display direction is selectable. The input longitude may be
%  Right Ascension or ecliptic longitude.
%
% =========================================================================

function [x,y] = ...
    mollweide_lonlat(lon_deg,lat_deg,centerLongitude_deg,longitudeDirection)

lonRelative_deg = ...
    mod( ...
        lon_deg - centerLongitude_deg + 180, ...
        360) - 180;

[x,y] = ...
    mollweide_relative_lonlat( ...
        lonRelative_deg, ...
        lat_deg, ...
        longitudeDirection);

end


%% =========================================================================
%
%  MOLLWEIDE FROM ALREADY-CENTERED RELATIVE LONGITUDE
%
% =========================================================================

function [x,y] = ...
    mollweide_relative_lonlat(lonRelative_deg,lat_deg,longitudeDirection)

switch longitudeDirection

    case "increasing_right"
        directionSign = +1;

    case "increasing_left"
        directionSign = -1;

end

lon = ...
    directionSign .* deg2rad(lonRelative_deg);

lat = ...
    deg2rad(lat_deg);


theta = ...
    lat;

poleMask = ...
    abs(abs(lat)-pi/2) < 1e-12;

theta(poleMask) = ...
    sign(lat(poleMask))*pi/2;


for iter = 1:15

    good = ...
        ~poleMask;

    if ~any(good)
        break
    end

    f = ...
        2*theta(good) + ...
        sin(2*theta(good)) - ...
        pi*sin(lat(good));

    df = ...
        2 + ...
        2*cos(2*theta(good));

    theta(good) = ...
        theta(good) - f./df;

end


x = ...
    2*sqrt(2)/pi .* ...
    lon .* ...
    cos(theta);

y = ...
    sqrt(2) .* ...
    sin(theta);

end


%% =========================================================================
%
%  DRAW MOLLWEIDE GRID
%
% =========================================================================

function draw_mollweide_grid( ...
    ax,centerLongitude_deg,dLongitude,dLatitude,coordinateFrame,RAUnits,longitudeDirection)

gridColor = ...
    [0.78 0.78 0.78];

edgeColor = ...
    [0.22 0.22 0.22];


%% ========================================================================
%  BOUNDARY
% =========================================================================

q = ...
    linspace(0,2*pi,1200);

plot( ...
    ax, ...
    2*sqrt(2)*cos(q), ...
    sqrt(2)*sin(q), ...
    'Color',edgeColor, ...
    'LineWidth',1.1, ...
    'HandleVisibility','off');


%% ========================================================================
%  LATITUDE GRID
% =========================================================================

lonFine = ...
    linspace(0,360,1400);

for lat = (-90+dLatitude):dLatitude:(90-dLatitude)

    [x,y] = ...
        mollweide_lonlat( ...
            lonFine, ...
            lat*ones(size(lonFine)), ...
            centerLongitude_deg, ...
            longitudeDirection);

    [x,y] = ...
        break_map_seam(x,y);

    plot( ...
        ax, ...
        x, ...
        y, ...
        'Color',gridColor, ...
        'LineWidth',0.60, ...
        'HandleVisibility','off');

end


%% ========================================================================
%  LONGITUDE GRID
% =========================================================================

latFine = ...
    linspace(-89.9,89.9,700);

lonValues = ...
    0:dLongitude:(360-dLongitude);

for lon = lonValues

    [x,y] = ...
        mollweide_lonlat( ...
            lon*ones(size(latFine)), ...
            latFine, ...
            centerLongitude_deg, ...
            longitudeDirection);

    plot( ...
        ax, ...
        x, ...
        y, ...
        'Color',gridColor, ...
        'LineWidth',0.60, ...
        'HandleVisibility','off');

end


%% ========================================================================
%  LONGITUDE LABELS
% =========================================================================

for lon = lonValues

    [x0,~] = ...
        mollweide_lonlat( ...
            lon, ...
            0, ...
            centerLongitude_deg, ...
            longitudeDirection);

    if abs(x0) > 2.76
        continue
    end


    switch coordinateFrame

        case "equatorial"

            switch RAUnits

                case "hours"
                    lonLabel = ...
                        sprintf('$%g^{\\mathrm{h}}$',mod(lon/15,24));

                case "degrees"
                    lonLabel = ...
                        sprintf('$%g^{\\circ}$',mod(lon,360));

            end

        case "ecliptic"

            lonLabel = ...
                sprintf('$%g^{\\circ}$',mod(lon,360));


        case "galactic"

            lonLabel = ...
                sprintf('$%g^{\\circ}$',mod(lon,360));

    end


    text( ...
        ax, ...
        x0, ...
        -1.49, ...
        lonLabel, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize',9, ...
        'Color',[0.33 0.33 0.33], ...
        'Interpreter','latex');

end


%% ========================================================================
%  LATITUDE LABELS
% =========================================================================

for lat = (-90+dLatitude):dLatitude:(90-dLatitude)

    [x0,y0] = ...
        mollweide_lonlat( ...
            centerLongitude_deg, ...
            lat, ...
            centerLongitude_deg, ...
            longitudeDirection);

    text( ...
        ax, ...
        x0 + 0.045, ...
        y0, ...
        sprintf('$%+d^{\\circ}$',lat), ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','middle', ...
        'FontSize',8, ...
        'Color',[0.38 0.38 0.38], ...
        'Interpreter','latex');

end


%% ========================================================================
%  AXIS LABELS
% =========================================================================

switch coordinateFrame

    case "equatorial"

        switch RAUnits

            case "hours"
                xLabel = '$\alpha\;[\mathrm{h}]$';

            case "degrees"
                xLabel = '$\alpha\;[^{\circ}]$';

        end

        yLabel = '$\delta\;[^{\circ}]$';

    case "ecliptic"
        xLabel = '$\lambda\;[^{\circ}]$';
        yLabel = '$\beta\;[^{\circ}]$';

    case "galactic"
        xLabel = '$l\;[^{\circ}]$';
        yLabel = '$b\;[^{\circ}]$';

end


text( ...
    ax, ...
    0, ...
    -1.68, ...
    xLabel, ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold', ...
    'FontSize',10, ...
    'Interpreter','latex');


text( ...
    ax, ...
    -3.08, ...
    0, ...
    yLabel, ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold', ...
    'FontSize',10, ...
    'Rotation',90, ...
    'Interpreter','latex');

end


%% =========================================================================
%
%  DRAW MILKY WAY / GALACTIC-PLANE REGION
%
%  Two representations are supported:
%
%    "legacy_footprint" (default)
%       Uses the historical irregular Milky Way footprint supplied by
%       milkyway.m in J2000 equatorial RA/Dec.  Its original vertex order
%       is preserved.  The boundary is first densified in the ORIGINAL
%       RA/Dec chart coordinates, then each point is transformed as a 3-D
%       celestial unit vector into the requested frame, and only then is
%       the Mollweide projection applied.
%
%       The fill is computed independently: a regular grid in the OUTPUT
%       celestial frame is inverse-transformed to equatorial RA/Dec and
%       tested against the ORIGINAL RA/Dec footprint.  This avoids trying
%       to rotate a planar polyshape and avoids false seam-closing chords.
%
%    "galactic_band"
%       Draws the idealized spherical band |b| <= MilkyWayHalfWidth_deg.
%       This is useful as a frame-validation reference but is not the
%       historical irregular visible-Milky-Way footprint.
%
%  The optional Galactic-equator centerline b=0 is generated independently
%  from the standard ICRS/J2000 <-> Galactic rotation.  Therefore it is a
%  useful check on the transformed legacy footprint.
%
% =========================================================================

function h = draw_milkyway_mollweide( ...
    ax,centerLongitude_deg,coordinateFrame,longitudeDirection, ...
    source,halfWidth_deg,showGalacticEquator, ...
    doFill,fillAlpha,lineWidth,gridStep_deg)

switch source

    case "galactic_band"

        h = draw_galactic_band_mollweide( ...
            ax, ...
            centerLongitude_deg, ...
            coordinateFrame, ...
            longitudeDirection, ...
            halfWidth_deg, ...
            showGalacticEquator, ...
            doFill, ...
            fillAlpha, ...
            lineWidth, ...
            gridStep_deg);


    case {"legacy_footprint","legacy_polygon"}

        if exist('milkyway','file') ~= 2

            warning( ...
                ['MilkyWaySource="legacy_footprint" requires milkyway.m ' ...
                 'on the MATLAB path. Continuing without the Milky Way overlay.']);

            h = struct( ...
                'fill',gobjects(0), ...
                'contour',gobjects(0), ...
                'centerline',gobjects(0));

            return
        end

        h = draw_legacy_milkyway_mollweide( ...
            ax, ...
            centerLongitude_deg, ...
            coordinateFrame, ...
            longitudeDirection, ...
            showGalacticEquator, ...
            doFill, ...
            fillAlpha, ...
            lineWidth, ...
            gridStep_deg);

end

end


%% =========================================================================
%
%  PHYSICALLY DEFINED CONSTANT-GALACTIC-LATITUDE BAND
%
% =========================================================================

function h = draw_galactic_band_mollweide( ...
    ax,centerLongitude_deg,coordinateFrame,longitudeDirection, ...
    halfWidth_deg,showGalacticEquator, ...
    doFill,fillAlpha,lineWidth,gridStep_deg)

h = struct( ...
    'fill',gobjects(0), ...
    'contour',gobjects(0), ...
    'centerline',gobjects(0));


l_deg = (0:gridStep_deg:360).';

if l_deg(end) < 360
    l_deg(end+1,1) = 360;
end


if doFill

    nLat = max(3,ceil(2*halfWidth_deg/gridStep_deg)+1);
    b_deg = linspace(-halfWidth_deg,+halfWidth_deg,nLat).';

    [L,B] = meshgrid(l_deg,b_deg);

    [lonPlot,latPlot] = ...
        galactic_lonlat_to_selected_sky(L,B,coordinateFrame);

    [xFill,yFill] = ...
        mollweide_lonlat( ...
            lonPlot,latPlot,centerLongitude_deg,longitudeDirection);

    [xFill,yFill] = break_surface_map_seam(xFill,yFill);

    h.fill = surface( ...
        ax,xFill,yFill,zeros(size(xFill)),zeros(size(xFill)), ...
        'FaceColor',[0.68 0.68 0.68], ...
        'FaceAlpha',fillAlpha, ...
        'EdgeColor','none', ...
        'HandleVisibility','off');

    view(ax,2);

end


boundaryB = [+halfWidth_deg,-halfWidth_deg];
hBoundary = gobjects(2,1);

for k = 1:2

    [lonPlot,latPlot] = ...
        galactic_lonlat_to_selected_sky( ...
            l_deg,boundaryB(k).*ones(size(l_deg)),coordinateFrame);

    [x,y] = project_sky_curve_mollweide( ...
        lonPlot,latPlot,centerLongitude_deg,longitudeDirection);

    hBoundary(k) = plot( ...
        ax,x,y,'-', ...
        'Color',[0.42 0.42 0.42], ...
        'LineWidth',lineWidth, ...
        'HandleVisibility','off');

end

h.contour = hBoundary;


if showGalacticEquator

    h.centerline = draw_galactic_equator( ...
        ax,centerLongitude_deg,coordinateFrame,longitudeDirection,lineWidth);

end

end


%% =========================================================================
%
%  LEGACY IRREGULAR MILKY WAY FOOTPRINT
%
% =========================================================================

function h = draw_legacy_milkyway_mollweide( ...
    ax,centerLongitude_deg,coordinateFrame,longitudeDirection, ...
    showGalacticEquator,doFill,fillAlpha,lineWidth,gridStep_deg)

% milkyway.m returns:
%   P             planar reference polygon in equatorial RA/Dec
%   northBoundary ordered northern footprint edge [RA_deg, Dec_deg]
%   southBoundary ordered southern footprint edge [RA_deg, Dec_deg]
%
% The boundary arrays intentionally contain local RA backtracking.  NEVER
% sort them by RA and NEVER reduce them to Dec=f(RA).
[P,northBoundary,southBoundary] = milkyway(); %#ok<ASGLU>

validate_milkyway_boundaries(northBoundary,southBoundary);

h = struct( ...
    'fill',gobjects(0), ...
    'contour',gobjects(0), ...
    'centerline',gobjects(0));


%% ------------------------------------------------------------------------
%  FILLED FOOTPRINT
%
%  The historical footprint is defined in the ordinary equatorial RA/Dec
%  chart.  We therefore classify points there, not after projection.
% -------------------------------------------------------------------------

if doFill

    relativeLon = -180:gridStep_deg:180;
    latitude = -90:gridStep_deg:90;

    if relativeLon(end) < 180
        relativeLon(end+1) = 180;
    end

    if latitude(end) < 90
        latitude(end+1) = 90;
    end

    [relativeLonGrid,latGrid] = meshgrid(relativeLon,latitude);

    lonGrid = mod(centerLongitude_deg + relativeLonGrid,360);

    [RAeq,Deceq] = ...
        selected_sky_to_equatorial_lonlat( ...
            lonGrid,latGrid,coordinateFrame);

    % Build the exact original planar footprint.  This is precisely the
    % geometry represented by milkyway.m: north edge from RA 0->360 and
    % the reversed south edge from RA 360->0.  The only closing segments
    % are at the RA=0/360 chart seam.
    polygonRA = [ ...
        northBoundary(:,1); ...
        flipud(southBoundary(:,1)) ...
        ];

    polygonDec = [ ...
        northBoundary(:,2); ...
        flipud(southBoundary(:,2)) ...
        ];

    inside = reshape( ...
        inpolygon( ...
            mod(RAeq(:),360), ...
            Deceq(:), ...
            polygonRA, ...
            polygonDec), ...
        size(RAeq));

    [xFill,yFill] = ...
        mollweide_relative_lonlat( ...
            relativeLonGrid,latGrid,longitudeDirection);

    alphaData = fillAlpha .* double(inside);

    h.fill = surface( ...
        ax, ...
        xFill, ...
        yFill, ...
        zeros(size(xFill)), ...
        zeros(size(xFill)), ...
        'FaceColor',[0.70 0.70 0.70], ...
        'FaceAlpha','flat', ...
        'AlphaData',alphaData, ...
        'AlphaDataMapping','none', ...
        'EdgeColor','none', ...
        'HandleVisibility','off');

    view(ax,2);

    % Keep the fill behind grid/reference curves and the asteroid cloud.
    try
        uistack(h.fill,'bottom');
    catch
        % uistack behavior can differ between MATLAB graphics releases.
    end

end


%% ------------------------------------------------------------------------
%  TRUE IRREGULAR BOUNDARIES
%
%  Densification is performed in the SOURCE equatorial chart before any
%  3-D rotation.  This preserves the intended piecewise-linear footprint
%  instead of letting sparse transformed vertices be joined by unrelated
%  straight lines in the final projection.
% -------------------------------------------------------------------------

boundaryStep_deg = min(0.25,max(0.05,0.5*gridStep_deg));

northDense = densify_ordered_lonlat_curve(northBoundary,boundaryStep_deg);
southDense = densify_ordered_lonlat_curve(southBoundary,boundaryStep_deg);

sourceCurves = {northDense,southDense};
contourHandles = gobjects(2,1);

for k = 1:2

    curve = sourceCurves{k};

    [lonPlot,latPlot] = ...
        equatorial_lonlat_to_selected_sky( ...
            curve(:,1),curve(:,2),coordinateFrame);

    [x,y] = project_sky_curve_mollweide( ...
        lonPlot,latPlot,centerLongitude_deg,longitudeDirection);

    contourHandles(k) = plot( ...
        ax,x,y,'-', ...
        'Color',[0.40 0.40 0.40], ...
        'LineWidth',lineWidth, ...
        'HandleVisibility','off');

end

h.contour = contourHandles;


%% ------------------------------------------------------------------------
%  INDEPENDENT GALACTIC EQUATOR VALIDATION CURVE
% -------------------------------------------------------------------------

if showGalacticEquator

    h.centerline = draw_galactic_equator( ...
        ax,centerLongitude_deg,coordinateFrame,longitudeDirection,lineWidth);

end

end


%% =========================================================================
%
%  VALIDATE LEGACY MILKY WAY BOUNDARIES
%
% =========================================================================

function validate_milkyway_boundaries(northBoundary,southBoundary)

for q = {northBoundary,southBoundary}

    curve = q{1};

    if size(curve,2) ~= 2 || size(curve,1) < 3
        error('milkyway.m must return N-by-2 north/south boundary arrays.');
    end

    if any(~isfinite(curve),'all')
        error('milkyway.m boundary arrays contain non-finite values.');
    end

    if any(curve(:,2) < -90 | curve(:,2) > 90)
        error('milkyway.m contains declinations outside [-90,90] deg.');
    end

end

% Both supplied historical curves explicitly run from RA=0 to RA=360.
% Do not require strict monotonicity; local RA backtracking is intentional.
if abs(northBoundary(1,1)) > 1e-8 || ...
   abs(southBoundary(1,1)) > 1e-8 || ...
   abs(northBoundary(end,1)-360) > 1e-6 || ...
   abs(southBoundary(end,1)-360) > 1e-6

    warning( ...
        ['Milky Way boundary endpoints are not at the expected RA=0/360 ' ...
         'chart seam. The overlay will still be attempted.']);

end

end


%% =========================================================================
%
%  DENSIFY AN ORDERED LONGITUDE/LATITUDE POLYLINE WITHOUT REORDERING IT
%
% =========================================================================

function dense = densify_ordered_lonlat_curve(curve,maxStep_deg)

N = size(curve,1);

pieces = cell(N-1,1);

for k = 1:N-1

    p0 = curve(k,:);
    p1 = curve(k+1,:);

    maxComponentChange = max(abs(p1-p0));
    nSub = max(1,ceil(maxComponentChange/maxStep_deg));

    t = (0:nSub-1).' ./ nSub;

    pieces{k} = p0 + t.*(p1-p0);

end

dense = vertcat(pieces{:},curve(end,:));

end


%% =========================================================================
%
%  DRAW THE TRUE GALACTIC EQUATOR b=0 IN ANY SELECTED FRAME
%
% =========================================================================

function h = draw_galactic_equator( ...
    ax,centerLongitude_deg,coordinateFrame,longitudeDirection,lineWidth)

l_deg = linspace(0,360,2001).';

[lonPlot,latPlot] = ...
    galactic_lonlat_to_selected_sky( ...
        l_deg,zeros(size(l_deg)),coordinateFrame);

[x,y] = project_sky_curve_mollweide( ...
    lonPlot,latPlot,centerLongitude_deg,longitudeDirection);

h = plot( ...
    ax,x,y,'--', ...
    'Color',[0.22 0.22 0.22], ...
    'LineWidth',max(0.8,0.9*lineWidth), ...
    'HandleVisibility','off');

end


%% =========================================================================
%
%  PROJECT A SKY CURVE WITH EXPLICIT PROJECTION-SEAM BREAKING
%
% =========================================================================

function [x,y] = project_sky_curve_mollweide( ...
    lon_deg,lat_deg,centerLongitude_deg,longitudeDirection)

lon_deg = lon_deg(:);
lat_deg = lat_deg(:);

relativeLon = mod(lon_deg-centerLongitude_deg+180,360)-180;

[x,y] = mollweide_relative_lonlat( ...
    relativeLon,lat_deg,longitudeDirection);

% A physical curve crossing the map seam jumps from approximately +180 to
% -180 deg (or vice versa).  Detect the longitude discontinuity itself,
% rather than using an arbitrary jump threshold in projected x.
seamJump = [false; abs(diff(relativeLon)) > 180];

x(seamJump) = NaN;
y(seamJump) = NaN;

end


%% =========================================================================
%
%  BREAK SEAMS IN A PROJECTED SURFACE
%
% =========================================================================

function [x,y] = break_surface_map_seam(x,y)

for iRow = 1:size(x,1)

    jump = find(abs(diff(x(iRow,:))) > 1.0);

    if isempty(jump)
        continue
    end

    cut = unique([jump jump+1]);

    x(iRow,cut) = NaN;
    y(iRow,cut) = NaN;

end

end


%% =========================================================================
%
%  GALACTIC (l,b) -> SELECTED PLOTTING FRAME
%
% =========================================================================

function [lon_deg,lat_deg] = ...
    galactic_lonlat_to_selected_sky(l_deg,b_deg,coordinateFrame)

inputSize = size(l_deg);

if ~isequal(size(l_deg),size(b_deg))
    error('Galactic longitude and latitude inputs must have matching sizes.');
end

rGal = lonlat_to_cartesian(l_deg(:),b_deg(:));
rEq = galactic_to_equatorial_J2000(rGal);

switch coordinateFrame

    case "equatorial"
        rPlot = rEq;

    case "ecliptic"
        rPlot = equatorial_to_ecliptic_J2000(rEq);

    case "galactic"
        rPlot = rGal;

end

[lonVector,latVector] = cartesian_to_lonlat(rPlot);

lon_deg = reshape(lonVector,inputSize);
lat_deg = reshape(latVector,inputSize);

end


%% =========================================================================
%
%  STANDARD ICRS/J2000 EQUATORIAL <-> GALACTIC ROTATION
%
% =========================================================================

function rEq = galactic_to_equatorial_J2000(rGal)

R_ICRS_TO_GAL = [ ...
   -0.0548755604162154, -0.8734370902348850, -0.4838350155487132
    0.4941094278755837, -0.4448296299600112,  0.7469822444972189
   -0.8676661490190047, -0.1980763734312015,  0.4559837761750669
    ];

% Row-vector convention:
%   rGal = rEq * R^T  ->  rEq = rGal * R
rEq = rGal * R_ICRS_TO_GAL;

end


function rGal = equatorial_to_galactic_J2000(rEq)

R_ICRS_TO_GAL = [ ...
   -0.0548755604162154, -0.8734370902348850, -0.4838350155487132
    0.4941094278755837, -0.4448296299600112,  0.7469822444972189
   -0.8676661490190047, -0.1980763734312015,  0.4559837761750669
    ];

rGal = rEq * R_ICRS_TO_GAL.';

end


%% =========================================================================
%
%  EQUATORIAL RA/DEC -> SELECTED SKY LONGITUDE/LATITUDE
%
% =========================================================================

function [lon_deg,lat_deg] = ...
    equatorial_lonlat_to_selected_sky(RA_deg,Dec_deg,coordinateFrame)

inputSize = size(RA_deg);

if ~isequal(size(RA_deg),size(Dec_deg))
    error('Longitude and latitude inputs must have matching sizes.');
end

rEq = lonlat_to_cartesian(RA_deg(:),Dec_deg(:));

switch coordinateFrame

    case "equatorial"
        rPlot = rEq;

    case "ecliptic"
        rPlot = equatorial_to_ecliptic_J2000(rEq);

    case "galactic"
        rPlot = equatorial_to_galactic_J2000(rEq);

end

[lonVector,latVector] = cartesian_to_lonlat(rPlot);

lon_deg = reshape(lonVector,inputSize);
lat_deg = reshape(latVector,inputSize);

end


%% =========================================================================
%
%  SELECTED SKY LONGITUDE/LATITUDE -> EQUATORIAL RA/DEC
%
% =========================================================================

function [RA_deg,Dec_deg] = ...
    selected_sky_to_equatorial_lonlat(lon_deg,lat_deg,coordinateFrame)

inputSize = size(lon_deg);

if ~isequal(size(lon_deg),size(lat_deg))
    error('Longitude and latitude inputs must have matching sizes.');
end

rSelected = lonlat_to_cartesian(lon_deg(:),lat_deg(:));

switch coordinateFrame

    case "equatorial"
        rEq = rSelected;

    case "ecliptic"
        rEq = ecliptic_to_equatorial_J2000(rSelected);

    case "galactic"
        rEq = galactic_to_equatorial_J2000(rSelected);

end

[RAvec,Decvec] = cartesian_to_lonlat(rEq);

RA_deg = reshape(RAvec,inputSize);
Dec_deg = reshape(Decvec,inputSize);

end


%% =========================================================================
%
%  DRAW J2000 ECLIPTIC
%
% =========================================================================

function draw_ecliptic(ax,centerLongitude_deg,coordinateFrame,longitudeDirection)

lambda = ...
    linspace(0,2*pi,1400).';


rEcl = [ ...
    cos(lambda), ...
    sin(lambda), ...
    zeros(size(lambda)) ...
    ];


rPlot = ...
    ecliptic_vectors_to_selected_frame( ...
        rEcl, ...
        coordinateFrame);


[lon,lat] = ...
    cartesian_to_lonlat(rPlot);


[x,y] = ...
    mollweide_lonlat( ...
        lon, ...
        lat, ...
        centerLongitude_deg, ...
        longitudeDirection);


[x,y] = ...
    break_map_seam(x,y);


plot( ...
    ax, ...
    x, ...
    y, ...
    '--', ...
    'Color',[0.30 0.30 0.30], ...
    'LineWidth',1.05, ...
    'HandleVisibility','off');

end


%% =========================================================================
%
%  PLOT LINEAR FORMATION-NORMAL SKY TRACK(S)
%
%  The overlay vectors are supplied in the J2000 ecliptic inertial frame.
%  They are transformed into the same coordinate frame and Mollweide
%  convention used by the asteroid population.
%
% =========================================================================

function hAll = plot_normal_overlays_mollweide( ...
    ax,overlays,centerLongitude_deg,coordinateFrame,longitudeDirection)

if isempty(overlays)
    hAll = gobjects(0);
    return
end

if ~isstruct(overlays)
    error('NormalOverlay must be a struct or struct array.');
end

hAll = gobjects(0);

for iOverlay = 1:numel(overlays)

    ov = overlays(iOverlay);

    if ~isfield(ov,'vectorsEclipticJ2000') || isempty(ov.vectorsEclipticJ2000)
        warning(['Skipping NormalOverlay element %d because it does not contain ' ...
                 'a nonempty vectorsEclipticJ2000 field.'],iOverlay);
        continue
    end

    nEcl = double(ov.vectorsEclipticJ2000);

    if size(nEcl,2) ~= 3
        error('NormalOverlay.vectorsEclipticJ2000 must be N x 3.');
    end

    good = all(isfinite(nEcl),2);
    nEcl = nEcl(good,:);

    if isempty(nEcl)
        continue
    end

    nNorm = vecnorm(nEcl,2,2);
    nEcl = nEcl ./ nNorm;

    nPlot = ecliptic_vectors_to_selected_frame(nEcl,coordinateFrame);
    [lon,lat] = cartesian_to_lonlat(nPlot);

    [x,y] = project_sky_curve_mollweide( ...
        lon,lat,centerLongitude_deg,longitudeDirection);

    % Style defaults.
    lineColor = [0.85 0.10 0.10];
    lineWidth = 2.2;
    lineStyle = '-';
    showAntipode = false;
    antipodeLineStyle = '--';
    showEndpoints = true;
    showLabel = true;
    labelText = "Linear Floquet normal";

    if isfield(ov,'color') && numel(ov.color) == 3
        lineColor = double(ov.color(:).');
    end
    if isfield(ov,'lineWidth') && isfinite(ov.lineWidth)
        lineWidth = double(ov.lineWidth);
    end
    if isfield(ov,'lineStyle') && strlength(string(ov.lineStyle)) > 0
        lineStyle = char(string(ov.lineStyle));
    end
    if isfield(ov,'showAntipode')
        showAntipode = logical(ov.showAntipode);
    end
    if isfield(ov,'antipodeLineStyle') && strlength(string(ov.antipodeLineStyle)) > 0
        antipodeLineStyle = char(string(ov.antipodeLineStyle));
    end
    if isfield(ov,'showEndpoints')
        showEndpoints = logical(ov.showEndpoints);
    end
    if isfield(ov,'showLabel')
        showLabel = logical(ov.showLabel);
    end
    if isfield(ov,'label') && strlength(string(ov.label)) > 0
        labelText = string(ov.label);
    end

    h = plot( ...
        ax,x,y, ...
        'LineStyle',lineStyle, ...
        'Color',lineColor, ...
        'LineWidth',lineWidth, ...
        'Tag',sprintf('JPL_NormalOverlay_%d_Main',iOverlay), ...
        'HandleVisibility','off');

    hAll(end+1,1) = h; %#ok<AGROW>

    if showAntipode

        nAnti = -nEcl;
        nAntiPlot = ecliptic_vectors_to_selected_frame(nAnti,coordinateFrame);
        [lonA,latA] = cartesian_to_lonlat(nAntiPlot);
        [xA,yA] = project_sky_curve_mollweide( ...
            lonA,latA,centerLongitude_deg,longitudeDirection);

        hA = plot( ...
            ax,xA,yA, ...
            'LineStyle',antipodeLineStyle, ...
            'Color',lineColor, ...
            'LineWidth',0.85*lineWidth, ...
            'Tag',sprintf('JPL_NormalOverlay_%d_Antipode',iOverlay), ...
            'HandleVisibility','off');

        hAll(end+1,1) = hA; %#ok<AGROW>

    end

    if showEndpoints

        finiteIdx = find(isfinite(x) & isfinite(y));

        if ~isempty(finiteIdx)
            i0 = finiteIdx(1);
            i1 = finiteIdx(end);

            plot(ax,x(i0),y(i0),'o', ...
                'MarkerSize',6, ...
                'MarkerFaceColor',lineColor, ...
                'MarkerEdgeColor','k', ...
                'LineWidth',0.7, ...
                'Tag',sprintf('JPL_NormalOverlay_%d_Start',iOverlay), ...
                'HandleVisibility','off');

            plot(ax,x(i1),y(i1),'s', ...
                'MarkerSize',6, ...
                'MarkerFaceColor',lineColor, ...
                'MarkerEdgeColor','k', ...
                'LineWidth',0.7, ...
                'Tag',sprintf('JPL_NormalOverlay_%d_End',iOverlay), ...
                'HandleVisibility','off');
        end

    end

    if showLabel

        finiteIdx = find(isfinite(x) & isfinite(y));

        if ~isempty(finiteIdx)
            imid = finiteIdx(max(1,round(numel(finiteIdx)/2)));

            text( ...
                ax, ...
                x(imid)+0.035, ...
                y(imid)+0.025, ...
                latex_text_local(labelText), ...
                'FontSize',9, ...
                'FontWeight','bold', ...
                'Color',lineColor, ...
                'Interpreter','latex', ...
                'Tag',sprintf('JPL_NormalOverlay_%d_Label',iOverlay), ...
                'Clipping','on');
        end

    end

end

end


%% =========================================================================
%
%  PLOT MAJOR PLANETS IN HELIOCENTRIC SKY VIEW
%
% =========================================================================

function h = plot_major_planets_mollweide( ...
    ax,planetTable,centerLongitude_deg,coordinateFrame,longitudeDirection,showLabels)

R_ecl = [ ...
    planetTable.x_AU, ...
    planetTable.y_AU, ...
    planetTable.z_AU ...
    ];


R_plot = ...
    ecliptic_vectors_to_selected_frame( ...
        R_ecl, ...
        coordinateFrame);


[lon,lat] = ...
    cartesian_to_lonlat(R_plot);


[x,y] = ...
    mollweide_lonlat( ...
        lon, ...
        lat, ...
        centerLongitude_deg, ...
        longitudeDirection);


% Readable conventional display colors.  These colors are independent of
% the asteroid spin-rate colorbar and are used only to distinguish planets.
planetColors = [ ...
    0.42 0.42 0.42   % Mercury
    0.85 0.65 0.28   % Venus
    0.12 0.45 0.90   % Earth
    0.75 0.28 0.18   % Mars
    0.70 0.52 0.34   % Jupiter
    0.82 0.70 0.43   % Saturn
    0.42 0.75 0.82   % Uranus
    0.20 0.35 0.80   % Neptune
    ];


h = ...
    gobjects(height(planetTable),1);


for k = 1:height(planetTable)

    h(k) = plot( ...
        ax, ...
        x(k), ...
        y(k), ...
        'o', ...
        'MarkerSize',8.5, ...
        'MarkerFaceColor',planetColors(k,:), ...
        'MarkerEdgeColor',[0.10 0.10 0.10], ...
        'LineWidth',1.0, ...
        'HandleVisibility','off');


    if showLabels

        text( ...
            ax, ...
            x(k) + 0.035, ...
            y(k), ...
            '$\mathrm{' + latex_text_local(planetTable.Name(k)) + '}$', ...
            'FontSize',8.5, ...
            'FontWeight','bold', ...
            'Color',[0.12 0.12 0.12], ...
            'Interpreter','latex', ...
            'Clipping','on');

    end

end

end


%% =========================================================================
%
%  BREAK MOLLWEIDE MAP SEAM
%
% =========================================================================

function [x,y] = break_map_seam(x,y)

originalSize = ...
    size(x);


xCol = ...
    x(:);

yCol = ...
    y(:);


jump = ...
    [false; abs(diff(xCol)) > 1.0];


xCol(jump) = ...
    NaN;

yCol(jump) = ...
    NaN;


x = ...
    reshape(xCol,originalSize);

y = ...
    reshape(yCol,originalSize);

end


%% =========================================================================
%
%  CHOOSE DIAMETER LEGEND VALUES
%
% =========================================================================

function out = latex_text_local(in)
out = string(in);
out = replace(out,"_","\_");
out = replace(out,"%","\%");
out = replace(out,"&","\&");
out = replace(out,"#","\#");
end


function out = latex_roman_text_local(in)
out = latex_text_local(in);
out = replace(out," ","\;");
out = "\mathrm{" + out + "}";
end



function out = build_overlay_family_summary_latex(normalOverlay)
out = "";

if isempty(normalOverlay) || ~isstruct(normalOverlay)
    return
end

nShown = 0;
familyNames = strings(0,1);
for k = 1:numel(normalOverlay)
    ov = normalOverlay(k);
    if ~isfield(ov,'vectorsEclipticJ2000') || isempty(ov.vectorsEclipticJ2000)
        continue
    end
    nShown = nShown + 1;
    if isfield(ov,'familyName') && strlength(string(ov.familyName)) > 0
        familyNames(end+1,1) = string(ov.familyName); %#ok<AGROW>
    elseif isfield(ov,'familyFile') && strlength(string(ov.familyFile)) > 0
        [~,nm,~] = fileparts(char(string(ov.familyFile)));
        familyNames(end+1,1) = string(nm); %#ok<AGROW>
    end
end

if nShown == 0
    return
end

familyNames = unique(familyNames,'stable');
out = "\quad|\quad N_{\mathrm{orbits}}=" + string(nShown);
if ~isempty(familyNames)
    if numel(familyNames) == 1
        famText = latex_roman_text_local(familyNames(1));
    else
        famText = latex_roman_text_local(join(familyNames,', '));
    end
    out = out + "\quad(\mathrm{family=}\;" + famText + ")";
end
end


function Dexample = choose_diameter_legend(D)

Dmax = ...
    max(D);


candidate = ...
    [1 2 5 10 20 50 100 200 500 1000];


candidate = ...
    candidate(candidate <= Dmax);


if isempty(candidate)

    Dexample = ...
        Dmax;

    return

end


if numel(candidate) <= 4

    Dexample = ...
        candidate;

else

    idx = ...
        round( ...
            linspace( ...
                1, ...
                numel(candidate), ...
                4));


    Dexample = ...
        candidate(idx);

end


Dexample = ...
    unique(Dexample);

end


%% =========================================================================
%
%  VERBOSE PRINT HELPER
%
% =========================================================================

function vfprintf(verbose,varargin)

if verbose
    fprintf(varargin{:});
end

end


%% =========================================================================
%
%  SIMPLE PERCENTILE WITHOUT TOOLBOX DEPENDENCY
%
% =========================================================================

function q = percentile_simple(x,p)

x = ...
    x(isfinite(x));


if isempty(x)

    q = ...
        NaN;

    return

end


x = ...
    sort(x(:));


if p <= 0

    q = ...
        x(1);

    return

elseif p >= 100

    q = ...
        x(end);

    return

end


position = ...
    1 + ...
    (numel(x)-1)*(p/100);


iLo = ...
    floor(position);

iHi = ...
    ceil(position);


if iLo == iHi

    q = ...
        x(iLo);

else

    f = ...
        position - iLo;


    q = ...
        (1-f)*x(iLo) + ...
        f*x(iHi);

end

end


%% =========================================================================
%
%  APPROXIMATE UTC -> TDB JULIAN DATE
%
%  For present-era observations (2017 through the current 2026 epoch),
%
%      TAI - UTC = 37 s
%      TT  - TAI = 32.184 s
%
%  The periodic TDB-TT term is then added at millisecond level.
%
%  This is far more precise than needed for this population map.
%
% =========================================================================

function jdTDB = utc_to_jdtdb_approx(epochUTC)

jdUTC = ...
    posixtime(epochUTC)/86400 + 2440587.5;


referenceDate = ...
    datetime( ...
        2017,1,1, ...
        0,0,0, ...
        'TimeZone','UTC');


if epochUTC < referenceDate

    warning( ...
        ['UTC-to-TDB helper assumes TAI-UTC = 37 s, valid from ' ...
         '2017-01-01 through the current 2026 epoch.']);

end


TAIminusUTC_s = ...
    37.0;


TTminusUTC_s = ...
    TAIminusUTC_s + 32.184;


jdTT = ...
    jdUTC + ...
    TTminusUTC_s/86400;


% Approximate TDB-TT periodic correction.
g_deg = ...
    357.53 + ...
    0.9856003*(jdTT - 2451545.0);


g = ...
    deg2rad(g_deg);


TDBminusTT_s = ...
    0.001657*sin(g) + ...
    0.000022*sin(2*g);


jdTDB = ...
    jdTT + ...
    TDBminusTT_s/86400;

end