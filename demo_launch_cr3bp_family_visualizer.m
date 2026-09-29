% DEMO_LAUNCH_CR3BP_FAMILY_VISUALIZER
% Simple FINESST demo: choose one JPL periodic-orbit family, compute its
% Floquet-normal overlays, and launch the live visualizer.
%
% This intentionally stays as a plain script: no local helper functions.

clear;
clc;

epochUTC = datetime('now','TimeZone','UTC');

%% USER SELECTIONS
% Choose the CR3BP system whose JPL CSV files you have downloaded.
systemName = "SunEarth";      % "SunEarth" or "EarthMoon"

% Representative case for each JPL Three-Body Periodic Orbits family type:
%
%   "halo_l1_n"
%   "vertical_l1"
%   "axial_l1"
%   "lyapunov_l1"
%   "longp_l4"
%   "short_l4"
%   "butterfly_n"
%   "dragonfly_n"
%   "resonant_21"
%   "dro"
%   "dpo"
%   "lpo_e"
%
demoCase = "halo_l1_n";

% Keep the demo deliberately simple. Change this directly if you want a
% denser sampling of a large continuation family.
parentRows = 1:10:1e3;

centerMode = "planar";      % change to "planar" when desired
NPrimaryPeriods = 10;
SamplesPerPeriod = 600;

% For the current idealized family-comparison study, heliocentric/ecliptic
% is the common fixed asteroid-reference sky. You can still change this to
% "geocentric" for the complementary observer-view plot.
observer = "heliocentric";
population = "all_asteroids";

%% FAMILY FILE
% These filenames follow the same FINESST export convention already used by
% Butterfly, northern Halo, and DRO examples. If one of your local exports
% has a slightly different filename, just edit that one line in the case.

if systemName == "SunEarth"
    systemSuffix = "Sun_Earth";
    showPrimary = false;      % hide Sun in the local orbit panel
    showSecondary = true;     % Earth, physical 1:1 scale
elseif systemName == "EarthMoon"
    systemSuffix = "Earth_Moon";
    showPrimary = true;       % Earth
    showSecondary = true;     % Moon
else
    error('systemName must be "SunEarth" or "EarthMoon".');
end

switch demoCase

    case "halo_l1_n"
        % JPL family = halo, libration point L1, northern branch
        poFile = "L1_northern_halo_orbits_JPL_IC_" + systemSuffix + ".csv";

    case "vertical_l1"
        % JPL family = vertical, libration point L1
        poFile = "L1_vertical_orbits_JPL_IC_" + systemSuffix + ".csv";

    case "axial_l1"
        % JPL family = axial, libration point L1
        poFile = "L1_axial_orbits_JPL_IC_" + systemSuffix + ".csv";

    case "lyapunov_l1"
        % JPL family = lyapunov, libration point L1
        poFile = "L1_lyapunov_orbits_JPL_IC_" + systemSuffix + ".csv";

    case "longp_l4"
        % JPL family = longp, libration point L4
        poFile = "L4_long_period_orbits_JPL_IC_" + systemSuffix + ".csv";

    case "short_l4"
        % JPL family = short, libration point L4
        poFile = "L4_short_period_orbits_JPL_IC_" + systemSuffix + ".csv";

    case "butterfly_n"
        % JPL family = butterfly, northern branch
        poFile = "Butterfly_orbits_JPL_IC_" + systemSuffix + ".csv";

    case "dragonfly_n"
        % JPL family = dragonfly, northern branch
        poFile = "Dragonfly_orbits_JPL_IC_" + systemSuffix + ".csv";

    case "resonant_21"
        % JPL family = resonant, 2:1 branch
        poFile = "2_1_resonant_orbits_JPL_IC_" + systemSuffix + ".csv";

    case "dro"
        % JPL family = distant retrograde orbit
        poFile = "DRO_orbits_JPL_IC_" + systemSuffix + ".csv";

    case "dpo"
        % JPL family = distant prograde orbit
        poFile = "DPO_orbits_JPL_IC_" + systemSuffix + ".csv";

    case "lpo_e"
        % JPL family = LPO, eastern branch
        poFile = "LPO_east_orbits_JPL_IC_" + systemSuffix + ".csv";

    otherwise
        error('Unknown demoCase: %s',demoCase);
end

if ~isfile(poFile)
    error(['Could not find "' char(poFile) '". ' ...
        'Check the filename of your local JPL export for this family.']);
end

%% COMPUTE FLOQUET-NORMAL OVERLAYS
overlays = make_cr3bp_floquet_normal_overlays( ...
    systemName,poFile,parentRows,epochUTC, ...
    CenterMode=centerMode, ...
    NPrimaryPeriods=NPrimaryPeriods, ...
    SamplesPerPeriod=SamplesPerPeriod, ...
    MassRatioSource="catalog", ...
    ShowAntipode=true, ...
    ColorMode="gradient", ...
    ColorMapName="parula", ...
    ColorParameter="parentRow", ...
    Verbose=true);

% The visualizer can now filter the selected center eigenvectors by their
% actual character: All / Vertical / Planar / Undecided. The default
% ambiguity band is 0.40 < f_v < 0.60.
%
%% LIVE VISUALIZER
handles = launch_cr3bp_family_visualizer( ...
    systemName,poFile,overlays,epochUTC, ...
    Observer=observer, ...
    Population=population, ...
    CoordinateFrame="ecliptic", ...
    ShowMilkyWay=true, ...
    ShowPlanets=false, ...
    ShowPrimaryBody=showPrimary, ...
    ShowSecondaryBody=showSecondary, ...
    Verbose=false);
