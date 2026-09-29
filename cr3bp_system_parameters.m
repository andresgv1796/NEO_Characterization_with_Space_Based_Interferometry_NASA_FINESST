function sys = cr3bp_system_parameters(systemName)
%CR3BP_SYSTEM_PARAMETERS Canonical CR3BP scales for supported systems.
%
% sys = cr3bp_system_parameters(systemName)
%
% Supported names:
%   'EarthMoon', 'EM'
%   'SunEarth',  'SE'
%
% Units:
%   Lstar_km     : canonical distance unit [km]
%   Tstar_s      : canonical time unit [s]
%   Vstar_km_s   : canonical velocity unit [km/s]
%   Astar_km_s2  : canonical acceleration unit [km/s^2]
%   primaryPeriod_days = 2*pi*Tstar_s/day
%
% Notes:
%   The nondimensional CR3BP dynamics must use the same mass parameter mu
%   as the periodic-orbit catalog used to select the parent orbit. For that
%   reason, sys.mu is set as the catalog/dynamical mass parameter, not
%   recomputed from the dimensional GM values.
%
%   Dimensional scales are computed from Lstar_km and GMtot_km3_s2:
%
%       Tstar_s = sqrt(Lstar_km^3/GMtot_km3_s2).
%
%   GM values are used for dimensional unit conversion only. Effective
%   primary/secondary GM values consistent with the catalog mu are also
%   stored as GM1_effective_km3_s2 and GM2_effective_km3_s2.
%
%   The nondimensional CR3BP has angular rate n = 1, so one dimensional
%   revolution of the primaries is 2*pi nondimensional time units.

    key = lower(strrep(strtrim(systemName), '-', ''));
    key = strrep(key, ' ', '');

    switch key
        case {'earthmoon', 'em', 'earth_moon'}
            sys.name = 'Earth-Moon';
            sys.shortName = 'EM';
            sys.primaryName = 'Earth';
            sys.secondaryName = 'Moon';

            % JPL three-body periodic-orbit catalog / CR3BP dynamical mu.
            % Keep this exact for compatibility with catalog initial states.
            sys.mu = 1.215058560962404e-2;
            sys.muSource = 'JPL three-body periodic-orbit catalog compatible value';

            % Dimensional scale choices, km and seconds.
            % GM values are directly measured gravitational parameters and
            % are used here to compute dimensional unit conversions.
            sys.GM1_km3_s2 = 398600.4418;       % Earth GM
            sys.GM2_km3_s2 = 4902.800066;       % Moon GM
            sys.GMtot_km3_s2 = sys.GM1_km3_s2 + sys.GM2_km3_s2;
            sys.Lstar_km = 384400.0;
            sys.defaultTimeUnit = 'days';

        case {'sunearth', 'se', 'sun_earth'}
            sys.name = 'Sun-Earth';
            sys.shortName = 'SE';
            sys.primaryName = 'Sun';
            sys.secondaryName = 'Earth';

            % JPL three-body periodic-orbit catalog / CR3BP dynamical mu.
            % Keep this exact for compatibility with catalog initial states.
            sys.mu = 3.0034806e-6;
            sys.muSource = 'JPL three-body periodic-orbit catalog compatible value';

            % Dimensional scale choices, km and seconds.
            % GM values are directly measured gravitational parameters and
            % are used here to compute dimensional unit conversions.
            sys.GM1_km3_s2 = 1.32712440018e11;  % Sun GM
            sys.GM2_km3_s2 = 3.986004418e5;     % Earth GM
            sys.GMtot_km3_s2 = sys.GM1_km3_s2 + sys.GM2_km3_s2;
            sys.Lstar_km = 149597870.7;
            sys.defaultTimeUnit = 'years';

        otherwise
            error('Unsupported systemName "%s". Use EarthMoon or SunEarth.', systemName);
    end

    day_s = 86400.0;
    julianYear_s = 365.25 * day_s;

    % Catalog/dynamics mu remains authoritative. Store the GM-derived value
    % only as a diagnostic, because small differences may exist between the
    % catalog CR3BP value and the chosen dimensional GM constants.
    sys.mu_from_GM = sys.GM2_km3_s2/sys.GMtot_km3_s2;
    sys.muDifference_from_GM = sys.mu - sys.mu_from_GM;
    sys.muRelativeDifference_from_GM = sys.muDifference_from_GM/sys.mu;

    % Effective GM split consistent with the catalog/dynamical mu while
    % preserving the selected total GM used for time scaling.
    sys.GM1_effective_km3_s2 = (1.0 - sys.mu)*sys.GMtot_km3_s2;
    sys.GM2_effective_km3_s2 = sys.mu*sys.GMtot_km3_s2;

    % Canonical dimensional scales.
    sys.Tstar_s = sqrt(sys.Lstar_km^3/sys.GMtot_km3_s2);
    sys.Tstar_days = sys.Tstar_s / day_s;
    sys.Tstar_years = sys.Tstar_s / julianYear_s;

    sys.Vstar_km_s = sys.Lstar_km / sys.Tstar_s;
    sys.Astar_km_s2 = sys.Lstar_km / sys.Tstar_s^2;
    sys.omega_dim_rad_s = 1.0 / sys.Tstar_s;

    % One nondimensional CR3BP revolution corresponds to 2*pi canonical
    % time units.
    sys.primaryPeriod_s = 2*pi*sys.Tstar_s;
    sys.primaryPeriod_days = sys.primaryPeriod_s / day_s;
    sys.primaryPeriod_years = sys.primaryPeriod_s / julianYear_s;

    % Convenience aliases for pipeline dimensional conversions.
    sys.day_s = day_s;
    sys.julianYear_s = julianYear_s;
    sys.km_per_LU = sys.Lstar_km;
    sys.s_per_TU = sys.Tstar_s;
    sys.days_per_TU = sys.Tstar_days;
    sys.years_per_TU = sys.Tstar_years;
    sys.km_s_per_VU = sys.Vstar_km_s;
    sys.km_s2_per_AU = sys.Astar_km_s2;

    sys.defaultLengthUnit = 'km';
    sys.defaultVelocityUnit = 'km/s';
    sys.defaultAccelerationUnit = 'km/s^2';
    sys.defaultAreaUnit = 'km^2';

    % Primary locations in the synodic nondimensional CR3BP frame.
    sys.r1_synodic_nd = [-sys.mu, 0.0, 0.0];
    sys.r2_synodic_nd = [1.0 - sys.mu, 0.0, 0.0];
    sys.r1_synodic_km = sys.Lstar_km * sys.r1_synodic_nd;
    sys.r2_synodic_km = sys.Lstar_km * sys.r2_synodic_nd;
end
