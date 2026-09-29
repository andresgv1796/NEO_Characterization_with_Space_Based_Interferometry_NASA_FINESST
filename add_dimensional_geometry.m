function geom = add_dimensional_geometry(geom, sys, timePlotUnit)
%ADD_DIMENSIONAL_GEOMETRY Add km, km^2, seconds, days, years to geom.

    if nargin < 3 || isempty(timePlotUnit) || strcmpi(timePlotUnit, 'auto')
        timePlotUnit = sys.defaultTimeUnit;
    end
    timePlotUnit = lower(timePlotUnit);

    L = sys.Lstar_km;
    T = sys.Tstar_s;

    dim = struct();
    dim.system = sys;
    dim.lengthUnit = 'km';
    dim.areaUnit = 'km^2';
    dim.velocityUnit = 'km/s';
    dim.time_s = geom.t * T;
    dim.time_days = dim.time_s / sys.day_s;
    dim.time_years = dim.time_s / sys.julianYear_s;
    dim.TPrimary_s = geom.config.TPrimary * T;
    dim.TPrimary_days = dim.TPrimary_s / sys.day_s;
    dim.TPrimary_years = dim.TPrimary_s / sys.julianYear_s;

    switch timePlotUnit
        case {'s','sec','second','seconds'}
            dim.timePlot = dim.time_s;
            dim.timePlotUnit = 's';
            dim.timePlotLabel = '$t$ [s]';
        case {'day','days','d'}
            dim.timePlot = dim.time_days;
            dim.timePlotUnit = 'days';
            dim.timePlotLabel = '$t$ [days]';
        case {'year','years','yr','yrs'}
            dim.timePlot = dim.time_years;
            dim.timePlotUnit = 'years';
            dim.timePlotLabel = '$t$ [yr]';
        case {'tau','periods','primaryperiods'}
            dim.timePlot = geom.tau;
            dim.timePlotUnit = 'primary periods';
            dim.timePlotLabel = '$t/T_C$';
        otherwise
            error('Unknown timePlotUnit "%s".', timePlotUnit);
    end

    dim.r_km = geom.r * L;
    dim.rC_km = geom.rC * L;
    dim.v_km_s = geom.v * sys.Vstar_km_s;
    dim.vC_km_s = geom.vC * sys.Vstar_km_s;
    dim.YC_km_kms = [geom.YC(:, 1:3) * L, geom.YC(:, 4:6) * sys.Vstar_km_s];
    dim.YCollectors_km_kms = geom.YCollectors;
    dim.YCollectors_km_kms(:, 1:3, :) = geom.YCollectors(:, 1:3, :) * L;
    dim.YCollectors_km_kms(:, 4:6, :) = geom.YCollectors(:, 4:6, :) * sys.Vstar_km_s;

    names = fieldnames(geom.baselines);
    dim.baselines = struct();
    for k = 1:numel(names)
        name = names{k};
        val = geom.baselines.(name);
        if isnumeric(val)
            if startsWith(name, 'b')
                dim.baselines.([name '_km']) = val * L;
            else
                dim.baselines.([name '_km']) = val * L;
            end
        end
    end

    dim.centroid_km = geom.centroid * L;
    dim.centroidMinusCombiner_km = geom.centroidMinusCombiner * L;
    dim.centroidOffset_km = geom.centroidOffset * L;
    dim.collectorMinusCentroid_km = geom.collectorMinusCentroid * L;
    dim.collectorRadius_km = geom.collectorRadius * L;
    dim.rmsCollectorRadius_km = geom.rmsCollectorRadius * L;
    dim.inPlaneOffsetVec_km = geom.inPlaneOffsetVec * L;
    dim.inPlaneOffset_km = geom.inPlaneOffset * L;
    dim.absPlaneOffset_km = geom.absPlaneOffset * L;
    dim.combinerRangeMean_km = geom.combinerRangeMean * L;
    dim.combinerRangeSpread_km = geom.combinerRangeSpread * L;
    dim.area_km2 = geom.area * L^2;
    dim.twoArea_km2 = geom.twoArea * L^2;

    geom.dim = dim;
end
