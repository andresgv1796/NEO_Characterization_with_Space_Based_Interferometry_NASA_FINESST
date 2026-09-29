function po = choose_periodic_orbit(poFile, orbitIndex, sys)
%CHOOSE_PERIODIC_ORBIT Load one periodic orbit from a JPL-style CSV file.

    if ~isfile(poFile)
        located = which(poFile);
        if isempty(located)
            error('Could not find periodic orbit file "%s" on path.', poFile);
        end
        poFile = located;
    end

    po = pick_po_from_JPL(poFile, orbitIndex);
    po.file = poFile;
    po.index = orbitIndex;
    po.mu = sys.mu;
    po.T_days = po.T * sys.Tstar_days;
    po.T_years = po.T * sys.Tstar_years;
    po.T_seconds = po.T * sys.Tstar_s;
    po.x0_dim_km_kms = [po.x0(1:3).' * sys.Lstar_km, po.x0(4:6).' * sys.Vstar_km_s];
end
