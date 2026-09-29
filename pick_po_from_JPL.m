function po = pick_po_from_JPL(filename, idx)
% pick_po_from_JPL  Convenience wrapper around load_JPL_orbit_data
%
% po.x0 : 6x1
% po.T  : scalar

    [IC, tF] = load_JPL_orbit_data(filename);

    if idx < 1 || idx > size(IC,1)
        error('idx out of range. File has %d orbits.', size(IC,1));
    end

    po = struct();
    po.x0 = IC(idx,:).';
    po.T  = tF(idx);
end