function z = get_family_member_z(fam, idx)
%GET_FAMILY_MEMBER_Z Extract one corrected GMOS solution vector from a family struct.

    if ~isstruct(fam) || ~isfield(fam, 'z')
        error('Expected a family struct with field .z.');
    end

    nMembers = numel(fam.z);
    if isempty(idx)
        idx = nMembers;
    end

    if idx < 1 || idx > nMembers
        error('qpoIndex out of range. Family has %d members.', nMembers);
    end

    if iscell(fam.z)
        z = fam.z{idx};
    else
        if isvector(fam.z)
            z = fam.z(:);
        else
            z = fam.z(:, idx);
        end
    end

    z = z(:);
end
