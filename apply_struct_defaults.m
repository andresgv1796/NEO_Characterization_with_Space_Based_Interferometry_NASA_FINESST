function s = apply_struct_defaults(s, defaults)
%APPLY_STRUCT_DEFAULTS Fill missing fields in struct s from defaults.

    if isempty(s)
        s = struct();
    end
    names = fieldnames(defaults);
    for k = 1:numel(names)
        name = names{k};
        if ~isfield(s, name) || isempty(s.(name))
            s.(name) = defaults.(name);
        end
    end
end
