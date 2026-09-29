function normal = enforce_normal_sign_continuity(normalRaw, normalReference)
%ENFORCE_NORMAL_SIGN_CONTINUITY Resolve n versus -n ambiguity continuously.

    normal = normalRaw;
    nGrid = size(normal, 1);

    firstGood = find(all(isfinite(normal), 2), 1, 'first');
    if isempty(firstGood)
        warning('No valid normals found. Collector triangle may be degenerate.');
        return;
    end

    for k = firstGood+1:nGrid
        if any(~isfinite(normal(k, :)))
            normal(k, :) = normal(k-1, :);
            continue;
        end

        if dot(normal(k, :), normal(k-1, :)) < 0
            normal(k, :) = -normal(k, :);
        end
    end

    if nargin >= 2 && ~isempty(normalReference)
        ref = normalize_row(normalReference);
        if dot(normal(firstGood, :), ref) < 0
            normal = -normal;
        end
    end
end
