function x = normalize_row(x)
%NORMALIZE_ROW Normalize one 1x3 row vector.

    x = x(:).';
    nx = norm(x);
    if nx <= eps
        error('Cannot normalize a zero vector.');
    end
    x = x / nx;
end
