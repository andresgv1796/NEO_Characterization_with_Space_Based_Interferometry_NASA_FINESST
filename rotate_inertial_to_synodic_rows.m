function qS = rotate_inertial_to_synodic_rows(qI, t)
%ROTATE_INERTIAL_TO_SYNODIC_ROWS Rotate row vectors from inertial to synodic.
%
% This wrapper uses the user's inertialToSynodic utility when possible. The
% velocity input is set to zero and only the transformed position output is
% used, which is correct for directions and baselines. A small explicit
% fallback is retained for single/two-row inputs and for compatibility.
%
% Convention used by inertialToSynodic:
%   q_S(t) = R3(-t) q_I(t), with nondimensional CR3BP rotation rate n = 1.

    qI = ensure_row_vectors(qI);
    t = t(:);

    if numel(t) == 1 && size(qI, 1) > 1
        t = repmat(t, size(qI, 1), 1);
    end
    if size(qI, 1) ~= numel(t)
        error('Length of t must match number of row vectors in qI.');
    end

    % The uploaded inertialToSynodic utility has a length check based on
    % max(size(X)), so use it for three or more row vectors. For one or two
    % row vectors, use the algebraically identical direct rotation below.
    if exist('inertialToSynodic', 'file') == 2 && size(qI, 1) >= 3
        vZero = zeros(size(qI));
        try
            [qS, ~] = inertialToSynodic(t, qI, vZero, 1.0);
            return;
        catch ME
            warning('rotate_inertial_to_synodic_rows:UserTransformFailed', ...
                ['inertialToSynodic failed, so using direct position ', ...
                 'rotation. Original message: %s'], ME.message);
        end
    end

    c = cos(t);
    s = sin(t);

    qS = zeros(size(qI));
    qS(:, 1) = c .* qI(:, 1) + s .* qI(:, 2);
    qS(:, 2) = -s .* qI(:, 1) + c .* qI(:, 2);
    qS(:, 3) = qI(:, 3);
end

function q = ensure_row_vectors(q)
    if size(q, 2) == 3
        return;
    elseif size(q, 1) == 3
        q = q.';
    else
        error('Input must be N x 3 or 3 x N.');
    end
end
