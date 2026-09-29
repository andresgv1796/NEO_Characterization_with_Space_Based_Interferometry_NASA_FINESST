function qI = rotate_synodic_to_inertial_rows(qS, t)
%ROTATE_SYNODIC_TO_INERTIAL_ROWS Rotate row vectors from synodic to inertial.
%
% This wrapper uses the user's synodicToInertial utility when possible. The
% velocity input is set to zero and only the transformed position output is
% used, which is correct for directions and baselines. A small explicit
% fallback is retained for single/two-row inputs and for compatibility.
%
% Convention used by synodicToInertial:
%   q_I(t) = R3(+t) q_S(t), with nondimensional CR3BP rotation rate n = 1.

    qS = ensure_row_vectors(qS);
    t = t(:);

    if numel(t) == 1 && size(qS, 1) > 1
        t = repmat(t, size(qS, 1), 1);
    end
    if size(qS, 1) ~= numel(t)
        error('Length of t must match number of row vectors in qS.');
    end

    % The uploaded synodicToInertial utility has a length check based on
    % max(size(X)), so use it for three or more row vectors. For one or two
    % row vectors, use the algebraically identical direct rotation below.
    if exist('synodicToInertial', 'file') == 2 && size(qS, 1) >= 3
        vZero = zeros(size(qS));
        try
            [qI, ~] = synodicToInertial(t, qS, vZero, 1.0);
            return;
        catch ME
            warning('rotate_synodic_to_inertial_rows:UserTransformFailed', ...
                ['synodicToInertial failed, so using direct position ', ...
                 'rotation. Original message: %s'], ME.message);
        end
    end

    c = cos(t);
    s = sin(t);

    qI = zeros(size(qS));
    qI(:, 1) = c .* qS(:, 1) - s .* qS(:, 2);
    qI(:, 2) = s .* qS(:, 1) + c .* qS(:, 2);
    qI(:, 3) = qS(:, 3);
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
