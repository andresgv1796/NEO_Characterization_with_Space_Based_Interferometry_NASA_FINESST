function [X_syn, V_syn] = inertialToSynodic(t, X_in, V_in, n)
% INERTIALTOSYNODIC  Rotate inertial states into the synodic CR3BP frame
%
%   [X_syn, V_syn] = inertialToSynodic(t, X_in, V_in, n)
%   Accepts X_in and V_in as either 3×N or N×3 arrays for vectorized use.
%
%   Inputs:
%     t    – scalar or length‑N vector of times (Nx1 or 1×N)
%     X_in – 3×N or N×3 inertial positions [x; y; z]
%     V_in – 3×N or N×3 inertial velocities [vx; vy; vz]
%     n    – scalar rotation rate about z (rad/s)
%
%   Outputs:
%     X_syn – same orientation as X_in input, synodic positions
%     V_syn – same orientation as V_in input, synodic velocities
%
%   Implements:
%     X_syn = R(–θ) * X_in
%     V_syn = R(–θ) * V_in – Ω × X_syn
%   where θ = n·t, Ω = [0; 0; n]^T.

    % Determine input orientation
    [r,c] = size(X_in);
    if r == 3
        N = c;
        flip = false;
    elseif c == 3
        N = r;
        X_in = X_in.';
        V_in = V_in.';
        flip = true;
    else
        error('X_in must be 3×N or N×3.');
    end

    % Reshape time to 1×N
    t = reshape(t, 1, []);
    if numel(t) ~= N, error('Length of t must match states.'); end

    % Precompute rotation
    theta = n .* t;          % 1×N
    c     = cos(theta);      % 1×N
    s     = sin(theta);      % 1×N

    % Unpack inertial
    xi = X_in(1, :); yi = X_in(2, :); zi = X_in(3, :);
    vxi = V_in(1, :); vyi = V_in(2, :); vzi = V_in(3, :);

    % Rotate position: X_syn = R(-θ)*X_in
    xs =  c .* xi +  s .* yi;
    ys = -s .* xi +  c .* yi;
    zs =  zi;
    X_out = [xs; ys; zs];

    % Rotate velocity: R(-θ)*V_in
    vxr =  c .* vxi +  s .* vyi;
    vyr = -s .* vxi +  c .* vyi;
    vzr =  vzi;

    % Subtract Coriolis: -Ω×X_syn = [ +n*ys; -n*xs; 0 ]
    vx_syn = vxr + n .* ys;
    vy_syn = vyr - n .* xs;
    vz_syn = vzr;
    V_out = [vx_syn; vy_syn; vz_syn];

    % Restore original orientation
    if flip
        X_syn = X_out.';
        V_syn = V_out.';
    else
        X_syn = X_out;
        V_syn = V_out;
    end
end
