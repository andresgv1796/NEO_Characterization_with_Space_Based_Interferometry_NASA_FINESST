function [X_in, V_in] = synodicToInertial(t, X_syn, V_syn, n)
% SYNODICTOINERTIAL  Rotate synodic states into the inertial CR3BP frame
%
%   [X_in, V_in] = synodicToInertial(t, X_syn, V_syn, n)
%   Accepts X_syn and V_syn as either 3×N or N×3 arrays for vectorized use.
%
%   Inputs:
%     t     – scalar or length‑N vector of times (Nx1 or 1×N)
%     X_syn – 3×N or N×3 synodic positions [x; y; z]
%     V_syn – 3×N or N×3 synodic velocities [vx; vy; vz]
%     n     – scalar rotation rate about z (rad/s)
%
%   Outputs:
%     X_in  – same orientation as X_syn input, inertial positions
%     V_in  – same orientation as V_syn input, inertial velocities  %Non
%     dimensional
%
%   Implements:
%     X_in = R(θ) * X_syn
%     V_in = R(θ) * V_syn + Ω × X_in
%   where θ = n·t, Ω = [0; 0; n]^T.

    % Detect input orientation and convert to 3×N
    [rows, cols] = size(X_syn);
    if rows == 3
        N = cols;
        transposeFlag = false;
    elseif cols == 3
        N = rows;
        X_syn = X_syn.';   % now 3×N
        V_syn = V_syn.';   % now 3×N
        transposeFlag = true;
    else
        error('X_syn must be 3×N or N×3.');
    end


    % Ensure t is row vector of length N
    t = reshape(t, 1, []);
    if numel(t) ~= max(size(X_syn))
        error('Length of t must match length of states.');
    end

    % 1) Compute rotation angles and trig terms
    theta = n .* t;          % 1×N
    c     = cos(theta);      % 1×N
    s     = sin(theta);      % 1×N

    % 2) Unpack synodic arrays for vector operations
    xs  = X_syn(1, :);
    ys  = X_syn(2, :);
    zs  = X_syn(3, :);
    vxs = V_syn(1, :);
    vys = V_syn(2, :);
    vzs = V_syn(3, :);

    % 3) Rotate position: [ c -s; s  c ] * [xs; ys]
    xi =  c .* xs  -  s .* ys;
    yi =  s .* xs  +  c .* ys;
    zi =  zs;              % z unchanged by rotation about z
    X_out = [xi; yi; zi];  % 3×N

    % 4) Rotate velocity: same 2×2 block on (vxs, vys)
    vxr =  c .* vxs  -  s .* vys;
    vyr =  s .* vxs  +  c .* vys;
    vzr =  vzs;

    % 5) Transport term: Ω×X_in = [-n*yi; +n*xi; 0]
    vx_in = vxr - n .* yi;
    vy_in = vyr + n .* xi;
    vz_in = vzr;
    V_out = [vx_in; vy_in; vz_in];

    % Convert back to original orientation
    if transposeFlag
        X_in = X_out.';
        V_in = V_out.';
    else
        X_in = X_out;
        V_in = V_out;
    end
end
