function JC = Jacobi_Constant(state, varargin)
%JACOBI_CONSTANT Compute the Jacobi constant for one or many CR3BP states.
%   JC = JACOBI_CONSTANT(state, mu) accepts:
%     • state: 6×1, 1×6 vector or N×6 matrix of [x y z vx vy vz]
%     • mu   : CR3BP mass ratio (scalar)
%
%   JC = JACOBI_CONSTANT(x, y, z, vx, vy, vz, mu) accepts each
%   component as a separate input.  x, y, z, vx, vy, vz can be
%   scalars or N×1 vectors (all same length), mu is scalar.
%
%   Output:
%     JC    – N×1 column vector of Jacobi constants
%
%   Examples:
%     % Single state as vector:
%     JC1 = Jacobi_Constant([1, 0.1, 0, 0, 0.01, 0], 0.01215);
%
%     % Many states as N×6 matrix:
%     S = [1,0.1,0,0,0.01,0;  1.1,0.0,0,0,0.02,0];
%     JCv = Jacobi_Constant(S, 0.01215);
%
%     % Separate inputs, vectorized:
%     x  = [1; 1.1];
%     y  = [0.1; 0];
%     z  = [0; 0];
%     vx = [0; 0];
%     vy = [0.01; 0.02];
%     vz = [0; 0];
%     JCv = Jacobi_Constant(x, y, z, vx, vy, vz, 0.01215);

    % --- Parse inputs ---
    if nargin == 2
        % state + mu
        mu    = varargin{1};
        S     = state;
        % allow 6×1 or 1×6
        if numel(S) == 6
            S = reshape(S, 1, 6);
        end
        if size(S,2) ~= 6
            error('State must be a 6-element vector or an N×6 matrix.');
        end
        x  = S(:,1);
        y  = S(:,2);
        z  = S(:,3);
        vx = S(:,4);
        vy = S(:,5);
        vz = S(:,6);

    elseif nargin == 7
        % x, y, z, vx, vy, vz, mu
        x  = state;
        y  = varargin{1};
        z  = varargin{2};
        vx = varargin{3};
        vy = varargin{4};
        vz = varargin{5};
        mu = varargin{6};
        % ensure column vectors
        x  = x(:);
        y  = y(:);
        z  = z(:);
        vx = vx(:);
        vy = vy(:);
        vz = vz(:);
        % check lengths
        N = numel(x);
        if any([numel(y), numel(z), numel(vx), numel(vy), numel(vz)] ~= N)
            error('All state components must have the same number of elements.');
        end

    else
        error(['Invalid call. Use either Jacobi_Constant(state, mu) ', ...
               'or Jacobi_Constant(x,y,z,vx,vy,vz,mu).']);
    end

    % --- Compute Jacobi constant ---
    n = 1;  % non‑dimensional mean motion

    r13 = sqrt((x + mu).^2 + y.^2 + z.^2);
    r23 = sqrt((x - 1 + mu).^2 + y.^2 + z.^2);

    U  = (1 - mu)./r13 + mu./r23 + 0.5 * n^2 .* (x.^2 + y.^2);

    JC = 2 .* U - (vx.^2 + vy.^2 + vz.^2);
end
