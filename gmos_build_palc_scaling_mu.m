function scale = gmos_build_palc_scaling_mu(y0, y1, N, opts)
%GMOS_BUILD_PALC_SCALING_MU Build a diagonal scaling for weighted PALC in [z; mu].
%
% The weighted pseudo arclength formulation is applied to the scaled variable
%
%   ys = D * y,
%
% where y = [X0(:); T; rho; mu]. This prevents the continuation metric from
% being dominated by the 6N state components and makes the mu component
% comparable to the other blocks.
%
% Default behavior
%   The block scales are inferred from the secant between the two supplied
%   converged seed solutions y0 and y1:
%
%     xScale   = rms block scale of X0(:)
%     TScale   = abs seed difference in T
%     rhoScale = abs seed difference in rho
%     muScale  = abs seed difference in mu
%
%   The returned weight vector is w = diag(D), that is
%
%     ys = w .* y.
%
% Inputs
%   y0, y1 : converged extended solutions [X0(:); T; rho; mu]
%   N      : number of curve nodes
%   opts   : optional struct with any of the fields
%              xScale
%              TScale
%              rhoScale
%              muScale
%              floorX   (default 1e-12)
%              floorT   (default 1e-12)
%              floorRho (default 1e-12)
%              floorMu  (default 1e-16)
%
% Output
%   scale  : struct with fields
%              w        : diagonal weights as a column vector
%              xScale   : physical X block scale
%              TScale   : physical T scale
%              rhoScale : physical rho scale
%              muScale  : physical mu scale
%              Ddesc    : short text description

    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    nY = 6 * N + 3;
    if numel(y0) ~= nY || numel(y1) ~= nY
        error('y0 and y1 must both have length 6*N + 3.');
    end

    if ~isfield(opts, 'floorX');   opts.floorX = 1e-12; end
    if ~isfield(opts, 'floorT');   opts.floorT = 1e-12; end
    if ~isfield(opts, 'floorRho'); opts.floorRho = 1e-12; end
    if ~isfield(opts, 'floorMu');  opts.floorMu = 1e-16; end

    nX = 6 * N;
    dy = y1(:) - y0(:);

    if isfield(opts, 'xScale') && ~isempty(opts.xScale)
        xScale = opts.xScale;
    else
        dx = dy(1:nX);
        xScale = norm(dx) / sqrt(nX);
        if ~isfinite(xScale) || xScale < opts.floorX
            xScale = max(norm(y0(1:nX)) / sqrt(nX), 1.0);
        end
    end

    if isfield(opts, 'TScale') && ~isempty(opts.TScale)
        TScale = opts.TScale;
    else
        TScale = abs(dy(nX + 1));
        if ~isfinite(TScale) || TScale < opts.floorT
            TScale = max(abs(y0(nX + 1)), 1.0);
        end
    end

    if isfield(opts, 'rhoScale') && ~isempty(opts.rhoScale)
        rhoScale = opts.rhoScale;
    else
        rhoScale = abs(dy(nX + 2));
        if ~isfinite(rhoScale) || rhoScale < opts.floorRho
            rhoScale = max(abs(y0(nX + 2)), 1.0);
        end
    end

    if isfield(opts, 'muScale') && ~isempty(opts.muScale)
        muScale = opts.muScale;
    else
        muScale = abs(dy(nX + 3));
        if ~isfinite(muScale) || muScale < opts.floorMu
            error(['muScale inferred from the seed secant is too small. ', ...
                   'Supply two distinct seed solutions with different mu, or ', ...
                   'set opts.muScale explicitly.']);
        end
    end

    w = [repmat(1.0 / xScale, nX, 1); ...
         1.0 / TScale; ...
         1.0 / rhoScale; ...
         1.0 / muScale];

    scale = struct();
    scale.w = w;
    scale.xScale = xScale;
    scale.TScale = TScale;
    scale.rhoScale = rhoScale;
    scale.muScale = muScale;
    scale.Ddesc = 'Weighted PALC in ys = D*y with block scales inferred from the seed secant.';
end
