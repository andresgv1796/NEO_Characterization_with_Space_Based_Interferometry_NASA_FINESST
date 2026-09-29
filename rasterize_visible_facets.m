function img = rasterize_visible_facets( ...
    proj, vis, nPix, halfWidth, facetRadiance)
%RASTERIZE_VISIBLE_FACETS
% Rasterize observer-facing triangular facets onto the target sky plane.
%
%   img = rasterize_visible_facets( ...
%       proj, vis, nPix, halfWidth, facetRadiance)
%
% The projected geometry is assumed to already be expressed in the
% target-centered interferometric sky frame:
%
%       +x_S = +eU  -> +l
%       +y_S = +eV  -> +m
%       +z_S = +eW  -> observer-to-target phase-center direction
%
% Therefore the observer, as seen from the target, lies along -eW.
% This direction is stored in
%
%       proj.observerDirection
%
% by project_shape_orthographic.
%
% Hidden-surface removal is performed with an observer-directed depth
%
%       d = r dot observerDirection.
%
% The facet having the largest d at a given pixel is closest to the
% observer and therefore determines that pixel's radiance.
%
%
% Inputs
% ------
% proj.vertices3D :
%       Nv x 3 mesh vertices expressed in the sky frame [km].
%
% proj.x, proj.y :
%       Projected vertex coordinates along +eU and +eV [km].
%
% proj.faces :
%       Nf x 3 triangular connectivity.
%
% proj.observerDirection :
%       Direction from target toward observer in the sky frame.
%       Under the current FINESST convention this is [0 0 -1].
%
% vis.isVisible :
%       Nf x 1 logical mask from find_visible_facets.
%
% nPix :
%       Number of pixels along each image dimension.
%
% halfWidth :
%       Image half-width in physical sky-plane coordinates [km].
%       The raster spans [-halfWidth,+halfWidth] in both x and y.
%
% facetRadiance :
%       Either:
%
%           scalar       -> same radiance for every facet
%
%       or
%
%           Nf x 1       -> one radiance value per facet
%
%       For the current morphology-only model, use facetRadiance = 1.
%
%
% Outputs
% -------
% img.I :
%       Normalized image used by the visibility model.
%       sum(img.I(:)) = 1.
%
% img.Iraw :
%       Unnormalized rasterized radiance image.
%
% img.x, img.y :
%       Physical sky-plane image axes [km].
%
% img.X, img.Y :
%       Meshgrid arrays corresponding to img.x and img.y.
%
% img.depthBuffer :
%       Observer-directed depth of the visible surface at each pixel.
%
% img.facetIndex :
%       Index of the visible facet occupying each pixel.
%       Zero indicates no asteroid surface.
%
% img.mask :
%       Logical mask of pixels occupied by the target.
%
% img.pixelArea :
%       Physical area represented by one image pixel [km^2].
%
% img.totalFlux :
%       Integrated unnormalized projected flux-like quantity
%
%           sum(Iraw) * pixelArea.
%
%       For constant radiance this is approximately the projected area
%       times that radiance.
%
% Notes
% -----
% This is a pixel-center rasterizer. Boundary pixels are treated as either
% fully occupied or unoccupied. Therefore the silhouette converges with
% increasing image resolution but is not fractionally antialiased.

arguments
    proj struct
    vis struct
    nPix (1,1) double {mustBeInteger, mustBeGreaterThan(nPix,1)}
    halfWidth (1,1) double {mustBePositive}
    facetRadiance double = 1
end

%% Validate required projected geometry

requiredProjFields = { ...
    'vertices3D', ...
    'faces', ...
    'observerDirection'};

for k = 1:numel(requiredProjFields)
    if ~isfield(proj, requiredProjFields{k})
        error('proj.%s is required.', requiredProjFields{k});
    end
end

if ~isfield(vis, 'isVisible')
    error('vis.isVisible is required.');
end

V = proj.vertices3D;
F = proj.faces;

Nf = size(F,1);

if numel(vis.isVisible) ~= Nf
    error('vis.isVisible must contain one entry per mesh facet.');
end


%% Observer direction

observerDirection = proj.observerDirection(:);

observerNorm = norm(observerDirection);

if ~isfinite(observerNorm) || observerNorm == 0
    error('proj.observerDirection must be a finite nonzero vector.');
end

observerDirection = observerDirection / observerNorm;

% In the current sky-frame formulation the observer direction must be
% perpendicular to the x-y image plane.
%
% Normally:
%
%       observerDirection = [0; 0; -1]
%
observerPlaneTolerance = 1e-12;

if abs(observerDirection(1)) > observerPlaneTolerance || ...
   abs(observerDirection(2)) > observerPlaneTolerance

    error(['proj.observerDirection must be perpendicular to the ', ...
           'projected x-y sky plane.']);
end


%% Radiance definition

if isscalar(facetRadiance)

    facetRadiance = repmat(facetRadiance, Nf, 1);

else

    facetRadiance = facetRadiance(:);

    if numel(facetRadiance) ~= Nf
        error(['facetRadiance must be either scalar or contain ', ...
               'one value per facet.']);
    end

end

if any(~isfinite(facetRadiance))
    error('facetRadiance contains nonfinite values.');
end


%% Image grid

x = linspace(-halfWidth, halfWidth, nPix);
y = linspace(-halfWidth, halfWidth, nPix);

[X,Y] = meshgrid(x,y);

dx = x(2) - x(1);
dy = y(2) - y(1);

pixelArea = dx*dy;


%% Initialize raster

Iraw = zeros(nPix,nPix);

% Larger observer-directed depth means closer to the observer.
depthBuffer = -inf(nPix,nPix);

% Zero means no facet occupies the pixel.
facetIndex = zeros(nPix,nPix);


%% Rasterize visible facets

visibleIdx = find(vis.isVisible);

baryTolerance = 1e-12;

for q = 1:numel(visibleIdx)

    k = visibleIdx(q);

    face = F(k,:);

    % ---------------------------------------------------------------------
    % Triangle vertices in the target-centered sky frame.
    % ---------------------------------------------------------------------
    rv = V(face,:);

    xv = rv(:,1);
    yv = rv(:,2);

    % Observer-directed depth of each triangle vertex.
    %
    % With observerDirection = [0 0 -1],
    %
    %       depth = -z_S.
    %
    % Therefore smaller z_S corresponds to a surface closer to the
    % observer, exactly as required by the chosen +eW convention.
    dv = rv * observerDirection;


    % ---------------------------------------------------------------------
    % Projected triangle bounding box.
    % ---------------------------------------------------------------------
    xmin = min(xv);
    xmax = max(xv);

    ymin = min(yv);
    ymax = max(yv);

    % Triangle completely outside raster.
    if xmax < x(1) || xmin > x(end) || ...
       ymax < y(1) || ymin > y(end)

        continue

    end


    % ---------------------------------------------------------------------
    % Convert physical bounding box to pixel-index bounding box.
    % ---------------------------------------------------------------------
    iMin = max(1, floor((xmin - x(1))/dx) + 1);
    iMax = min(nPix, ceil((xmax - x(1))/dx) + 1);

    jMin = max(1, floor((ymin - y(1))/dy) + 1);
    jMax = min(nPix, ceil((ymax - y(1))/dy) + 1);

    Xsub = X(jMin:jMax, iMin:iMax);
    Ysub = Y(jMin:jMax, iMin:iMax);


    % ---------------------------------------------------------------------
    % Barycentric coordinates of the pixel centers.
    %
    % The projected triangle is described by
    %
    %       P = lambda1*P1 + lambda2*P2 + lambda3*P3
    %
    % with
    %
    %       lambda1 + lambda2 + lambda3 = 1.
    %
    % A point lies inside the triangle when all three lambda values are
    % nonnegative, to numerical tolerance.
    % ---------------------------------------------------------------------
    denom = ...
        (yv(2)-yv(3))*(xv(1)-xv(3)) + ...
        (xv(3)-xv(2))*(yv(1)-yv(3));

    % Edge-on facets can collapse to zero projected area.
    if abs(denom) <= eps(max(1,abs(denom)))
        continue
    end

    lambda1 = ...
        ((yv(2)-yv(3)).*(Xsub-xv(3)) + ...
         (xv(3)-xv(2)).*(Ysub-yv(3))) ./ denom;

    lambda2 = ...
        ((yv(3)-yv(1)).*(Xsub-xv(3)) + ...
         (xv(1)-xv(3)).*(Ysub-yv(3))) ./ denom;

    lambda3 = 1 - lambda1 - lambda2;

    inside = ...
        lambda1 >= -baryTolerance & ...
        lambda2 >= -baryTolerance & ...
        lambda3 >= -baryTolerance;

    if ~any(inside,'all')
        continue
    end


    % ---------------------------------------------------------------------
    % Interpolate observer-directed depth across the projected triangle.
    %
    % Since the original 3-D facet is planar, barycentric interpolation
    % gives the correct depth at each projected pixel center.
    % ---------------------------------------------------------------------
    Dsub = ...
        lambda1*dv(1) + ...
        lambda2*dv(2) + ...
        lambda3*dv(3);


    % ---------------------------------------------------------------------
    % Z/depth-buffer test.
    %
    % The largest observer-directed depth is the frontmost surface.
    % ---------------------------------------------------------------------
    oldDepth = depthBuffer(jMin:jMax, iMin:iMax);

    update = inside & (Dsub > oldDepth);

    if ~any(update,'all')
        continue
    end


    % Update depth buffer.
    oldDepth(update) = Dsub(update);

    depthBuffer(jMin:jMax, iMin:iMax) = oldDepth;


    % ---------------------------------------------------------------------
    % Assign radiance of frontmost facet.
    % ---------------------------------------------------------------------
    Iblock = Iraw(jMin:jMax, iMin:iMax);

    Iblock(update) = facetRadiance(k);

    Iraw(jMin:jMax, iMin:iMax) = Iblock;


    % ---------------------------------------------------------------------
    % Record frontmost facet index.
    % ---------------------------------------------------------------------
    Fblock = facetIndex(jMin:jMax, iMin:iMax);

    Fblock(update) = k;

    facetIndex(jMin:jMax, iMin:iMax) = Fblock;

end


%% Final image mask

mask = facetIndex > 0;


%% Preserve unnormalized projected flux

totalFlux = sum(Iraw(:))*pixelArea;


%% Normalize for interferometric visibility calculations

imageSum = sum(Iraw(:));

if imageSum <= 0
    error(['Rasterization produced zero image intensity. Check the ', ...
           'observer direction, facet visibility, and image field of view.']);
end

I = Iraw / imageSum;


%% Package output

img.I = I;
img.Iraw = Iraw;

img.x = x;
img.y = y;

img.X = X;
img.Y = Y;

img.depthBuffer = depthBuffer;

img.facetIndex = facetIndex;
img.mask = mask;

img.pixelArea = pixelArea;
img.totalFlux = totalFlux;

img.halfWidth = halfWidth;
img.nPix = nPix;

img.observerDirection = observerDirection.';

end