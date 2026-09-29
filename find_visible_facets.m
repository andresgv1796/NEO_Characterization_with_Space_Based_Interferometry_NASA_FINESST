function vis = find_visible_facets(proj)
%FIND_VISIBLE_FACETS Compute facet geometry and observer visibility.
%
% The observer direction is supplied by the projection structure and is
% expressed in the same sky frame as proj.vertices3D.

V = proj.vertices3D;
F = proj.faces;

observerDirection = proj.observerDirection(:);
observerDirection = observerDirection / norm(observerDirection);

v1 = V(F(:,1),:);
v2 = V(F(:,2),:);
v3 = V(F(:,3),:);

% Facet centers
centers = (v1 + v2 + v3)/3;

% Oriented area vectors
areaVector = 0.5 * cross(v2-v1, v3-v1, 2);

% Facet areas
areas = vecnorm(areaVector,2,2);

if any(areas <= 0)
    error('Mesh contains zero-area facets.');
end

% Outward unit normals
normals = areaVector ./ areas;

% Cosine between outward normal and target-to-observer direction
cosView = normals * observerDirection;

visibilityTolerance = 1e-12;

isVisible = cosView > visibilityTolerance;

% Physical projected area as seen by observer
projectedAreas = areas .* max(cosView,0);

vis.faceCenters    = centers;
vis.faceNormals    = normals;
vis.faceAreas      = areas;
vis.cosView        = cosView;
vis.projectedAreas = projectedAreas;
vis.isVisible      = isVisible;

vis.observerDirection = observerDirection.';

end