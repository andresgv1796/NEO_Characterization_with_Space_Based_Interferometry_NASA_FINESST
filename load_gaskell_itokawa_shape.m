function shape = load_gaskell_itokawa_shape(tabFile)
%LOAD_GASKELL_ITOKAWA_SHAPE
%
% Load a Gaskell Itokawa vertex-facet shape model from the PDS
% verXXXq.tab format.
%
% The PDS vertex-facet file contains:
%
%   Line 1:
%
%       nVertices   nFaces
%
%   Next nVertices lines:
%
%       vertexID    x    y    z
%
%   Final nFaces lines:
%
%       faceID      vertex1    vertex2    vertex3
%
% The vertex coordinates are retained in the native Gaskell body-fixed
% frame:
%
%       +z_B : Itokawa rotation pole
%       +x_B : zero longitude
%
% No recentering or rotation is applied by this loader.
%
% OUTPUT
%
%   shape.vertices
%       N x 3 vertex coordinates [km]
%
%   shape.faces
%       M x 3 triangular facet connectivity
%
%   shape.type
%       'Gaskell Itokawa'
%
%   shape.params
%       Metadata and mesh diagnostics
%
%
% IMPORTANT
%
% The native reference frame is preserved deliberately. In particular,
% this function does NOT subtract mean(vertices) or otherwise force the
% shape origin to the geometric center.

arguments
    tabFile (1,:) char
end


%% Open file

fid = fopen(tabFile, 'r');

if fid < 0
    error('Could not open Gaskell shape file:\n%s', tabFile);
end

cleanupObj = onCleanup(@() fclose(fid));


%% Read model dimensions

counts = fscanf(fid, '%d', 2);

if numel(counts) ~= 2
    error('Could not read vertex and facet counts from %s.', tabFile);
end

nVertices = counts(1);
nFaces    = counts(2);

if nVertices <= 0 || nFaces <= 0
    error('Invalid vertex/facet counts in Gaskell shape file.');
end


%% Read vertex table
%
% Each row:
%
%   vertexID   x   y   z

vertexData = fscanf( ...
    fid, ...
    '%d %f %f %f', ...
    [4, nVertices]).';

if size(vertexData,1) ~= nVertices
    error([ ...
        'Expected %d vertex records but read %d.'], ...
        nVertices, ...
        size(vertexData,1));
end

vertexID = vertexData(:,1);

V = vertexData(:,2:4);


%% Validate vertex numbering

expectedVertexID = (1:nVertices).';

if ~isequal(vertexID, expectedVertexID)
    error([ ...
        'Vertex IDs are not sequential from 1 to nVertices. ', ...
        'An explicit ID-to-row mapping would be required.']);
end


%% Read triangular facet table
%
% Each row:
%
%   faceID   vertex1   vertex2   vertex3

faceData = fscanf( ...
    fid, ...
    '%d %d %d %d', ...
    [4, nFaces]).';

if size(faceData,1) ~= nFaces
    error([ ...
        'Expected %d facet records but read %d.'], ...
        nFaces, ...
        size(faceData,1));
end

faceID = faceData(:,1);

F = faceData(:,2:4);


%% Validate facet numbering

expectedFaceID = (1:nFaces).';

if ~isequal(faceID, expectedFaceID)
    error([ ...
        'Facet IDs are not sequential from 1 to nFaces.']);
end


%% Validate connectivity

if any(F(:) < 1) || any(F(:) > nVertices)
    error('Facet table contains invalid vertex indices.');
end


%% Remove vertices not referenced by any facet
%
% The PDS documentation notes that the vertex-facet conversion contains
% duplicate cube-edge/corner vertices that are not referenced by the facet
% table. They do not affect the surface, so remove only these UNUSED
% vertices.
%
% Do NOT merge vertices merely because they have identical coordinates.

usedVertex = false(nVertices,1);

usedVertex(unique(F(:))) = true;

nUnusedVertices = nnz(~usedVertex);


if nUnusedVertices > 0

    oldToNew = zeros(nVertices,1);

    oldToNew(usedVertex) = ...
        1:nnz(usedVertex);

    V = V(usedVertex,:);

    F = oldToNew(F);

end


%% Check for degenerate facets

v1 = V(F(:,1),:);
v2 = V(F(:,2),:);
v3 = V(F(:,3),:);

areaVector = ...
    0.5 * cross( ...
        v2 - v1, ...
        v3 - v1, ...
        2);

faceArea = ...
    vecnorm(areaVector,2,2);

areaTolerance = ...
    1e-14 * max(faceArea);

degenerateFacet = ...
    faceArea <= areaTolerance;

nDegenerateFaces = ...
    nnz(degenerateFacet);


if nDegenerateFaces > 0

    warning([ ...
        'Removing %d degenerate facets from Gaskell model.'], ...
        nDegenerateFaces);

    F = F(~degenerateFacet,:);

end


%% Signed polyhedron volume
%
% For outward-oriented triangular facets,
%
%   V = (1/6) sum r1 dot (r2 x r3)
%
% should be positive.

v1 = V(F(:,1),:);
v2 = V(F(:,2),:);
v3 = V(F(:,3),:);

signedSixVolume = ...
    dot( ...
        v1, ...
        cross(v2,v3,2), ...
        2);

signedVolume = ...
    sum(signedSixVolume)/6;


%% Correct global winding only if necessary

windingFlipped = false;

if signedVolume < 0

    warning([ ...
        'Gaskell facet winding was inward. ', ...
        'Reversing all triangle orientations.']);

    F(:,[2 3]) = F(:,[3 2]);

    windingFlipped = true;

    % Recompute with corrected orientation.

    v1 = V(F(:,1),:);
    v2 = V(F(:,2),:);
    v3 = V(F(:,3),:);

    signedSixVolume = ...
        dot( ...
            v1, ...
            cross(v2,v3,2), ...
            2);

    signedVolume = ...
        sum(signedSixVolume)/6;

end


if signedVolume <= 0
    error('Unable to obtain a positive enclosed polyhedron volume.');
end


%% Uniform-density polyhedron centroid
%
% Each oriented surface triangle and the origin form a signed tetrahedron.
%
% Tetrahedron centroid:
%
%   c_i = (r1 + r2 + r3)/4
%
% Weight by signed tetrahedral volume.

tetraVolume = ...
    signedSixVolume/6;

tetraCentroid = ...
    (v1 + v2 + v3)/4;

volumeCentroid = ...
    sum(tetraVolume .* tetraCentroid, 1) / ...
    signedVolume;


%% Surface area

v1 = V(F(:,1),:);
v2 = V(F(:,2),:);
v3 = V(F(:,3),:);

faceArea = ...
    0.5 * vecnorm( ...
        cross( ...
            v2 - v1, ...
            v3 - v1, ...
            2), ...
        2, ...
        2);

surfaceArea = ...
    sum(faceArea);


%% Bounding box

minXYZ = min(V,[],1);
maxXYZ = max(V,[],1);

dimensionsXYZ = ...
    maxXYZ - minXYZ;


%% Package shape

shape = struct();

shape.vertices = V;
shape.faces = F;

shape.type = 'Gaskell Itokawa';


%% Metadata

shape.params.source = ...
    'NASA PDS Gaskell Itokawa Shape Model';

shape.params.file = ...
    tabFile;

shape.params.frame = ...
    'Gaskell body-fixed';

shape.params.coordinateUnits = ...
    'km';

shape.params.nVerticesOriginal = ...
    nVertices;

shape.params.nFacesOriginal = ...
    nFaces;

shape.params.nVertices = ...
    size(V,1);

shape.params.nFaces = ...
    size(F,1);

shape.params.nUnusedVerticesRemoved = ...
    nUnusedVertices;

shape.params.nDegenerateFacesRemoved = ...
    nDegenerateFaces;

shape.params.windingFlipped = ...
    windingFlipped;

shape.params.signedVolume_km3 = ...
    signedVolume;

shape.params.surfaceArea_km2 = ...
    surfaceArea;

shape.params.volumeCentroid_km = ...
    volumeCentroid;

shape.params.boundingBoxMin_km = ...
    minXYZ;

shape.params.boundingBoxMax_km = ...
    maxXYZ;

shape.params.dimensionsXYZ_km = ...
    dimensionsXYZ;


%% Print summary

fprintf('\n');
fprintf('Gaskell Itokawa shape model\n');
fprintf('----------------------------\n');

fprintf('  file                       = %s\n', ...
    tabFile);

fprintf('  original vertices          = %d\n', ...
    nVertices);

fprintf('  original triangular facets = %d\n', ...
    nFaces);

fprintf('  unused vertices removed    = %d\n', ...
    nUnusedVertices);

fprintf('  degenerate facets removed  = %d\n', ...
    nDegenerateFaces);

fprintf('  retained vertices          = %d\n', ...
    size(V,1));

fprintf('  retained facets            = %d\n', ...
    size(F,1));

fprintf('  facet winding flipped      = %d\n', ...
    windingFlipped);

fprintf('\n');

fprintf('  dimensions [x y z]         = [%.6f %.6f %.6f] km\n', ...
    dimensionsXYZ);

fprintf('  enclosed volume            = %.8e km^3\n', ...
    signedVolume);

fprintf('  surface area               = %.8e km^2\n', ...
    surfaceArea);

fprintf('  volume centroid            = [%+.8e %+.8e %+.8e] km\n', ...
    volumeCentroid);

end