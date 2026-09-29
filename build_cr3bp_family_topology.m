function topology = build_cr3bp_family_topology( ...
    ICfamily,Tfamily,spectrumFeatures,cfg)
%BUILD_CR3BP_FAMILY_TOPOLOGY Recover continuation topology from unordered rows.
%
% JPL family tables can be ordered by a scalar invariant while containing
% several interleaved continuation branches.  Consequently, adjacent file
% rows are not assumed to be adjacent periodic orbits.  This routine builds
% a weighted complete graph in normalized state-period and Floquet-spectrum
% space, extracts its minimum spanning tree, and cuts only edges that are
% anomalously long relative to the local nearest-neighbor scale.
%
% Each connected component of the cut tree is an explicitly identified
% numerical continuation segment. It is a coverage/tracking domain, not an
% asserted physical subfamily. A richer within-segment neighbor graph is
% returned for invariant-center-subspace tracking. No row-order fallback is
% used.

    narginchk(4,4);

    validateattributes(ICfamily,{'numeric'}, ...
        {'real','finite','2d','ncols',6},mfilename,'ICfamily',1);
    validateattributes(Tfamily,{'numeric'}, ...
        {'real','finite','vector','positive'},mfilename,'Tfamily',2);
    validateattributes(spectrumFeatures,{'numeric'}, ...
        {'real','finite','2d'},mfilename,'spectrumFeatures',3);

    Tfamily = Tfamily(:);
    nRows = size(ICfamily,1);
    if numel(Tfamily) ~= nRows || size(spectrumFeatures,1) ~= nRows
        error('FINESST:FamilyTopology:DimensionMismatch', ...
            'ICfamily, Tfamily, and spectrumFeatures must have matching rows.');
    end

    required = {'topologyStateWeight','topologySpectrumWeight', ...
        'topologyEdgeRatioMax','topologyNeighborCount', ...
        'topologyMinSegmentRows'};
    if ~isstruct(cfg) || ~all(isfield(cfg,required))
        error('FINESST:FamilyTopology:InvalidConfiguration', ...
            'cfg must explicitly contain: %s.',strjoin(required,', '));
    end
    validateattributes(cfg.topologyStateWeight,{'numeric'}, ...
        {'real','finite','scalar','positive'});
    validateattributes(cfg.topologySpectrumWeight,{'numeric'}, ...
        {'real','finite','scalar','nonnegative'});
    validateattributes(cfg.topologyEdgeRatioMax,{'numeric'}, ...
        {'real','finite','scalar','>',1});
    validateattributes(cfg.topologyNeighborCount,{'numeric'}, ...
        {'real','finite','scalar','integer','>=',2,'<',nRows});
    validateattributes(cfg.topologyMinSegmentRows,{'numeric'}, ...
        {'real','finite','scalar','integer','positive'});

    statePeriod = [ICfamily,Tfamily];
    Zstate = local_normalize_columns(statePeriod);
    Zspectrum = local_normalize_columns(spectrumFeatures);
    Z = [cfg.topologyStateWeight*Zstate, ...
        cfg.topologySpectrumWeight*Zspectrum];

    squaredNorm = sum(Z.^2,2);
    D2 = max(0,squaredNorm+squaredNorm.'-2*(Z*Z.'));
    D = sqrt(D2);
    D(1:nRows+1:end) = inf;

    nearestDistance = min(D,[],2);
    if any(~isfinite(nearestDistance)) || any(nearestDistance <= 0)
        bad = find(~isfinite(nearestDistance) | nearestDistance <= 0,1);
        error('FINESST:FamilyTopology:DuplicateOrIsolatedDescriptor', ...
            ['Row %d has zero or nonfinite nearest-neighbor distance in ' ...
             'the state-period/Floquet descriptor.'],bad);
    end

    % The complete-graph MST is deterministic and does not depend on a
    % guessed file ordering or a k-nearest-neighbor connectivity fallback.
    [iUpper,jUpper] = find(triu(true(nRows),1));
    wUpper = D(sub2ind([nRows nRows],iUpper,jUpper));
    completeGraph = graph(iUpper,jUpper,wUpper,nRows);
    tree = minspantree(completeGraph,'Method','sparse');
    treeEnds = tree.Edges.EndNodes;
    treeWeight = tree.Edges.Weight;
    localScale = max(nearestDistance(treeEnds(:,1)), ...
        nearestDistance(treeEnds(:,2)));
    edgeRatio = treeWeight./localScale;
    keepTreeEdge = edgeRatio <= cfg.topologyEdgeRatioMax;

    cutTree = graph( ...
        treeEnds(keepTreeEdge,1),treeEnds(keepTreeEdge,2), ...
        treeWeight(keepTreeEdge),nRows);
    rawComponent = conncomp(cutTree).';
    rawIds = unique(rawComponent,'stable');

    % Stable component labels: largest component first, then lowest source row.
    componentSize = zeros(numel(rawIds),1);
    componentFirstRow = zeros(numel(rawIds),1);
    for k = 1:numel(rawIds)
        rows = find(rawComponent == rawIds(k));
        componentSize(k) = numel(rows);
        componentFirstRow(k) = min(rows);
    end
    [~,componentOrder] = sortrows([-componentSize,componentFirstRow],[1 2]);
    continuationSegmentId = zeros(nRows,1);
    for newId = 1:numel(componentOrder)
        oldId = rawIds(componentOrder(newId));
        continuationSegmentId(rawComponent == oldId) = newId;
    end

    nContinuationSegments = max(continuationSegmentId);
    segmentSize = accumarray(continuationSegmentId,1,[nContinuationSegments 1]);
    segmentAdmissible = segmentSize >= cfg.topologyMinSegmentRows;
    rowSegmentAdmissible = segmentAdmissible(continuationSegmentId);

    % A graph-geodesic coordinate provides an ordering diagnostic inside each
    % component without pretending that source-row order is continuation order.
    segmentCoordinate = nan(nRows,1);
    for iSegment = 1:nContinuationSegments
        rows = find(continuationSegmentId == iSegment);
        if numel(rows) == 1
            segmentCoordinate(rows) = 0;
            continue
        end
        Gsub = subgraph(cutTree,rows);
        d0 = distances(Gsub,1);
        [~,a] = max(d0);
        da = distances(Gsub,a);
        [diameterLength,~] = max(da);
        if ~(isfinite(diameterLength) && diameterLength > 0)
            error('FINESST:FamilyTopology:InvalidComponentDiameter', ...
                'Segment %d has invalid graph diameter.',iSegment);
        end
        segmentCoordinate(rows) = da(:)/diameterLength;
    end

    % Retain local graph edges for center-subspace continuation.  These edges
    % never cross a detected segment boundary.
    kNeighbor = cfg.topologyNeighborCount;
    nTreeRetained = nnz(keepTreeEdge);
    neighborPairs = zeros(nRows*kNeighbor+nTreeRetained,2);
    neighborWeight = zeros(nRows*kNeighbor+nTreeRetained,1);
    q = 0;
    for i = 1:nRows
        same = find(continuationSegmentId == continuationSegmentId(i));
        same(same == i) = [];
        [w,order] = sort(D(i,same),'ascend');
        nKeep = min(kNeighbor,numel(order));
        for k = 1:nKeep
            q = q+1;
            neighborPairs(q,:) = sort([i,same(order(k))]);
            neighborWeight(q) = w(k);
        end
    end

    % Always include the retained spanning-tree edges. The richer k-nearest
    % graph improves matching redundancy, while the tree guarantees that the
    % tracker sees a connected row graph inside every recovered component.
    retainedTreeEnds = treeEnds(keepTreeEdge,:);
    retainedTreeWeight = treeWeight(keepTreeEdge);
    for k = 1:nTreeRetained
        q = q+1;
        neighborPairs(q,:) = sort(retainedTreeEnds(k,:));
        neighborWeight(q) = retainedTreeWeight(k);
    end
    neighborPairs = neighborPairs(1:q,:);
    neighborWeight = neighborWeight(1:q);
    [neighborPairs,~,group] = unique(neighborPairs,'rows','sorted');
    neighborWeight = accumarray(group,neighborWeight,[],@min);

    neighborScale = max(nearestDistance(neighborPairs(:,1)), ...
        nearestDistance(neighborPairs(:,2)));
    neighborRatio = neighborWeight./neighborScale;
    keepNeighbor = neighborRatio <= cfg.topologyEdgeRatioMax;
    neighborPairs = neighborPairs(keepNeighbor,:);
    neighborWeight = neighborWeight(keepNeighbor);
    neighborRatio = neighborRatio(keepNeighbor);

    topology = struct();
    topology.descriptor = Z;
    topology.nearestNeighborDistance = nearestDistance;
    topology.continuationSegmentId = continuationSegmentId;
    topology.segmentCoordinate = segmentCoordinate;
    topology.nContinuationSegments = nContinuationSegments;
    topology.segmentSize = segmentSize;
    topology.segmentAdmissible = segmentAdmissible;
    topology.rowSegmentAdmissible = rowSegmentAdmissible;
    topology.treeEdges = treeEnds;
    topology.treeWeight = treeWeight;
    topology.treeEdgeRatio = edgeRatio;
    topology.treeEdgeRetained = keepTreeEdge;
    topology.neighborEdges = neighborPairs;
    topology.neighborWeight = neighborWeight;
    topology.neighborEdgeRatio = neighborRatio;
    topology.method = [ ...
        'complete-graph MST in normalized state-period/Floquet space; ' ...
        'local-scale edge cut; graph-geodesic component coordinate; ' ...
        'components are numerical coverage domains, not physical subfamilies'];
end


function Z = local_normalize_columns(A)
    columnMin = min(A,[],1);
    columnRange = max(A,[],1)-columnMin;
    active = columnRange > 100*eps(max(1,max(abs(A),[],1)));
    Z = zeros(size(A));
    Z(:,active) = (A(:,active)-columnMin(active))./columnRange(active);
    if any(active)
        Z = Z/sqrt(nnz(active));
    end
end
