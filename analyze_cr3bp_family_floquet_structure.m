function analysis = analyze_cr3bp_family_floquet_structure( ...
    ICfamily,Tfamily,mu,odeOpts,floquetOpts,cfg)
%ANALYZE_CR3BP_FAMILY_FLOQUET_STRUCTURE Analyze and sample one PO family.
%
% analysis = analyze_cr3bp_family_floquet_structure(...)
%
% Every row of the supplied periodic-orbit family is analyzed. Candidate
% parents are selected only after the family-wide Floquet structure is known.
% The selection contains:
%   1. the first and last admissible center-bearing family members;
%   2. farthest-point samples in combined continuation / multiplier-spectrum
%      space until the requested coverage is met or maxParents is reached.
% Stability transitions are retained as diagnostics and annotated whenever a
% selected sample lies on a transition boundary; they are not allowed to
% create an unbounded number of parents in a numerically noisy family.
%
% There are no fixed family fractions and no planar/vertical mode request.
% Every nontrivial center mode at every selected parent is returned.
% Source-row order is never interpreted as continuation order. A graph in
% normalized state-period/Floquet-spectrum space first identifies connected
% continuation segments, then tracks center invariant subspaces over graph
% edges. These segments are numerical coverage domains, not physical
% subfamily classifications.
%
% Required cfg fields
%   familyName
%   analysisVersion
%   stateClosureTol
%   branchTrackMaxCost
%   branchTrackUniquenessTol
%   signatureExcursionMaxRows
%   topologyStateWeight
%   topologySpectrumWeight
%   topologyEdgeRatioMax
%   topologyNeighborCount
%   topologyMinSegmentRows
%   minParentsPerSegment
%   maxParents
%   parentCoverageTol
%   continuationWeight
%   spectrumWeight
%   printProgress
%
% Existing dependencies
%   STM_vec, ode89, enumerate_cr3bp_center_modes

    narginchk(6,6);

    validateattributes(ICfamily,{'numeric'}, ...
        {'real','finite','2d','ncols',6},mfilename,'ICfamily',1);
    validateattributes(Tfamily,{'numeric'}, ...
        {'real','finite','vector','positive'},mfilename,'Tfamily',2);
    validateattributes(mu,{'numeric'}, ...
        {'real','finite','scalar','>',0,'<',0.5},mfilename,'mu',3);

    Tfamily = Tfamily(:);
    nRows = size(ICfamily,1);
    if numel(Tfamily) ~= nRows || nRows < 2
        error('FINESST:FamilyFloquet:InvalidFamilyDimensions', ...
            'ICfamily and Tfamily must describe at least two matching rows.');
    end

    if exist('ode89','file') ~= 2
        error('FINESST:FamilyFloquet:MissingOde89', ...
            'ode89 is required; no alternate solver is substituted.');
    end
    if ~isstruct(odeOpts)
        error('FINESST:FamilyFloquet:InvalidOdeOptions', ...
            'odeOpts must be an explicit ODE options structure.');
    end

    requiredCfg = { ...
        'familyName', ...
        'analysisVersion', ...
        'stateClosureTol', ...
        'branchTrackMaxCost', ...
        'branchTrackUniquenessTol', ...
        'signatureExcursionMaxRows', ...
        'topologyStateWeight', ...
        'topologySpectrumWeight', ...
        'topologyEdgeRatioMax', ...
        'topologyNeighborCount', ...
        'topologyMinSegmentRows', ...
        'minParentsPerSegment', ...
        'maxParents', ...
        'parentCoverageTol', ...
        'continuationWeight', ...
        'spectrumWeight', ...
        'printProgress'};

    if ~isstruct(cfg) || ~all(isfield(cfg,requiredCfg))
        error('FINESST:FamilyFloquet:InvalidConfiguration', ...
            'cfg must explicitly contain: %s.',strjoin(requiredCfg,', '));
    end

    familyName = string(cfg.familyName);
    if ~isscalar(familyName) || strlength(familyName) == 0
        error('FINESST:FamilyFloquet:InvalidFamilyName', ...
            'cfg.familyName must be a nonempty scalar string or character vector.');
    end

    validateattributes(cfg.stateClosureTol,{'numeric'}, ...
        {'real','finite','scalar','positive'});
    validateattributes(cfg.branchTrackMaxCost,{'numeric'}, ...
        {'real','finite','scalar','positive'});
    validateattributes(cfg.branchTrackUniquenessTol,{'numeric'}, ...
        {'real','finite','scalar','nonnegative'});
    validateattributes(cfg.signatureExcursionMaxRows,{'numeric'}, ...
        {'real','finite','scalar','integer','nonnegative'});
    validateattributes(cfg.minParentsPerSegment,{'numeric'}, ...
        {'real','finite','scalar','integer','positive'});
    validateattributes(cfg.maxParents,{'numeric'}, ...
        {'real','finite','scalar','integer','positive'});
    validateattributes(cfg.parentCoverageTol,{'numeric'}, ...
        {'real','finite','scalar','positive'});
    validateattributes(cfg.continuationWeight,{'numeric'}, ...
        {'real','finite','scalar','positive'});
    validateattributes(cfg.spectrumWeight,{'numeric'}, ...
        {'real','finite','scalar','positive'});

    if cfg.maxParents < 2
        error('FINESST:FamilyFloquet:InvalidParentLimits', ...
            'cfg.maxParents must be at least 2 to retain both family endpoints.');
    end
    if ~(islogical(cfg.printProgress) && isscalar(cfg.printProgress))
        error('FINESST:FamilyFloquet:InvalidPrintProgress', ...
            'cfg.printProgress must be a scalar logical.');
    end

    spectra = cell(nRows,1);
    monodromy = zeros(6,6,nRows);
    stateClosureNorm = nan(nRows,1);
    stateClosureAdmissible = false(nRows,1);

    Y0 = zeros(42,1);
    Y0(1:36) = reshape(eye(6),36,1);

    if cfg.printProgress
        fprintf('\nAnalyzing all %d family members for %s\n',nRows,familyName);
    end
    progressStride = max(1,floor(nRows/20));

    for iRow = 1:nRows
        if cfg.printProgress && (iRow == 1 || iRow == nRows || mod(iRow-1,progressStride) == 0)
            fprintf('  Floquet row %d / %d\n',iRow,nRows);
        end

        x0 = ICfamily(iRow,:).';
        T = Tfamily(iRow);
        Y0(37:42) = x0;

        [tSol,YSol] = ode89(@(t,Y) STM_vec(t,Y,mu),[0 T],Y0,odeOpts);

        if numel(tSol) < 2 || size(YSol,2) ~= 42 || any(~isfinite(YSol(:)))
            error('FINESST:FamilyFloquet:InvalidPropagation', ...
                'Invalid state/STM propagation at family row %d.',iRow);
        end

        PhiT = reshape(YSol(end,1:36),6,6);
        xT = YSol(end,37:42).';
        closureNorm = norm(xT-x0);

        monodromy(:,:,iRow) = PhiT;
        stateClosureNorm(iRow) = closureNorm;
        stateClosureAdmissible(iRow) = closureNorm <= cfg.stateClosureTol;
        spectra{iRow} = enumerate_cr3bp_center_modes( ...
            PhiT,T,x0,mu,floquetOpts);
    end

    centerCount = zeros(nRows,1);
    nUnstableMultipliers = zeros(nRows,1);
    spectralRadius = zeros(nRows,1);
    stabilitySignature = strings(nRows,1);
    trivialPairResolved = false(nRows,1);
    trivialPairMaxDistance = nan(nRows,1);
    trivialFlowAlignment = nan(nRows,1);
    complexQuartetDetected = false(nRows,1);
    nFeatures = numel(spectra{1}.featureVector);
    spectrumFeatures = zeros(nRows,nFeatures);

    for iRow = 1:nRows
        centerCount(iRow) = spectra{iRow}.nCenterModes;
        nUnstableMultipliers(iRow) = spectra{iRow}.nUnstableMultipliers;
        spectralRadius(iRow) = spectra{iRow}.spectralRadius;
        stabilitySignature(iRow) = string(spectra{iRow}.stabilitySignature);
        trivialPairResolved(iRow) = spectra{iRow}.trivialPairResolved;
        trivialPairMaxDistance(iRow) = ...
            max(spectra{iRow}.trivialMultiplierDistance);
        trivialFlowAlignment(iRow) = spectra{iRow}.trivialFlowAlignment;
        complexQuartetDetected(iRow) = spectra{iRow}.complexQuartetDetected;
        spectrumFeatures(iRow,:) = spectra{iRow}.featureVector;
    end

    topology = build_cr3bp_family_topology( ...
        ICfamily,Tfamily,spectrumFeatures,cfg);
    continuationSegmentId = topology.continuationSegmentId;
    segmentCoordinate = topology.segmentCoordinate;

    % Center branches are connected over the recovered family-neighbor graph,
    % never over adjacent CSV rows. Connected components of unambiguous
    % invariant-subspace matches define persistent tracking IDs; neither the
    % IDs nor graph segments assert a physical subfamily taxonomy.
    [spectra,branchDiagnostics] = local_track_modes_on_family_graph( ...
        spectra,topology,cfg);

    % Signature filtering and transition detection are performed separately
    % inside each recovered segment. Interleaved source rows can therefore
    % never create a fictitious stability transition.
    selectionSignature = stabilitySignature;
    transitionBoundary = false(nRows,1);
    for iSegment = 1:topology.nContinuationSegments
        rows = find(continuationSegmentId == iSegment);
        [~,order] = sort(segmentCoordinate(rows),'ascend');
        rows = rows(order);
        filtered = local_remove_short_signature_excursions( ...
            stabilitySignature(rows),cfg.signatureExcursionMaxRows);
        selectionSignature(rows) = filtered;
        changed = filtered(2:end) ~= filtered(1:end-1);
        kChange = find(changed);
        transitionBoundary(rows(kChange)) = true;
        transitionBoundary(rows(kChange+1)) = true;
    end

    eligible = stateClosureAdmissible & trivialPairResolved & ...
        centerCount > 0 & topology.rowSegmentAdmissible;
    eligibleRows = find(eligible);
    selectionReason = strings(nRows,1);

    % Every center-bearing numerical segment is represented independently by
    % its accessible endpoints. This is a coverage rule, not a physical-family
    % classification and not a row-fraction rule.
    for iSegment = 1:topology.nContinuationSegments
        rows = find(eligible & continuationSegmentId == iSegment);
        if isempty(rows)
            continue
        end
        [~,iFirst] = min(segmentCoordinate(rows));
        [~,iLast] = max(segmentCoordinate(rows));
        firstRow = rows(iFirst);
        lastRow = rows(iLast);
        selectionReason(firstRow) = local_add_reason( ...
            selectionReason(firstRow),sprintf( ...
            'segment %d first center-bearing member',iSegment));
        selectionReason(lastRow) = local_add_reason( ...
            selectionReason(lastRow),sprintf( ...
            'segment %d last center-bearing member',iSegment));
    end

    selectedRows = find(strlength(selectionReason) > 0);
    if numel(selectedRows) > cfg.maxParents
        error('FINESST:FamilyFloquet:ParentBudgetBelowTopologyMinimum', ...
            ['The %d mandatory segment endpoints exceed maxParents=%d. ' ...
             'Increase maxParents; no segment endpoint is silently dropped.'], ...
            numel(selectedRows),cfg.maxParents);
    end

    featureMin = min(spectrumFeatures,[],1);
    featureRange = max(spectrumFeatures,[],1)-featureMin;
    activeFeature = featureRange > ...
        100*eps(max(1,max(abs(spectrumFeatures),[],1)));
    Fnormalized = zeros(size(spectrumFeatures));
    Fnormalized(:,activeFeature) = ...
        (spectrumFeatures(:,activeFeature)-featureMin(activeFeature)) ./ ...
        featureRange(activeFeature);
    if any(activeFeature)
        Fnormalized = Fnormalized/sqrt(nnz(activeFeature));
    end
    Zcoverage = [cfg.continuationWeight*segmentCoordinate, ...
        cfg.spectrumWeight*Fnormalized];

    % Meet the minimum representation independently in each segment.
    for iSegment = 1:topology.nContinuationSegments
        rows = find(eligible & continuationSegmentId == iSegment);
        nRequired = min(cfg.minParentsPerSegment,numel(rows));
        while nnz(ismember(rows,selectedRows)) < nRequired
            row = local_farthest_unselected_row(rows,selectedRows,Zcoverage);
            selectedRows(end+1,1) = row; %#ok<AGROW>
            selectionReason(row) = local_add_reason( ...
                selectionReason(row),"segment minimum coverage");
        end
    end

    % Continue bounded farthest-point coverage, always measuring distance to
    % representatives from the same numerical continuation segment.
    while numel(selectedRows) < cfg.maxParents
        bestDistance = -inf;
        bestRow = NaN;
        for iSegment = 1:topology.nContinuationSegments
            rows = find(eligible & continuationSegmentId == iSegment);
            if isempty(rows)
                continue
            end
            selectedSegment = selectedRows( ...
                continuationSegmentId(selectedRows) == iSegment);
            [distance,row] = local_farthest_distance( ...
                rows,selectedSegment,Zcoverage);
            if distance > bestDistance
                bestDistance = distance;
                bestRow = row;
            end
        end
        if ~(isfinite(bestDistance) && bestDistance > cfg.parentCoverageTol)
            break
        end
        selectedRows(end+1,1) = bestRow; %#ok<AGROW>
        if transitionBoundary(bestRow)
            reason = "Floquet-space coverage; transition vicinity";
        else
            reason = "Floquet-space coverage";
        end
        selectionReason(bestRow) = local_add_reason( ...
            selectionReason(bestRow),reason);
    end

    selectedRows = unique(selectedRows(:),'sorted');

    emptyParent = struct( ...
        'row',{}, ...
        'continuationSegmentId',{}, ...
        'segmentCoordinate',{}, ...
        'state',{}, ...
        'TPrimary',{}, ...
        'selectionReason',{}, ...
        'spectrum',{}, ...
        'centerModes',{});
    candidateParents = emptyParent;

    for k = 1:numel(selectedRows)
        row = selectedRows(k);
        parent = struct();
        parent.row = row;
        parent.continuationSegmentId = continuationSegmentId(row);
        parent.segmentCoordinate = segmentCoordinate(row);
        parent.state = ICfamily(row,:).';
        parent.TPrimary = Tfamily(row);
        parent.selectionReason = char(selectionReason(row));
        parent.spectrum = spectra{row};
        parent.centerModes = spectra{row}.centerModes;
        candidateParents(end+1) = parent; %#ok<AGROW>
    end

    selectedMask = false(nRows,1);
    selectedMask(selectedRows) = true;
    rowStatus = repmat("admissible",nRows,1);
    rowStatus(~topology.rowSegmentAdmissible) = ...
        "topology-segment-too-short";
    rowStatus(~stateClosureAdmissible) = "poor-state-closure";
    rowStatus(stateClosureAdmissible & ~trivialPairResolved) = ...
        "trivial-pair-quality-warning";
    rowStatus(stateClosureAdmissible & trivialPairResolved & centerCount == 0) = ...
        "no-nontrivial-center-mode";

    scanTable = table( ...
        (1:nRows).',continuationSegmentId,segmentCoordinate, ...
        topology.rowSegmentAdmissible,Tfamily,stateClosureNorm, ...
        stateClosureAdmissible,trivialPairResolved,trivialPairMaxDistance, ...
        trivialFlowAlignment,centerCount,nUnstableMultipliers, ...
        complexQuartetDetected,spectralRadius,stabilitySignature, ...
        selectionSignature,rowStatus,selectedMask,selectionReason, ...
        'VariableNames',{ ...
            'FamilyRow','ContinuationSegmentId','SegmentCoordinate', ...
            'TopologySegmentAdmissible','TPrimary_TU','StateClosureNorm', ...
            'StateClosureAdmissible','TrivialPairResolved', ...
            'TrivialPairMaxDistance','TrivialFlowAlignment', ...
            'CenterModeCount','UnstableMultiplierCount','ComplexQuartetDetected', ...
            'SpectralRadius','StabilitySignature','SelectionSignature','RowStatus', ...
            'SelectedParent','SelectionReason'});

    modeRows = cell(sum(centerCount),1);
    q = 0;
    for iRow = 1:nRows
        modes = spectra{iRow}.centerModes;
        for j = 1:numel(modes)
            q = q + 1;
            modeRows{q} = table( ...
                iRow,continuationSegmentId(iRow),segmentCoordinate(iRow), ...
                modes(j).localIndex,modes(j).segmentBranchId, ...
                modes(j).branchId,real(modes(j).lambda),imag(modes(j).lambda), ...
                modes(j).rho,abs(modes(j).lambda)-1, ...
                modes(j).radialLogDefect,modes(j).conjugacyDefect, ...
                modes(j).eigenpairResidualRel, ...
                modes(j).planarProjection,modes(j).verticalProjection, ...
                modes(j).positionRankRatioAtSection,selectedMask(iRow), ...
                'VariableNames',{ ...
                    'FamilyRow','ContinuationSegmentId','SegmentCoordinate', ...
                    'LocalModeIndex','SegmentBranchId','BranchId', ...
                    'LambdaReal','LambdaImag','Rho_rad','SignedUnitCircleDefect', ...
                    'RadialLogDefect','ConjugacyDefect', ...
                    'EigenpairResidualRel','PlanarProjection','VerticalProjection', ...
                    'PositionRankRatioAtSection','SelectedParent'});
        end
    end

    if q > 0
        modeTable = vertcat(modeRows{1:q});
    else
        modeTable = table();
    end

    analysis = struct();
    analysis.familyName = char(familyName);
    analysis.analysisVersion = char(string(cfg.analysisVersion));
    analysis.ICfamily = ICfamily;
    analysis.Tfamily = Tfamily;
    analysis.topology = topology;
    analysis.continuationSegmentId = continuationSegmentId;
    analysis.segmentCoordinate = segmentCoordinate;
    analysis.branchDiagnostics = branchDiagnostics;
    analysis.monodromy = monodromy;
    analysis.spectra = spectra;
    analysis.stateClosureNorm = stateClosureNorm;
    analysis.stateClosureAdmissible = stateClosureAdmissible;
    analysis.trivialPairResolved = trivialPairResolved;
    analysis.candidateParents = candidateParents;
    analysis.selectedRows = selectedRows;
    analysis.scanTable = scanTable;
    analysis.modeTable = modeTable;
    analysis.floquetOptions = floquetOpts;
    analysis.configuration = cfg;
    analysis.selectionMethod = ...
        ['all-row Floquet analysis; graph-recovered numerical coverage segments; ' ...
         'graph-connected center invariant subspaces; segment endpoints + ' ...
         'bounded farthest-point continuation/spectrum coverage'];
end


function [spectra,diagnostics] = local_track_modes_on_family_graph( ...
        spectra,topology,cfg)
    nRows = numel(spectra);
    modeCount = cellfun(@(s) numel(s.centerModes),spectra);
    nInstances = sum(modeCount);
    instanceIndex = cell(nRows,1);
    instanceRow = zeros(nInstances,1);
    instanceLocal = zeros(nInstances,1);
    q = 0;
    for iRow = 1:nRows
        nMode = modeCount(iRow);
        instanceIndex{iRow} = q+(1:nMode);
        for j = 1:nMode
            q = q+1;
            instanceRow(q) = iRow;
            instanceLocal(q) = j;
        end
    end

    if nInstances == 0
        diagnostics = struct();
        diagnostics.nModeInstances = 0;
        diagnostics.nMatchEdges = 0;
        diagnostics.nBranches = 0;
        diagnostics.nBranchesPerSegment = zeros(topology.nContinuationSegments,1);
        diagnostics.ambiguousMatchCount = 0;
        diagnostics.method = ...
            'no center-mode instances were available for graph tracking';
        return
    end

    edgeA = zeros(4*size(topology.neighborEdges,1),1);
    edgeB = zeros(size(edgeA));
    qEdge = 0;
    ambiguousMatchCount = 0;

    for iEdge = 1:size(topology.neighborEdges,1)
        rowA = topology.neighborEdges(iEdge,1);
        rowB = topology.neighborEdges(iEdge,2);
        modesA = spectra{rowA}.centerModes;
        modesB = spectra{rowB}.centerModes;
        if isempty(modesA) || isempty(modesB)
            continue
        end

        cost = zeros(numel(modesA),numel(modesB));
        for i = 1:numel(modesA)
            for j = 1:numel(modesB)
                overlap = norm(modesA(i).basis6.'*modesB(j).basis6,'fro')/sqrt(2);
                overlap = min(1,max(0,overlap));
                rhoDistance = abs(modesA(i).rho-modesB(j).rho)/pi;
                cost(i,j) = rhoDistance+1-overlap;
            end
        end

        pairList = zeros(numel(cost),3);
        p = 0;
        for i = 1:size(cost,1)
            for j = 1:size(cost,2)
                p = p+1;
                pairList(p,:) = [cost(i,j),i,j];
            end
        end
        pairList = sortrows(pairList,1);
        usedA = false(numel(modesA),1);
        usedB = false(numel(modesB),1);

        for p = 1:size(pairList,1)
            c = pairList(p,1);
            i = pairList(p,2);
            j = pairList(p,3);
            if c > cfg.branchTrackMaxCost
                break
            end
            if usedA(i) || usedB(j)
                continue
            end
            if c > min(cost(i,:))+10*eps(max(1,c)) || ...
                    c > min(cost(:,j))+10*eps(max(1,c))
                continue
            end

            rowAlternative = sort(cost(i,:),'ascend');
            columnAlternative = sort(cost(:,j),'ascend');
            if numel(rowAlternative) > 1
                rowMargin = rowAlternative(2)-rowAlternative(1);
            else
                rowMargin = inf;
            end
            if numel(columnAlternative) > 1
                columnMargin = columnAlternative(2)-columnAlternative(1);
            else
                columnMargin = inf;
            end
            if rowMargin <= cfg.branchTrackUniquenessTol || ...
                    columnMargin <= cfg.branchTrackUniquenessTol
                ambiguousMatchCount = ambiguousMatchCount+1;
                continue
            end

            qEdge = qEdge+1;
            edgeA(qEdge) = instanceIndex{rowA}(i);
            edgeB(qEdge) = instanceIndex{rowB}(j);
            usedA(i) = true;
            usedB(j) = true;
        end
    end

    edgeA = edgeA(1:qEdge);
    edgeB = edgeB(1:qEdge);
    modeGraph = graph(edgeA,edgeB,[],nInstances);
    rawBranch = conncomp(modeGraph).';
    rawIds = unique(rawBranch,'stable');

    branchSegment = zeros(numel(rawIds),1);
    branchMedianRho = zeros(numel(rawIds),1);
    branchFirstRow = zeros(numel(rawIds),1);
    for k = 1:numel(rawIds)
        instances = find(rawBranch == rawIds(k));
        rows = instanceRow(instances);
        segments = unique(topology.continuationSegmentId(rows));
        if numel(segments) ~= 1
            error('FINESST:FamilyFloquet:BranchCrossesContinuationSegments', ...
                'A tracked center branch crosses recovered segment boundaries.');
        end
        branchSegment(k) = segments;
        rho = zeros(numel(instances),1);
        for j = 1:numel(instances)
            rho(j) = spectra{rows(j)}.centerModes( ...
                instanceLocal(instances(j))).rho;
        end
        branchMedianRho(k) = median(rho);
        branchFirstRow(k) = min(rows);
    end

    [~,branchOrder] = sortrows( ...
        [branchSegment,branchMedianRho,branchFirstRow],[1 2 3]);
    globalBranchId = zeros(numel(rawIds),1);
    localBranchId = zeros(numel(rawIds),1);
    nextLocal = zeros(topology.nContinuationSegments,1);
    for newId = 1:numel(branchOrder)
        old = branchOrder(newId);
        globalBranchId(old) = newId;
        iSegment = branchSegment(old);
        nextLocal(iSegment) = nextLocal(iSegment)+1;
        localBranchId(old) = nextLocal(iSegment);
    end

    for iInstance = 1:nInstances
        raw = find(rawIds == rawBranch(iInstance),1);
        row = instanceRow(iInstance);
        local = instanceLocal(iInstance);
        modes = spectra{row}.centerModes;
        modes(local).continuationSegmentId = topology.continuationSegmentId(row);
        modes(local).segmentBranchId = localBranchId(raw);
        modes(local).branchId = globalBranchId(raw);
        spectra{row}.centerModes = modes;
    end

    diagnostics = struct();
    diagnostics.nModeInstances = nInstances;
    diagnostics.nMatchEdges = qEdge;
    diagnostics.nBranches = numel(rawIds);
    diagnostics.nBranchesPerSegment = nextLocal;
    diagnostics.ambiguousMatchCount = ambiguousMatchCount;
    diagnostics.method = [ ...
        'connected components of unique invariant-subspace matches over ' ...
        'the recovered family-neighbor graph'];
end


function row = local_farthest_unselected_row(rows,selectedRows,Z)
    [~,row] = local_farthest_distance(rows,selectedRows,Z);
end


function [maxDistance,row] = local_farthest_distance(rows,selectedRows,Z)
    rows = rows(:);
    available = ~ismember(rows,selectedRows);
    if ~any(available)
        maxDistance = 0;
        row = NaN;
        return
    end
    if isempty(selectedRows)
        error('FINESST:FamilyFloquet:MissingSegmentSeed', ...
            'Farthest-point coverage requires a selected seed in each segment.');
    end
    dmin = inf(numel(rows),1);
    for k = 1:numel(selectedRows)
        delta = Z(rows,:)-Z(selectedRows(k),:);
        dmin = min(dmin,sqrt(sum(delta.^2,2)));
    end
    dmin(~available) = -inf;
    [maxDistance,i] = max(dmin);
    row = rows(i);
end


function out = local_add_reason(existing,newReason)
    if strlength(existing) == 0
        out = newReason;
    elseif contains(existing,newReason)
        out = existing;
    else
        out = existing + "; " + newReason;
    end
end


function filtered = local_remove_short_signature_excursions(raw,maxRows)
    filtered = raw;
    if maxRows == 0 || numel(raw) < 3
        return
    end

    changed = true;
    while changed
        changed = false;
        runStart = [1;find(filtered(2:end) ~= filtered(1:end-1))+1];
        runEnd = [runStart(2:end)-1;numel(filtered)];

        for k = 2:numel(runStart)-1
            runLength = runEnd(k)-runStart(k)+1;
            if runLength <= maxRows && ...
                    filtered(runStart(k)-1) == filtered(runEnd(k)+1)
                filtered(runStart(k):runEnd(k)) = filtered(runStart(k)-1);
                changed = true;
            end
        end
    end
end
