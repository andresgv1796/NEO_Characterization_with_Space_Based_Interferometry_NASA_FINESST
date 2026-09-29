function closure = compute_array_closure_phases( ...
    obs, baselinePairs, minVisibilityAmplitude, referenceCollector)
%COMPUTE_ARRAY_CLOSURE_PHASES
%
% Compute closure phases for every available collector triangle in a
% general N-collector interferometric array.
%
% A triangle formed by collectors (i,j,k) satisfies
%
%   B_ij + B_jk + B_ki = 0
%
% and therefore has bispectrum
%
%   T_ijk = V_ij .* V_jk .* V_ki
%
% and closure phase
%
%   Phi_ijk = angle(T_ijk).
%
% When only the increasing-index baseline V_ik is stored,
%
%   V_ki = conj(V_ik),
%
% so this becomes
%
%   T_ijk = V_ij .* V_jk .* conj(V_ik).
%
% INPUTS
%
%   obs
%       Observation structure containing
%
%           obs.baselines(q).V
%
%       for q = 1,...,N_b.
%
%       Every visibility vector must use the same sample ordering and have
%       the same number of samples. Missing measurements may be represented
%       by complex NaN values.
%
%   baselinePairs
%       N_b-by-2 array giving the directed collector pair corresponding to
%       each element of obs.baselines.
%
%       Example for three collectors:
%
%           baselinePairs = [
%               1 2
%               1 3
%               2 3
%           ];
%
%       Row q corresponds to obs.baselines(q).V.
%
%       The orientation matters. For example, [2 1] means that the stored
%       visibility corresponds to
%
%           B_21 = r_1 - r_2.
%
%       The function automatically conjugates a visibility when the
%       opposite orientation is required.
%
%   minVisibilityAmplitude
%       Minimum required visibility amplitude on all three triangle
%       baselines. A closure phase is marked invalid when any participating
%       baseline satisfies
%
%           abs(V_ij) < minVisibilityAmplitude.
%
%   referenceCollector
%       Collector used to construct a nonredundant reference-triangle
%       closure-phase set.
%
%       For a complete array, choosing referenceCollector = 1 produces
%
%           (1,2,3), (1,2,4), ..., (1,N-1,N),
%
%       with
%
%           (N_c-1)(N_c-2)/2
%
%       independent closure phases.
%
%       Pass [] when no reference subset is needed.
%
% OUTPUT
%
%   closure.collectorIds
%       Collector identifiers found in baselinePairs.
%
%   closure.baselinePairs
%       Original directed baseline-pair list.
%
%   closure.triangles
%       N_tri-by-3 list of all available collector triangles.
%
%   closure.baselineIndices
%       N_tri-by-3 indices into obs.baselines for the directed edges
%
%           i -> j, j -> k, k -> i.
%
%   closure.baselineOrientation
%       N_tri-by-3 orientation factors.
%
%       +1 means the stored baseline has the requested direction.
%       -1 means the stored visibility was conjugated.
%
%   closure.bispectrum
%       N_sample-by-N_tri complex bispectrum array.
%
%   closure.bispectrumAmplitude
%       Absolute bispectrum amplitude.
%
%   closure.phase_rad
%       Closure phases in radians, wrapped to [-pi,pi].
%
%   closure.phase_deg
%       Closure phases in degrees, wrapped to [-180,180].
%
%   closure.valid
%       Logical validity mask with size N_sample-by-N_tri.
%
%   closure.reference
%       Reference-collector subset and completeness diagnostics.
%
% NOTES
%
%   1. baselinePairs must contain physical baselines only. Do not pass the
%      Hermitian samples added during dirty-image formation.
%
%   2. The visibility arrays must be aligned sample-by-sample. Element k of
%      every baseline visibility vector must correspond to the same
%      observation time.
%
%   3. The function permits incomplete arrays. A triangle is returned only
%      when all three required physical baselines are available.
%
%   4. A complete N-collector array has
%
%          N_b   = N_c(N_c-1)/2
%          N_tri = N_c(N_c-1)(N_c-2)/6
%
%      triangle closure phases, but only
%
%          N_ind = (N_c-1)(N_c-2)/2
%
%      independent closure phases.

    %% Validate primary inputs

    if ~isstruct(obs) || ~isfield(obs, 'baselines')
        error('obs must contain the field obs.baselines.');
    end

    nBaselines = numel(obs.baselines);

    if nBaselines == 0
        error('obs.baselines is empty.');
    end

    validateattributes( ...
        baselinePairs, ...
        {'numeric'}, ...
        {'2d', 'ncols', 2, 'integer', 'positive', 'finite'}, ...
        mfilename, ...
        'baselinePairs');

    if size(baselinePairs, 1) ~= nBaselines
        error(['baselinePairs must contain one row for every element of ', ...
            'obs.baselines.']);
    end

    if any(baselinePairs(:,1) == baselinePairs(:,2))
        error('A baseline cannot connect a collector to itself.');
    end

    validateattributes( ...
        minVisibilityAmplitude, ...
        {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'}, ...
        mfilename, ...
        'minVisibilityAmplitude');

    if ~isempty(referenceCollector)
        validateattributes( ...
            referenceCollector, ...
            {'numeric'}, ...
            {'scalar', 'integer', 'positive', 'finite'}, ...
            mfilename, ...
            'referenceCollector');
    end


    %% Check for duplicate physical baselines

    undirectedPairs = sort(baselinePairs, 2);
    uniqueUndirectedPairs = unique(undirectedPairs, 'rows');

    if size(uniqueUndirectedPairs, 1) ~= nBaselines
        error(['baselinePairs contains duplicate physical baselines. ', ...
            'Each collector pair must appear exactly once.']);
    end


    %% Assemble the visibility matrix

    if ~isfield(obs.baselines(1), 'V')
        error('Every obs.baselines element must contain a field named V.');
    end

    nSamples = numel(obs.baselines(1).V);

    if nSamples == 0
        error('The visibility vectors are empty.');
    end

    visibilityMatrix = complex(zeros(nSamples, nBaselines));

    for q = 1:nBaselines

        if ~isfield(obs.baselines(q), 'V')
            error('Every obs.baselines element must contain a field named V.');
        end

        Vq = obs.baselines(q).V(:);

        if numel(Vq) ~= nSamples
            error(['All baseline visibility arrays must have the same ', ...
                'number of samples and use the same sample ordering.']);
        end

        visibilityMatrix(:,q) = Vq;
    end


    %% Identify the collectors

    collectorIds = unique(baselinePairs(:)).';
    nCollectors = numel(collectorIds);

    if nCollectors < 3
        error('At least three distinct collectors are required.');
    end

    if ~isempty(referenceCollector) && ...
            ~ismember(referenceCollector, collectorIds)

        error('referenceCollector is not present in baselinePairs.');
    end


    %% Generate all possible collector triangles

    candidateTriangles = nchoosek(collectorIds, 3);
    nCandidateTriangles = size(candidateTriangles, 1);

    triangles = zeros(nCandidateTriangles, 3);
    baselineIndices = zeros(nCandidateTriangles, 3);
    baselineOrientation = zeros(nCandidateTriangles, 3);

    bispectrum = complex(zeros(nSamples, nCandidateTriangles));
    valid = false(nSamples, nCandidateTriangles);

    nTriangles = 0;


    %% Construct every available triangle closure

    for qTriangle = 1:nCandidateTriangles

        i = candidateTriangles(qTriangle, 1);
        j = candidateTriangles(qTriangle, 2);
        k = candidateTriangles(qTriangle, 3);

        % Directed triangle:
        %
        %   i -> j
        %   j -> k
        %   k -> i
        %
        % so that
        %
        %   B_ij + B_jk + B_ki = 0.

        [Vij, idxij, signij, foundij] = ...
            get_directed_visibility( ...
                visibilityMatrix, baselinePairs, i, j);

        [Vjk, idxjk, signjk, foundjk] = ...
            get_directed_visibility( ...
                visibilityMatrix, baselinePairs, j, k);

        [Vki, idxki, signki, foundki] = ...
            get_directed_visibility( ...
                visibilityMatrix, baselinePairs, k, i);

        % Skip incomplete triangles.

        if ~(foundij && foundjk && foundki)
            continue
        end

        nTriangles = nTriangles + 1;

        triangles(nTriangles,:) = [i, j, k];

        baselineIndices(nTriangles,:) = [
            idxij
            idxjk
            idxki
        ].';

        baselineOrientation(nTriangles,:) = [
            signij
            signjk
            signki
        ].';

        Tijk = Vij .* Vjk .* Vki;

        finiteTriangle = ...
            isfinite(real(Vij)) & isfinite(imag(Vij)) & ...
            isfinite(real(Vjk)) & isfinite(imag(Vjk)) & ...
            isfinite(real(Vki)) & isfinite(imag(Vki));

        amplitudeValid = ...
            abs(Vij) >= minVisibilityAmplitude & ...
            abs(Vjk) >= minVisibilityAmplitude & ...
            abs(Vki) >= minVisibilityAmplitude;

        validTriangle = finiteTriangle & amplitudeValid;

        bispectrum(:,nTriangles) = Tijk;
        valid(:,nTriangles) = validTriangle;
    end


    %% Remove unused preallocated columns

    if nTriangles == 0
        error(['No complete three-collector triangles could be formed ', ...
            'from baselinePairs.']);
    end

    triangles = triangles(1:nTriangles,:);
    baselineIndices = baselineIndices(1:nTriangles,:);
    baselineOrientation = baselineOrientation(1:nTriangles,:);

    bispectrum = bispectrum(:,1:nTriangles);
    valid = valid(:,1:nTriangles);


    %% Compute closure phase

    phaseRad = NaN(size(bispectrum));

    phaseRad(valid) = angle(bispectrum(valid));

    phaseDeg = rad2deg(phaseRad);


    %% Determine complete-array properties

    nCompleteBaselines = nCollectors*(nCollectors - 1)/2;
    nCompleteTriangles = nchoosek(nCollectors, 3);
    nIndependentClosurePhases = ...
        (nCollectors - 1)*(nCollectors - 2)/2;

    isCompletePhysicalArray = ...
        nBaselines == nCompleteBaselines;

    allTrianglesAvailable = ...
        nTriangles == nCompleteTriangles;


    %% Reference-collector closure basis

    reference = struct();

    reference.collector = referenceCollector;

    if isempty(referenceCollector)

        reference.triangleMask = false(nTriangles, 1);
        reference.triangleIndices = zeros(0,1);
        reference.triangles = zeros(0,3);

        reference.bispectrum = complex(zeros(nSamples,0));
        reference.bispectrumAmplitude = zeros(nSamples,0);

        reference.phase_rad = NaN(nSamples,0);
        reference.phase_deg = NaN(nSamples,0);
        reference.valid = false(nSamples,0);

        reference.expectedCount = 0;
        reference.availableCount = 0;
        reference.isCompleteBasis = false;

    else

        referenceMask = ...
            any(triangles == referenceCollector, 2);

        referenceIndices = find(referenceMask);

        referenceExpectedCount = ...
            (nCollectors - 1)*(nCollectors - 2)/2;

        reference.triangleMask = referenceMask;
        reference.triangleIndices = referenceIndices;
        reference.triangles = triangles(referenceMask,:);

        reference.bispectrum = bispectrum(:,referenceMask);
        reference.bispectrumAmplitude = ...
            abs(bispectrum(:,referenceMask));

        reference.phase_rad = phaseRad(:,referenceMask);
        reference.phase_deg = phaseDeg(:,referenceMask);
        reference.valid = valid(:,referenceMask);

        reference.expectedCount = referenceExpectedCount;
        reference.availableCount = nnz(referenceMask);

        reference.isCompleteBasis = ...
            isCompletePhysicalArray && ...
            allTrianglesAvailable && ...
            nnz(referenceMask) == referenceExpectedCount;
    end


    %% Output structure

    closure = struct();

    closure.collectorIds = collectorIds;
    closure.nCollectors = nCollectors;

    closure.baselinePairs = baselinePairs;
    closure.nPhysicalBaselines = nBaselines;

    closure.triangles = triangles;
    closure.nTriangles = nTriangles;

    closure.baselineIndices = baselineIndices;
    closure.baselineOrientation = baselineOrientation;

    closure.bispectrum = bispectrum;
    closure.bispectrumAmplitude = abs(bispectrum);

    closure.phase_rad = phaseRad;
    closure.phase_deg = phaseDeg;

    closure.valid = valid;

    closure.minVisibilityAmplitude = minVisibilityAmplitude;

    closure.isCompletePhysicalArray = isCompletePhysicalArray;
    closure.allTrianglesAvailable = allTrianglesAvailable;

    closure.nCompleteBaselines = nCompleteBaselines;
    closure.nCompleteTriangles = nCompleteTriangles;
    closure.nIndependentClosurePhases = ...
        nIndependentClosurePhases;

    closure.reference = reference;

    if isfield(obs, 'accessibleIndex')
        closure.accessibleIndex = obs.accessibleIndex;
    else
        closure.accessibleIndex = [];
    end
end


function [Vdirected, baselineIndex, orientationSign, found] = ...
    get_directed_visibility( ...
        visibilityMatrix, baselinePairs, collectorA, collectorB)
%GET_DIRECTED_VISIBILITY
%
% Return the visibility associated with the directed baseline
%
%   B_AB = r_B - r_A.
%
% If baselinePairs stores [A B], the visibility is used directly.
% If baselinePairs stores [B A], its complex conjugate is used because
%
%   V_BA = V(-u_AB,-v_AB) = conj(V_AB)
%
% for a real sky brightness distribution.

    directIndex = find( ...
        baselinePairs(:,1) == collectorA & ...
        baselinePairs(:,2) == collectorB);

    reverseIndex = find( ...
        baselinePairs(:,1) == collectorB & ...
        baselinePairs(:,2) == collectorA);

    if numel(directIndex) > 1 || numel(reverseIndex) > 1
        error('Duplicate directed baseline definitions were found.');
    end

    if ~isempty(directIndex)

        baselineIndex = directIndex;
        orientationSign = 1;
        Vdirected = visibilityMatrix(:,directIndex);
        found = true;

    elseif ~isempty(reverseIndex)

        baselineIndex = reverseIndex;
        orientationSign = -1;
        Vdirected = conj(visibilityMatrix(:,reverseIndex));
        found = true;

    else

        baselineIndex = 0;
        orientationSign = 0;
        Vdirected = complex(zeros(size(visibilityMatrix,1),1));
        found = false;
    end
end