function jList = qpo_equal_phase_indices(N, nCollectors, startIndex)
%QPO_EQUAL_PHASE_INDICES Select equally spaced samples on the GMOS curve.

    if nargin < 3 || isempty(startIndex)
        startIndex = 1;
    end

    if nCollectors < 2
        error('nCollectors must be at least 2.');
    end

    phaseFractions = (0:nCollectors-1) / nCollectors;
    jList = startIndex + round(phaseFractions * N);
    jList = 1 + mod(jList - 1, N);

    if numel(unique(jList)) ~= numel(jList)
        error('Repeated collector indices. Increase N or change nCollectors.');
    end
end
