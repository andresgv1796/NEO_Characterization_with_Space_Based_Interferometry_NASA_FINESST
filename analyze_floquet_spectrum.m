function spec = analyze_floquet_spectrum(V, D, opts)
%ANALYZE_FLOQUET_SPECTRUM Analyze a 6x6 CR3BP Floquet eigensystem.
%
%   spec = analyze_floquet_spectrum(V,D)
%   spec = analyze_floquet_spectrum(V,D,opts)
%
% Inputs
%   V    : 6x6 matrix of right eigenvectors
%   D    : 6x6 diagonal matrix of Floquet multipliers
%   opts : optional struct
%          .tolUnit  (default 1e-6)
%          .tolImag  (default 1e-10)
%
% This function contains the GENERIC Floquet-spectrum bookkeeping used by
% both the production GMOS center-mode selector and the whole-catalog
% survey.  It does not decide which center mode GMOS should use.
%
% In particular it:
%   - identifies the two trivial multipliers nearest +1;
%   - defines the remaining four as the non-trivial spectrum;
%   - computes spectral radius and the JPL-style stability index;
%   - checks determinant and reciprocal-pair structure;
%   - computes planar / vertical eigenvector projection scores;
%   - identifies genuine non-trivial complex unit-circle candidates;
%   - classifies the four non-trivial multipliers.
%
% Output fields include
%   lambda, absLambda
%   trivialIdx, trivialLambda, trivialPairDistance
%   nontrivialIdx, nontrivialLambda
%   spectralRadius, stabilityIndex
%   determinantError, reciprocalPairError
%   planarScore, verticalScore, verticalFraction
%   positiveUnitComplexIdx, nPositiveUnitComplexPairs
%   spectralClass, spectralClassCode
%   dominantUnstableIdx, dominantUnstableLambda
%   dominantStableIdx, dominantStableLambda

    if nargin < 3 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts,'tolUnit')
        opts.tolUnit = 1e-6;
    end

    if ~isfield(opts,'tolImag')
        opts.tolImag = 1e-10;
    end

    assert(isequal(size(V),[6 6]),'V must be 6x6.');
    assert(isequal(size(D),[6 6]),'D must be 6x6.');

    lam = diag(D);

    spec = struct();

    spec.lambda = lam;
    spec.absLambda = abs(lam);

    spec.tolUnit = opts.tolUnit;
    spec.tolImag = opts.tolImag;


    %% --------------------------------------------------------
    %  Trivial autonomous multipliers
    %  --------------------------------------------------------

    [~,idxNearOne] = sort(abs(lam - 1),'ascend');

    trivialIdx = idxNearOne(1:2);

    nontrivialMask = true(6,1);
    nontrivialMask(trivialIdx) = false;

    nontrivialIdx = find(nontrivialMask);

    spec.trivialIdx = trivialIdx(:);
    spec.trivialLambda = lam(trivialIdx);
    spec.trivialPairDistance = max(abs(lam(trivialIdx) - 1));

    spec.nontrivialIdx = nontrivialIdx(:);
    spec.nontrivialLambda = lam(nontrivialIdx);


    %% --------------------------------------------------------
    %  Global symplectic diagnostics
    %  --------------------------------------------------------

    spec.spectralRadius = max(abs(lam));

    spec.stabilityIndex = 0.5*( ...
        spec.spectralRadius + ...
        1/spec.spectralRadius);

    spec.determinantError = abs(prod(lam) - 1);

    spec.reciprocalPairError = local_reciprocal_pair_error(lam);


    %% --------------------------------------------------------
    %  Eigenvector mode-character projections
    %  --------------------------------------------------------
    %
    % Full state-space:
    %   S_p = ||[wx wy wvx wvy]||
    %   S_v = ||[wz wvz]||
    %   f_v,state = S_v^2/(S_p^2 + S_v^2)
    %
    % Position only:
    %   f_v,pos = |wz|^2/(|wx|^2 + |wy|^2 + |wz|^2)
    %
    % Velocity only:
    %   f_v,vel = |wvz|^2/(|wvx|^2 + |wvy|^2 + |wvz|^2)
    %
    % The full state-space fraction remains the production classification
    % quantity. The position-only fraction is retained because formation
    % geometry depends directly on the position part of the Floquet mode.

    planarScore = zeros(6,1);
    verticalScore = zeros(6,1);

    verticalFraction = nan(6,1);
    positionVerticalFraction = nan(6,1);
    velocityVerticalFraction = nan(6,1);

    for j = 1:6

        w = V(:,j);

        planarScore(j) = norm([w(1);w(2);w(4);w(5)]);
        verticalScore(j) = norm([w(3);w(6)]);

        denomState = planarScore(j)^2 + verticalScore(j)^2;
        if denomState > 0
            verticalFraction(j) = verticalScore(j)^2/denomState;
        end

        denomPos = abs(w(1))^2 + abs(w(2))^2 + abs(w(3))^2;
        if denomPos > 0
            positionVerticalFraction(j) = abs(w(3))^2/denomPos;
        end

        denomVel = abs(w(4))^2 + abs(w(5))^2 + abs(w(6))^2;
        if denomVel > 0
            velocityVerticalFraction(j) = abs(w(6))^2/denomVel;
        end

    end

    spec.planarScore = planarScore;
    spec.verticalScore = verticalScore;

    spec.verticalFraction = verticalFraction;
    spec.positionVerticalFraction = positionVerticalFraction;
    spec.velocityVerticalFraction = velocityVerticalFraction;


    %% --------------------------------------------------------
    %  Genuine non-trivial unit-circle center candidates
    %  --------------------------------------------------------

    centerMask = ...
        nontrivialMask & ...
        abs(abs(lam)-1) < opts.tolUnit & ...
        abs(imag(lam)) > opts.tolImag & ...
        imag(lam) > 0;

    spec.positiveUnitComplexIdx = find(centerMask);
    spec.nPositiveUnitComplexPairs = numel(spec.positiveUnitComplexIdx);


    %% --------------------------------------------------------
    %  Dominant stable / unstable representatives
    %  --------------------------------------------------------
    %
    % These are simply the non-trivial multipliers with smallest/largest
    % modulus.  In a saddle-center case they are the familiar real
    % reciprocal stable/unstable pair.  In a complex quartet they are
    % representatives of the inner/outer reciprocal pairs instead.

    nontrivialAbs = abs(spec.nontrivialLambda);

    [~,kMax] = max(nontrivialAbs);
    [~,kMin] = min(nontrivialAbs);

    spec.dominantUnstableIdx = nontrivialIdx(kMax);
    spec.dominantUnstableLambda = lam(spec.dominantUnstableIdx);

    spec.dominantStableIdx = nontrivialIdx(kMin);
    spec.dominantStableLambda = lam(spec.dominantStableIdx);


    %% --------------------------------------------------------
    %  Four-multiplier topology
    %  --------------------------------------------------------

    [spec.spectralClass,spec.spectralClassCode] = ...
        local_classify_nontrivial_spectrum( ...
            spec.nontrivialLambda, ...
            opts.tolUnit, ...
            opts.tolImag);

end


function err = local_reciprocal_pair_error(lam)
%LOCAL_RECIPROCAL_PAIR_ERROR

    n = numel(lam);
    bestResidual = nan(n,1);

    for i = 1:n

        residuals = abs(lam(i)*lam - 1);

        bestResidual(i) = min(residuals);

    end

    err = max(bestResidual);

end


function [name,code] = local_classify_nontrivial_spectrum( ...
    lam, ...
    tolUnit, ...
    tolImag)
%LOCAL_CLASSIFY_NONTRIVIAL_SPECTRUM
%
% Codes:
%   1 = saddle-center
%   2 = elliptic-elliptic
%   3 = saddle-saddle
%   4 = complex quartet
%   5 = transition/degenerate

    lam = lam(:);

    isUnit = abs(abs(lam)-1) < tolUnit;
    isComplex = abs(imag(lam)) > tolImag;

    nUnitPositiveComplex = nnz( ...
        isUnit & ...
        isComplex & ...
        imag(lam) > 0);

    isOff = ~isUnit;
    nOff = nnz(isOff);

    offComplex = isOff & isComplex;
    offReal = isOff & ~isComplex;

    if nUnitPositiveComplex == 1 && ...
       nOff == 2 && ...
       nnz(offReal) == 2

        name = "saddle-center";
        code = 1;

    elseif nUnitPositiveComplex == 2 && ...
           nOff == 0

        name = "elliptic-elliptic";
        code = 2;

    elseif nUnitPositiveComplex == 0 && ...
           nOff == 4 && ...
           nnz(offReal) == 4

        name = "saddle-saddle";
        code = 3;

    elseif nUnitPositiveComplex == 0 && ...
           nOff == 4 && ...
           nnz(offComplex) == 4

        name = "complex-quartet";
        code = 4;

    else

        name = "transition/degenerate";
        code = 5;

    end

end
