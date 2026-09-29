function clean = make_clean_image_from_dirty(dirty, cleanParams)
%MAKE_CLEAN_IMAGE_FROM_DIRTY
%
% Historically faithful discrete image-plane Högbom CLEAN.
%
% The algorithm:
%
%   1. Construct the dirty beam on the complete difference-coordinate grid.
%
%   2. Initialize
%
%          M^(0) = 0
%          R^(0) = dirty image.
%
%   3. At every iteration, find the maximum absolute residual over the
%      complete image.
%
%   4. Preserve the sign of the selected residual value.
%
%   5. Add loopGain times that value to the delta-component model.
%
%   6. Subtract the exact shifted dirty beam over the complete finite image.
%
%   7. Stop at stopThreshold or nIterMax.
%
%   8. Fit an elliptical Gaussian clean beam to the connected central lobe
%      of the dirty beam.
%
%   9. Restore the delta components using the clean beam and add the final
%      residual.
%
% Required dirty fields:
%
%   dirty.image
%   dirty.beam
%   dirty.l_mas
%   dirty.m_mas
%   dirty.l_rad
%   dirty.m_rad
%   dirty.L_mas
%   dirty.M_mas
%   dirty.L_rad
%   dirty.M_rad
%   dirty.u
%   dirty.v
%   dirty.weights
%
% Required cleanParams fields:
%
%   cleanParams.loopGain
%   cleanParams.nIterMax
%   cleanParams.stopThreshold
%   cleanParams.cleanBeamFitLevel
%   cleanParams.cleanBeamTruncateSigma
%
% The function intentionally has:
%
%   no CLEAN mask,
%   no positive-only mode,
%   no zero-padded dirty-beam shifting,
%   no major cycles,
%   no residual scaling,
%   no multiscale components.
%
% This is the direct discrete Högbom base case.

    %% Validate required structure fields

    require_structure_field(dirty, 'image', 'dirty');
    require_structure_field(dirty, 'beam', 'dirty');

    require_structure_field(dirty, 'l_mas', 'dirty');
    require_structure_field(dirty, 'm_mas', 'dirty');
    require_structure_field(dirty, 'l_rad', 'dirty');
    require_structure_field(dirty, 'm_rad', 'dirty');

    require_structure_field(dirty, 'L_mas', 'dirty');
    require_structure_field(dirty, 'M_mas', 'dirty');
    require_structure_field(dirty, 'L_rad', 'dirty');
    require_structure_field(dirty, 'M_rad', 'dirty');

    require_structure_field(dirty, 'u', 'dirty');
    require_structure_field(dirty, 'v', 'dirty');
    require_structure_field(dirty, 'weights', 'dirty');

    require_structure_field(cleanParams, 'loopGain', 'cleanParams');
    require_structure_field(cleanParams, 'nIterMax', 'cleanParams');
    require_structure_field(cleanParams, 'stopThreshold', 'cleanParams');

    require_structure_field( ...
        cleanParams, ...
        'cleanBeamFitLevel', ...
        'cleanParams');

    require_structure_field( ...
        cleanParams, ...
        'cleanBeamTruncateSigma', ...
        'cleanParams');


    %% Read image and measurement arrays

    dirtyImageComplex = dirty.image;
    dirtyBeamImageComplex = dirty.beam;

    u = dirty.u(:);
    v = dirty.v(:);
    visibilityWeights = dirty.weights(:);

    l_rad = dirty.l_rad(:).';
    m_rad = dirty.m_rad(:);

    l_mas = dirty.l_mas(:).';
    m_mas = dirty.m_mas(:);


    %% Validate dimensions

    [nM, nL] = size(dirtyImageComplex);

    if ~isequal(size(dirtyBeamImageComplex), [nM, nL])
        error('dirty.beam must have the same dimensions as dirty.image.');
    end

    if numel(l_rad) ~= nL || numel(l_mas) ~= nL
        error('The l-coordinate vectors must match the image columns.');
    end

    if numel(m_rad) ~= nM || numel(m_mas) ~= nM
        error('The m-coordinate vectors must match the image rows.');
    end

    if numel(u) ~= numel(v) || numel(u) ~= numel(visibilityWeights)
        error('dirty.u, dirty.v, and dirty.weights must have equal length.');
    end

    if any(~isfinite(u)) || any(~isfinite(v))
        error('The UV coordinates contain nonfinite values.');
    end

    if any(~isfinite(visibilityWeights))
        error('The visibility weights contain nonfinite values.');
    end

    if any(visibilityWeights < 0)
        error('Historical Högbom CLEAN requires nonnegative weights.');
    end

    sumWeights = sum(visibilityWeights);

    if sumWeights <= 0
        error('The visibility-weight sum must be positive.');
    end


    %% Validate image grids

    validate_uniform_coordinate_grid(l_rad, 'dirty.l_rad');
    validate_uniform_coordinate_grid(m_rad, 'dirty.m_rad');

    validate_uniform_coordinate_grid(l_mas, 'dirty.l_mas');
    validate_uniform_coordinate_grid(m_mas, 'dirty.m_mas');


    %% Validate CLEAN parameters

    if ~isscalar(cleanParams.loopGain) || ...
            ~isfinite(cleanParams.loopGain) || ...
            cleanParams.loopGain <= 0 || ...
            cleanParams.loopGain > 1

        error('cleanParams.loopGain must satisfy 0 < loopGain <= 1.');
    end

    if ~isscalar(cleanParams.nIterMax) || ...
            cleanParams.nIterMax < 1 || ...
            cleanParams.nIterMax ~= floor(cleanParams.nIterMax)

        error('cleanParams.nIterMax must be a positive integer.');
    end

    if ~isscalar(cleanParams.stopThreshold) || ...
            ~isfinite(cleanParams.stopThreshold) || ...
            cleanParams.stopThreshold < 0

        error('cleanParams.stopThreshold must be finite and nonnegative.');
    end

    if ~isscalar(cleanParams.cleanBeamFitLevel) || ...
            cleanParams.cleanBeamFitLevel <= 0 || ...
            cleanParams.cleanBeamFitLevel >= 1

        error([ ...
            'cleanParams.cleanBeamFitLevel must lie strictly ', ...
            'between zero and one.']);
    end

    if ~isscalar(cleanParams.cleanBeamTruncateSigma) || ...
            ~isfinite(cleanParams.cleanBeamTruncateSigma) || ...
            cleanParams.cleanBeamTruncateSigma <= 0

        error([ ...
            'cleanParams.cleanBeamTruncateSigma must be ', ...
            'finite and positive.']);
    end


    %% The historical dirty map must be real

    validate_numerically_real( ...
        dirtyImageComplex, ...
        'dirty.image');

    validate_numerically_real( ...
        dirtyBeamImageComplex, ...
        'dirty.beam');

    dirtyImageInput = real(dirtyImageComplex);
    dirtyBeamImageInput = real(dirtyBeamImageComplex);


    %% Construct the exact dirty beam on the complete difference grid
    %
    % For an nM-by-nL image, every source-image pixel pair can differ by
    %
    %   Delta i = -(nL-1), ..., +(nL-1)
    %   Delta j = -(nM-1), ..., +(nM-1).
    %
    % Therefore, the required dirty-beam support has dimensions
    %
    %   (2*nM-1)-by-(2*nL-1).

    [dirtyBeamDifference, deltaL_rad, deltaM_rad] = ...
        make_exact_dirty_beam_difference_grid( ...
            u, ...
            v, ...
            visibilityWeights, ...
            l_rad, ...
            m_rad);

    validate_numerically_real( ...
        dirtyBeamDifference, ...
        'full difference-grid dirty beam');

    dirtyBeamDifference = real(dirtyBeamDifference);


    %% Normalize at the zero-lag dirty-beam center

    jDifferenceZero = nM;
    iDifferenceZero = nL;

    dirtyBeamCenter = ...
        dirtyBeamDifference(jDifferenceZero, iDifferenceZero);

    if dirtyBeamCenter <= 0
        error('The zero-lag dirty-beam value must be positive.');
    end

    dirtyBeamDifference = ...
        dirtyBeamDifference/dirtyBeamCenter;

    dirtyImage = ...
        dirtyImageInput/dirtyBeamCenter;

    dirtyBeamImage = ...
        dirtyBeamImageInput/dirtyBeamCenter;

    stopThreshold = ...
        cleanParams.stopThreshold/dirtyBeamCenter;


    %% Fit the clean beam to the exact central dirty-beam lobe

    rad2mas = (180/pi)*3600*1000;

    deltaL_mas = deltaL_rad*rad2mas;
    deltaM_mas = deltaM_rad*rad2mas;

    cleanBeamFit = fit_clean_beam_to_dirty_main_lobe( ...
        dirtyBeamDifference, ...
        deltaL_mas, ...
        deltaM_mas, ...
        jDifferenceZero, ...
        iDifferenceZero, ...
        cleanParams.cleanBeamFitLevel);

    [cleanBeamKernel, cleanBeamImage] = ...
        make_fitted_clean_beam_kernel( ...
            l_mas, ...
            m_mas, ...
            cleanBeamFit, ...
            cleanParams.cleanBeamTruncateSigma);


    %% Initialize Högbom model and residual

    model = zeros(nM, nL);
    residual = dirtyImage;


    %% Allocate iteration history

    nIterMax = cleanParams.nIterMax;

    histPeak = zeros(nIterMax, 1);
    histAbsPeak = zeros(nIterMax, 1);
    histComponentAmplitude = zeros(nIterMax, 1);

    histLmas = zeros(nIterMax, 1);
    histMmas = zeros(nIterMax, 1);

    histResidualRms = zeros(nIterMax, 1);
    histResidualMaxAbs = zeros(nIterMax, 1);
    histModelFlux = zeros(nIterMax, 1);

    nDone = 0;
    stoppingReason = 'maximum iterations reached';


    %% Exact discrete Högbom subtraction loop

    for iter = 1:nIterMax

        % Search the complete residual map for its largest absolute value.
        [peakAbs, idxPeak] = max(abs(residual(:)));

        if peakAbs <= stopThreshold
            stoppingReason = 'threshold reached';
            break
        end

        peakValue = residual(idxPeak);

        [jPeak, iPeak] = ind2sub( ...
            size(residual), ...
            idxPeak);

        componentAmplitude = ...
            cleanParams.loopGain*peakValue;

        % Add the signed delta component to the CLEAN model.
        model(jPeak, iPeak) = ...
            model(jPeak, iPeak) + componentAmplitude;

        % Extract the exact dirty-beam pattern
        %
        %   B_D(l_i-l_iPeak, m_j-m_jPeak)
        %
        % for every image pixel. No zero padding and no circular wrapping
        % are used.

        shiftedDirtyBeam = ...
            extract_exact_shifted_dirty_beam( ...
                dirtyBeamDifference, ...
                nM, ...
                nL, ...
                jPeak, ...
                iPeak);

        % Historical Högbom residual update.
        residual = ...
            residual - componentAmplitude*shiftedDirtyBeam;

        nDone = iter;

        histPeak(iter) = peakValue;
        histAbsPeak(iter) = peakAbs;
        histComponentAmplitude(iter) = componentAmplitude;

        histLmas(iter) = l_mas(iPeak);
        histMmas(iter) = m_mas(jPeak);

        histResidualRms(iter) = ...
            sqrt(mean(residual(:).^2));

        histResidualMaxAbs(iter) = ...
            max(abs(residual(:)));

        histModelFlux(iter) = ...
            sum(model(:));
    end


    %% Historical restoration
    %
    % Replace every delta component by the fitted Gaussian clean beam.

    restoredComponents = conv2( ...
        model, ...
        cleanBeamKernel, ...
        'same');

    % Högbom's final map:
    %
    %   restored CLEAN image = restored components + final residual.

    cleanImage = restoredComponents + residual;


    %% Record the parameters actually used

    paramsUsed = cleanParams;

    paramsUsed.stopThreshold = stopThreshold;

    paramsUsed.cleanBeamFwhmMajor_mas = ...
        cleanBeamFit.fwhmMajor_mas;

    paramsUsed.cleanBeamFwhmMinor_mas = ...
        cleanBeamFit.fwhmMinor_mas;

    paramsUsed.cleanBeamAxisRatio = ...
        cleanBeamFit.axisRatio;

    paramsUsed.cleanBeamPaDeg = ...
        cleanBeamFit.paDeg;


    %% Output structure

    clean = struct();

    clean.l_mas = dirty.l_mas;
    clean.m_mas = dirty.m_mas;
    clean.L_mas = dirty.L_mas;
    clean.M_mas = dirty.M_mas;

    clean.l_rad = dirty.l_rad;
    clean.m_rad = dirty.m_rad;
    clean.L_rad = dirty.L_rad;
    clean.M_rad = dirty.M_rad;

    clean.params = paramsUsed;

    clean.dirtyImageNormalized = dirtyImage;
    clean.dirtyBeamNormalized = dirtyBeamImage;

    clean.dirtyBeamCenter = dirtyBeamCenter;

    clean.model = model;
    clean.residual = residual;

    clean.restoredComponents = restoredComponents;
    clean.image = cleanImage;

    clean.cleanBeamKernel = cleanBeamKernel;
    clean.cleanBeam = cleanBeamImage;
    clean.cleanBeamFit = cleanBeamFit;

    clean.nIter = nDone;
    clean.stoppingReason = stoppingReason;

    clean.history.iteration = (1:nDone).';

    clean.history.peak = histPeak(1:nDone);
    clean.history.absPeak = histAbsPeak(1:nDone);

    clean.history.componentAmplitude = ...
        histComponentAmplitude(1:nDone);

    clean.history.l_mas = histLmas(1:nDone);
    clean.history.m_mas = histMmas(1:nDone);

    clean.history.residualRms = ...
        histResidualRms(1:nDone);

    clean.history.residualMaxAbs = ...
        histResidualMaxAbs(1:nDone);

    clean.history.modelFlux = ...
        histModelFlux(1:nDone);
end


function [dirtyBeamDifference, deltaL_rad, deltaM_rad] = ...
    make_exact_dirty_beam_difference_grid( ...
        u, v, visibilityWeights, l_rad, m_rad)
%MAKE_EXACT_DIRTY_BEAM_DIFFERENCE_GRID
%
% Evaluate the weighted dirty beam at every possible image-pixel
% separation:
%
%   B_D(Delta l, Delta m)
%       = sum_k w_k exp[2*pi*i*(u_k*Delta l + v_k*Delta m)]
%         -----------------------------------------------------
%                           sum_k w_k.

    nL = numel(l_rad);
    nM = numel(m_rad);

    deltaLGrid = l_rad(2) - l_rad(1);
    deltaMGrid = m_rad(2) - m_rad(1);

    deltaL_rad = ...
        (-(nL-1):(nL-1))*deltaLGrid;

    deltaM_rad = ...
        (-(nM-1):(nM-1))*deltaMGrid;

    nDeltaL = numel(deltaL_rad);
    nDeltaM = numel(deltaM_rad);

    sumWeights = sum(visibilityWeights);

    % The l-dependent exponential is reused for every Delta m row.
    exponentialL = exp( ...
        1i*2*pi*(u*deltaL_rad));

    dirtyBeamDifference = ...
        complex(zeros(nDeltaM, nDeltaL));

    parfor jDifference = 1:nDeltaM

        exponentialM = exp( ...
            1i*2*pi*v*deltaM_rad(jDifference));

        weightedM = ...
            visibilityWeights.*exponentialM;

        dirtyBeamDifference(jDifference,:) = ...
            (weightedM.'*exponentialL)/sumWeights;
    end
end


function shiftedDirtyBeam = extract_exact_shifted_dirty_beam( ...
    dirtyBeamDifference, nM, nL, jPeak, iPeak)
%EXTRACT_EXACT_SHIFTED_DIRTY_BEAM
%
% Return
%
%   B_D(l_i-l_iPeak, m_j-m_jPeak)
%
% for every image pixel.
%
% The zero difference is located at
%
%   dirtyBeamDifference(nM,nL).

    rowIndices = ...
        (1:nM) - jPeak + nM;

    columnIndices = ...
        (1:nL) - iPeak + nL;

    shiftedDirtyBeam = ...
        dirtyBeamDifference(rowIndices, columnIndices);
end


function fit = fit_clean_beam_to_dirty_main_lobe( ...
    dirtyBeam, deltaL_mas, deltaM_mas, ...
    jBeam0, iBeam0, fitLevel)
%FIT_CLEAN_BEAM_TO_DIRTY_MAIN_LOBE
%
% Fit an elliptical Gaussian to the connected central dirty-beam lobe:
%
%   B(x,y) = exp[-0.5*[x y]*A*[x;y]].
%
% Hence
%
%   -2*log(B)
%       = A11*x^2 + 2*A12*x*y + A22*y^2.

    centralMask = dirtyBeam >= fitLevel;

    if ~centralMask(jBeam0, iBeam0)
        error('The zero-lag dirty-beam point is not in the fit mask.');
    end

    connectedMask = extract_connected_component( ...
        centralMask, ...
        jBeam0, ...
        iBeam0);

    [jFit, iFit] = find(connectedMask);

    lOffset = ...
        deltaL_mas(iFit).' - deltaL_mas(iBeam0);

    mOffset = ...
        deltaM_mas(jFit).' - deltaM_mas(jBeam0);

    linearIndex = ...
        sub2ind(size(dirtyBeam), jFit, iFit);

    beamValues = dirtyBeam(linearIndex);

    valid = ...
        beamValues > 0 & ...
        beamValues < 1 & ...
        isfinite(beamValues);

    lOffset = lOffset(valid);
    mOffset = mOffset(valid);
    beamValues = beamValues(valid);

    if numel(beamValues) < 6
        error([ ...
            'Too few pixels were available for the clean-beam fit. ', ...
            'Reduce cleanBeamFitLevel or refine the image grid.']);
    end

    designMatrix = [
        lOffset.^2, ...
        2*lOffset.*mOffset, ...
        mOffset.^2
    ];

    rightHandSide = ...
        -2*log(beamValues);

    if rank(designMatrix) < 3
        error('The clean-beam fit matrix is rank deficient.');
    end

    coefficients = ...
        designMatrix\rightHandSide;

    precisionMatrix = [
        coefficients(1), coefficients(2)
        coefficients(2), coefficients(3)
    ];

    precisionMatrix = ...
        0.5*(precisionMatrix + precisionMatrix.');

    [eigenvectors, eigenvalueMatrix] = ...
        eig(precisionMatrix);

    eigenvalues = diag(eigenvalueMatrix);

    if any(~isfinite(eigenvalues)) || any(eigenvalues <= 0)

    error([ ...
        'The fitted Gaussian precision matrix is not positive definite.\n', ...
        '  cleanBeamFitLevel = %.3f\n', ...
        '  fitted pixels     = %d\n', ...
        '  eigenvalues       = [%+.6e %+.6e]\n', ...
        'The selected dirty-beam region is not adequately described by ', ...
        'an elliptical Gaussian. Increase cleanBeamFitLevel so that the ', ...
        'fit is restricted more closely to the central beam peak.'], ...
        fitLevel, ...
        numel(beamValues), ...
        eigenvalues(1), ...
        eigenvalues(2));

    end

    [eigenvalues, order] = ...
        sort(eigenvalues, 'ascend');

    eigenvectors = eigenvectors(:,order);

    sigmaMajor_mas = ...
        1/sqrt(eigenvalues(1));

    sigmaMinor_mas = ...
        1/sqrt(eigenvalues(2));

    fwhmFactor = ...
        2*sqrt(2*log(2));

    fwhmMajor_mas = ...
        fwhmFactor*sigmaMajor_mas;

    fwhmMinor_mas = ...
        fwhmFactor*sigmaMinor_mas;

    majorAxisVector = ...
        eigenvectors(:,1);

    paDeg = atan2d( ...
        majorAxisVector(2), ...
        majorAxisVector(1));

    paDeg = ...
        mod(paDeg + 90, 180) - 90;

    fittedExponent = ...
        coefficients(1)*lOffset.^2 + ...
        2*coefficients(2)*lOffset.*mOffset + ...
        coefficients(3)*mOffset.^2;

    fittedValues = ...
        exp(-0.5*fittedExponent);

    rmsError = sqrt(mean( ...
        (fittedValues - beamValues).^2));

    fit = struct();

    fit.fitLevel = fitLevel;
    fit.nFitPixels = numel(beamValues);

    fit.precisionMatrix_mas2 = precisionMatrix;

    fit.sigmaMajor_mas = sigmaMajor_mas;
    fit.sigmaMinor_mas = sigmaMinor_mas;

    fit.fwhmMajor_mas = fwhmMajor_mas;
    fit.fwhmMinor_mas = fwhmMinor_mas;

    fit.axisRatio = ...
        fwhmMinor_mas/fwhmMajor_mas;

    fit.paDeg = paDeg;
    fit.rmsError = rmsError;

    fit.connectedMask = connectedMask;
end


function connectedMask = extract_connected_component( ...
    binaryMask, jSeed, iSeed)
%EXTRACT_CONNECTED_COMPONENT
%
% Eight-connected flood fill without Image Processing Toolbox functions.

    [nRows, nCols] = size(binaryMask);

    connectedMask = false(nRows, nCols);

    nQueueMax = nnz(binaryMask);

    queueJ = zeros(nQueueMax, 1);
    queueI = zeros(nQueueMax, 1);

    queueHead = 1;
    queueTail = 1;

    queueJ(1) = jSeed;
    queueI(1) = iSeed;

    connectedMask(jSeed, iSeed) = true;

    while queueHead <= queueTail

        jCurrent = queueJ(queueHead);
        iCurrent = queueI(queueHead);

        queueHead = queueHead + 1;

        for deltaJ = -1:1
            for deltaI = -1:1

                if deltaJ == 0 && deltaI == 0
                    continue
                end

                jNeighbor = jCurrent + deltaJ;
                iNeighbor = iCurrent + deltaI;

                insideImage = ...
                    jNeighbor >= 1 && ...
                    jNeighbor <= nRows && ...
                    iNeighbor >= 1 && ...
                    iNeighbor <= nCols;

                if ~insideImage
                    continue
                end

                shouldAdd = ...
                    binaryMask(jNeighbor, iNeighbor) && ...
                    ~connectedMask(jNeighbor, iNeighbor);

                if shouldAdd

                    queueTail = queueTail + 1;

                    queueJ(queueTail) = jNeighbor;
                    queueI(queueTail) = iNeighbor;

                    connectedMask(jNeighbor, iNeighbor) = true;
                end
            end
        end
    end
end


function [kernel, beamImage] = make_fitted_clean_beam_kernel( ...
    l_mas, m_mas, fit, truncateSigma)
%MAKE_FITTED_CLEAN_BEAM_KERNEL
%
% Construct the fitted elliptical Gaussian restoring beam.

    deltaL = abs(l_mas(2) - l_mas(1));
    deltaM = abs(m_mas(2) - m_mas(1));

    halfWidthL = ceil( ...
        truncateSigma*fit.sigmaMajor_mas/deltaL);

    halfWidthM = ceil( ...
        truncateSigma*fit.sigmaMajor_mas/deltaM);

    halfWidthL = max(halfWidthL, 1);
    halfWidthM = max(halfWidthM, 1);

    lLocal = ...
        (-halfWidthL:halfWidthL)*deltaL;

    mLocal = ...
        (-halfWidthM:halfWidthM)*deltaM;

    [Llocal, Mlocal] = ...
        meshgrid(lLocal, mLocal);

    phi = deg2rad(fit.paDeg);

    Xmajor = ...
        cos(phi)*Llocal + sin(phi)*Mlocal;

    Yminor = ...
        -sin(phi)*Llocal + cos(phi)*Mlocal;

    ellipticalRadiusSquared = ...
        (Xmajor/fit.sigmaMajor_mas).^2 + ...
        (Yminor/fit.sigmaMinor_mas).^2;

    kernel = ...
        exp(-0.5*ellipticalRadiusSquared);

    kernel( ...
        ellipticalRadiusSquared > truncateSigma^2) = 0;

    kernel = ...
        kernel/max(kernel(:));

    beamImage = ...
        zeros(numel(m_mas), numel(l_mas));

    [~, iImageCenter] = ...
        min(abs(l_mas));

    [~, jImageCenter] = ...
        min(abs(m_mas));

    beamImage = add_kernel_at_pixel( ...
        beamImage, ...
        kernel, ...
        jImageCenter, ...
        iImageCenter, ...
        1.0);
end


function imageOut = add_kernel_at_pixel( ...
    imageIn, kernel, jTarget, iTarget, amplitude)
%ADD_KERNEL_AT_PIXEL

    imageOut = imageIn;

    [nRows, nCols] = size(imageIn);
    [nKernelRows, nKernelCols] = size(kernel);

    jKernelCenter = ...
        floor((nKernelRows + 1)/2);

    iKernelCenter = ...
        floor((nKernelCols + 1)/2);

    destinationRows = ...
        jTarget + ((1:nKernelRows) - jKernelCenter);

    destinationCols = ...
        iTarget + ((1:nKernelCols) - iKernelCenter);

    validRows = ...
        destinationRows >= 1 & ...
        destinationRows <= nRows;

    validCols = ...
        destinationCols >= 1 & ...
        destinationCols <= nCols;

    imageOut( ...
        destinationRows(validRows), ...
        destinationCols(validCols)) = ...
        imageOut( ...
            destinationRows(validRows), ...
            destinationCols(validCols)) + ...
        amplitude*kernel(validRows, validCols);
end


function validate_uniform_coordinate_grid(coordinate, coordinateName)
%VALIDATE_UNIFORM_COORDINATE_GRID

    coordinate = coordinate(:);

    if numel(coordinate) < 2
        error('%s must contain at least two points.', coordinateName);
    end

    spacing = diff(coordinate);

    referenceSpacing = spacing(1);

    tolerance = ...
        100*eps(max(1, max(abs(coordinate))));

    if max(abs(spacing-referenceSpacing)) > tolerance
        error('%s must be uniformly spaced.', coordinateName);
    end
end


function validate_numerically_real(arrayIn, arrayName)
%VALIDATE_NUMERICALLY_REAL

    realScale = ...
        max(1, max(abs(real(arrayIn(:)))));

    imaginaryMaximum = ...
        max(abs(imag(arrayIn(:))));

    tolerance = ...
        1e-10*realScale;

    if imaginaryMaximum > tolerance
        error([ ...
            '%s is not numerically real. Historical image-plane ', ...
            'Högbom CLEAN requires Hermitian visibility sampling.'], ...
            arrayName);
    end
end


function require_structure_field(structureIn, fieldName, structureName)
%REQUIRE_STRUCTURE_FIELD

    if ~isstruct(structureIn)
        error('%s must be a structure.', structureName);
    end

    if ~isfield(structureIn, fieldName)
        error('%s.%s is required.', structureName, fieldName);
    end
end