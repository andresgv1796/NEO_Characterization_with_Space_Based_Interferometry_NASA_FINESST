function cen = gmos_center_eigpair_from_eigs(V, D, mode, opts)
%GMOS_CENTER_EIGPAIR_FROM_EIGS Select a non-trivial GMOS center eigenpair.
%
%   cen = gmos_center_eigpair_from_eigs(V,D,mode)
%   cen = gmos_center_eigpair_from_eigs(V,D,mode,opts)
%
% Inputs
%   V    : 6x6 right-eigenvector matrix
%   D    : 6x6 diagonal Floquet-multiplier matrix
%   mode : 'planar' (default) or 'vertical'
%   opts : optional struct
%          .tolUnit      default 1e-6
%          .tolImag      default 1e-10
%          .warnMismatch default true
%
% This is the CANONICAL GMOS center-mode selector.
%
% The generic spectral bookkeeping, including removal of the two trivial
% +1 multipliers, is delegated to analyze_floquet_spectrum.
%
% Mode character
% --------------
% For a normalized eigenvector w,
%
%   S_p = ||[wx wy wvx wvy]||
%   S_v = ||[wz wvz]||
%
% and
%
%   f_v = S_v^2/(S_p^2 + S_v^2).
%
% Thus:
%
%   f_v > 0.5  -> vertical-dominant
%   f_v < 0.5  -> planar-dominant
%   f_v = 0.5  -> exactly mixed
%
% This replaces the previous arbitrary "score > 0.2" interpretation.
% A requested mode is NOT rejected when it is not dominant; a warning is
% issued because mixed center modes may still be useful.
%
% Output
%   cen.lambda
%   cen.w
%   cen.rho0
%   cen.rho0_cycles
%   cen.idx
%   cen.proj.planar
%   cen.proj.vertical
%   cen.proj.verticalFraction
%   cen.proj.planarFraction
%   cen.proj.verticalDominant
%   cen.proj.planarDominant
%   cen.proj.requestedModeDominant
%   cen.proj.character
%   cen.trivialIdx
%   cen.trivialLambda
%   cen.candidateIdx
%   cen.candidateLambda
%   cen.allLambda

    if nargin < 3 || isempty(mode)
        mode = 'planar';
    end

    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts,'tolUnit')
        opts.tolUnit = 1e-6;
    end

    if ~isfield(opts,'tolImag')
        opts.tolImag = 1e-10;
    end

    if ~isfield(opts,'warnMismatch')
        opts.warnMismatch = true;
    end

    if ~(strcmpi(mode,'planar') || strcmpi(mode,'vertical'))
        error('mode must be ''planar'' or ''vertical''.');
    end


    %% Generic spectrum analysis

    spec = analyze_floquet_spectrum(V,D,opts);

    cand = spec.positiveUnitComplexIdx;

    if isempty(cand)

        error( ...
            ['No nontrivial complex unit-magnitude eigenvalue candidates ', ...
             'found after removing the two trivial +1 multipliers.']);

    end


    %% Select the candidate that best matches the requested character

    if strcmpi(mode,'vertical')

        score = spec.verticalScore(cand);

    else

        score = spec.planarScore(cand);

    end

    [~,kbest] = max(score);

    i0 = cand(kbest);

    lambda = spec.lambda(i0);
    w = V(:,i0);


    %% Rotation angle

    % lambda = exp(i*rho0)
    %
    % The positive-imaginary member is retained, so rho0 is normally in
    % (0,pi).  A later QPO-initialization convention check may flip its sign.

    rho0 = angle(lambda);


    %% Physically interpretable mode character

    Sp = spec.planarScore(i0);
    Sv = spec.verticalScore(i0);

    denom = Sp^2 + Sv^2;

    if denom > 0

        fv = Sv^2/denom;
        fp = Sp^2/denom;

    else

        fv = NaN;
        fp = NaN;

    end

    verticalDominant = isfinite(fv) && fv > 0.5;
    planarDominant   = isfinite(fp) && fp > 0.5;

    if verticalDominant

        character = "vertical-dominant";

    elseif planarDominant

        character = "planar-dominant";

    else

        character = "mixed";

    end

    if strcmpi(mode,'vertical')

        requestedModeDominant = verticalDominant;

    else

        requestedModeDominant = planarDominant;

    end


    %% Package output

    cen = struct();

    cen.lambda = lambda;
    cen.w = w;

    cen.rho0 = rho0;
    cen.rho0_cycles = rho0/(2*pi);

    cen.idx = i0;

    cen.proj = struct();

    cen.proj.planar = Sp;
    cen.proj.vertical = Sv;

    cen.proj.planarFraction = fp;
    cen.proj.verticalFraction = fv;

    cen.proj.positionVerticalFraction = ...
        spec.positionVerticalFraction(i0);

    cen.proj.velocityVerticalFraction = ...
        spec.velocityVerticalFraction(i0);

    cen.proj.planarDominant = planarDominant;
    cen.proj.verticalDominant = verticalDominant;

    cen.proj.requestedModeDominant = requestedModeDominant;
    cen.proj.character = character;

    cen.trivialIdx = spec.trivialIdx;
    cen.trivialLambda = spec.trivialLambda;

    cen.candidateIdx = cand(:);
    cen.candidateLambda = spec.lambda(cand);

    cen.allLambda = spec.lambda;

    cen.tolUnit = opts.tolUnit;
    cen.tolImag = opts.tolImag;


    %% Informative warning only

    if opts.warnMismatch && ~requestedModeDominant

        warning( ...
            ['Requested mode "%s" selected a %s center eigenvector ', ...
             '(planar fraction %.3f, vertical fraction %.3f). ', ...
             'The mode is retained, but the requested character is not dominant.'], ...
            mode, ...
            character, ...
            fp, ...
            fv);

    end

end
