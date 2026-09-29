function [branchMaster, meta] = gmos_extend_fixT_mu_branch_weighted_bidirectional( ...
    branchMaster, N, Ttarget, ode_opts, optsPALC, extOpts)
%GMOS_EXTEND_FIXT_MU_BRANCH_WEIGHTED_BIDIRECTIONAL
%  Extend an existing weighted mu branch and merge automatically.
%
% This wrapper takes an already ordered branch and continues it from one end
% using the weighted mu PALC solver. It then stitches the new segment back
% into the master branch automatically.
%
% Supported continuation sides
%   extOpts.continuationSide = 'tail'  : continue from the high index end.
%   extOpts.continuationSide = 'head'  : continue from the low index end.
%
% On the side being extended, tailSkip removes suspicious edge members before
% building the seed pair. This means those skipped members are not kept in
% the merged output.
%
% Required extOpts fields
%   tailSkip          nonnegative integer
%   deltaMuScale      nonzero scalar multiplier applied to the local delta mu
%   nSteps            nonnegative integer number of new continuation steps
%   saveFile          filename or ''
%   verbose           logical
%
% Optional extOpts fields
%   continuationSide  'tail' or 'head'   default 'tail'
%
% Inputs
%   branchIn          existing ordered branch struct with at least two members
%   N                 number of invariant curve nodes
%   Ttarget           fixed stroboscopic time
%   ode_opts          ODE options
%   optsPALC          PALC options passed to the continuation function
%   extOpts           extension options listed above
%
% Outputs
%   branchOut         merged branch in the same ordered format as branchIn
%   meta              metadata for the extension and merge

    narginchk(6, 6);

    validate_branch(branchMaster, N);
    validateattributes(N, {'numeric'}, {'scalar','integer','>=',3}, mfilename, 'N', 2);
    validateattributes(Ttarget, {'numeric'}, {'scalar','real','finite'}, mfilename, 'Ttarget', 3);

    if nargin < 4 || isempty(ode_opts)
        ode_opts = [];
    end
    if nargin < 5 || isempty(optsPALC)
        optsPALC = struct();
    end
    if nargin < 6 || isempty(extOpts)
        extOpts = struct();
    elseif ~isstruct(extOpts)
        error('extOpts must be a struct.');
    end

    if ~isfield(extOpts, 'tailSkip') || isempty(extOpts.tailSkip)
        extOpts.tailSkip = 1;
    end
    if ~isfield(extOpts, 'deltaMuScale') || isempty(extOpts.deltaMuScale)
        extOpts.deltaMuScale = 0.25;
    end
    if ~isfield(extOpts, 'nSteps') || isempty(extOpts.nSteps)
        extOpts.nSteps = 100;
    end
    if ~isfield(extOpts, 'saveFile') || isempty(extOpts.saveFile)
        extOpts.saveFile = '';
    end
    if ~isfield(extOpts, 'verbose') || isempty(extOpts.verbose)
        extOpts.verbose = true;
    end
    if ~isfield(extOpts, 'continuationSide') || isempty(extOpts.continuationSide)
        extOpts.continuationSide = 'tail';
    end

    validateattributes(extOpts.tailSkip,     {'numeric'}, {'scalar','integer','>=',0}, mfilename, 'extOpts.tailSkip');
    validateattributes(extOpts.deltaMuScale, {'numeric'}, {'scalar','real','finite','nonzero'}, mfilename, 'extOpts.deltaMuScale');
    validateattributes(extOpts.nSteps,       {'numeric'}, {'scalar','integer','>=',0}, mfilename, 'extOpts.nSteps');
    side = lower(string(extOpts.continuationSide));
    if side ~= "tail" && side ~= "head"
        error('extOpts.continuationSide must be ''tail'' or ''head''.');
    end

    nOld = numel(branchMaster.y);
    if extOpts.nSteps == 0
        meta = build_meta_empty(branchMaster, side, extOpts);
        maybe_save(branchMaster, meta, extOpts.saveFile);
        return
    end

    [k0, k1, skippedIdx, goodIdx] = select_seed_indices(branchMaster, extOpts.tailSkip, side);

    y0_sol = branchMaster.y{k0};
    y1_sol = branchMaster.y{k1};
    y0_sol = y0_sol(:);
    y1_sol = y1_sol(:);

    if numel(y0_sol) ~= 6*N + 3 || numel(y1_sol) ~= 6*N + 3
        error('Selected seed members do not have length 6*N + 3.');
    end

    scale = gmos_build_palc_scaling_mu(y0_sol, y1_sol, N, struct());
    dySeedScaled = scale.w .* (y1_sol - y0_sol);
    dySeedNorm = norm(dySeedScaled);
    if dySeedNorm <= 0
        error('The weighted seed secant has zero norm.');
    end
    tHat0 = dySeedScaled / dySeedNorm;

    deltaMuSeed = y1_sol(end) - y0_sol(end);
    if abs(deltaMuSeed) <= eps(max(abs([y0_sol(end), y1_sol(end), 1])))
        error('deltaMuSeed is zero. The selected seed pair must have distinct mu values.');
    end

    if abs(tHat0(end)) <= eps(max(abs(tHat0)))
        error('The weighted tangent has zero mu component. Cannot build ds from delta mu.');
    end

    deltaMuTarget = extOpts.deltaMuScale * deltaMuSeed;
    ds = (scale.w(end) * deltaMuTarget) / tHat0(end);
    if ~isfinite(ds) || ds == 0
        error('Computed ds is invalid.');
    end

    optsUsed = optsPALC;
    optsUsed.continuationSide = char(side);

    if extOpts.verbose
        fprintf('Extending mu branch from the %s side\n', char(side));
        fprintf('Using seed indices k0 = %d, k1 = %d\n', k0, k1);
        fprintf('mu0 = %.16g\n', y0_sol(end));
        fprintf('mu1 = %.16g\n', y1_sol(end));
        fprintf('deltaMuSeed   = %.6e\n', deltaMuSeed);
        fprintf('deltaMuTarget = %.6e\n', deltaMuTarget);
        fprintf('ds            = %.6e\n', ds);
    end

    continueFcn = pick_continue_function();
    branchSeg = continueFcn(y0_sol, y1_sol, extOpts.nSteps, ds, N, Ttarget, ode_opts, optsUsed);

    branchMaster = merge_branch_segment(branchMaster, branchSeg, k0, k1, side);

    meta = struct();
    meta.side = char(side);
    meta.seedIdx = [k0, k1];
    meta.skippedIdx = skippedIdx;
    meta.deltaMuSeed = deltaMuSeed;
    meta.deltaMuTarget = deltaMuTarget;
    meta.ds = ds;
    meta.scale = scale;
    meta.optsUsed = optsUsed;
    meta.segment = branchSeg;
    meta.goodIdx = goodIdx;
    meta.nOld = nOld;
    meta.nSegment = numel(branchSeg.y);
    meta.nAdded = max(numel(branchSeg.y) - 2, 0);
    meta.nNew = numel(branchMaster.y);

    maybe_save(branchMaster, meta, extOpts.saveFile);
end

function validate_branch(branchIn, N)
    if ~isstruct(branchIn)
        error('branchIn must be a struct.');
    end
    needed = {'y','z','mu','T','rho','info','success','dsUsed'};
    for i = 1:numel(needed)
        if ~isfield(branchIn, needed{i})
            error('branchIn is missing field ''%s''.', needed{i});
        end
    end
    if numel(branchIn.y) < 2
        error('branchIn must contain at least two members.');
    end
    if nargin > 1 && ~isempty(branchIn.y)
        y1 = branchIn.y{1};
        if ~isempty(y1) && numel(y1(:)) ~= 6*N + 3
            error('branchIn.y{1} does not have length 6*N + 3.');
        end
    end
end

function [k0, k1, skippedIdx, goodIdx] = select_seed_indices(branchMaster, tailSkip, side)
    if isfield(branchMaster, 'success') && ~isempty(branchMaster.success)
        goodIdx = find(branchMaster.success(:));
    else
        goodIdx = (1:numel(branchMaster.y)).';
    end

    if numel(goodIdx) < tailSkip + 2
        error('Not enough converged branch members to restart continuation.');
    end

    if side == "tail"
        k1 = goodIdx(end - tailSkip);
        k0 = goodIdx(end - tailSkip - 1);
        skippedIdx = goodIdx((end - tailSkip + 1):end);
    else
        k0 = goodIdx(1 + tailSkip);
        k1 = goodIdx(2 + tailSkip);
        skippedIdx = goodIdx(1:tailSkip);
    end
end

function continueFcn = pick_continue_function()
    if exist('gmos_continue_fixT_mu_palc_weighted_bidirectional', 'file') == 2
        continueFcn = @gmos_continue_fixT_mu_palc_weighted_bidirectional;
    elseif exist('gmos_continue_fixT_mu_palc_weighted', 'file') == 2
        continueFcn = @gmos_continue_fixT_mu_palc_weighted;
    else
        error(['Could not find gmos_continue_fixT_mu_palc_weighted_bidirectional ', ...
               'or gmos_continue_fixT_mu_palc_weighted on the MATLAB path.']);
    end
end

function branchOut = merge_branch_segment(branchIn, branchSeg, k0, k1, side)
    if numel(branchSeg.y) < 2
        error('branchSeg must contain at least its two seed members.');
    end

    branchOut = struct();

    if side == "tail"
        branchOut.y = [tocol(branchIn.y(1:k1)); tocol(branchSeg.y(3:end))];
        branchOut.z = [tocol(branchIn.z(1:k1)); tocol(branchSeg.z(3:end))];
        branchOut.mu = [tocol(branchIn.mu(1:k1)); tocol(branchSeg.mu(3:end))];
        branchOut.T = [tocol(branchIn.T(1:k1)); tocol(branchSeg.T(3:end))];
        branchOut.rho = [tocol(branchIn.rho(1:k1)); tocol(branchSeg.rho(3:end))];
        branchOut.info = [tocol(branchIn.info(1:k1)); tocol(branchSeg.info(3:end))];
        branchOut.success = [tocol(branchIn.success(1:k1)); tocol(branchSeg.success(3:end))];
        branchOut.dsUsed = [tocol(branchIn.dsUsed(1:k1)); tocol(branchSeg.dsUsed(3:end))];
        if isfield(branchIn, 'EGMOS') && ~isempty(branchIn.EGMOS)
            branchOut.EGMOS = [tocol(branchIn.EGMOS(1:k1)); tocol(branchSeg.EGMOS(3:end))];
        else
            branchOut.EGMOS = [];
        end
    else
        branchOut.y = [tocol(branchSeg.y(1:end-2)); tocol(branchIn.y(k0:end))];
        branchOut.z = [tocol(branchSeg.z(1:end-2)); tocol(branchIn.z(k0:end))];
        branchOut.mu = [tocol(branchSeg.mu(1:end-2)); tocol(branchIn.mu(k0:end))];
        branchOut.T = [tocol(branchSeg.T(1:end-2)); tocol(branchIn.T(k0:end))];
        branchOut.rho = [tocol(branchSeg.rho(1:end-2)); tocol(branchIn.rho(k0:end))];
        branchOut.info = [tocol(branchSeg.info(1:end-2)); tocol(branchIn.info(k0:end))];
        branchOut.success = [tocol(branchSeg.success(1:end-2)); tocol(branchIn.success(k0:end))];
        branchOut.dsUsed = [tocol(branchSeg.dsUsed(1:end-2)); tocol(branchIn.dsUsed(k0:end))];
        if isfield(branchIn, 'EGMOS') && ~isempty(branchIn.EGMOS)
            branchOut.EGMOS = [tocol(branchSeg.EGMOS(1:end-2)); tocol(branchIn.EGMOS(k0:end))];
        else
            branchOut.EGMOS = [];
        end
    end

    if isfield(branchSeg, 'scale')
        branchOut.scale = branchSeg.scale;
    elseif isfield(branchIn, 'scale')
        branchOut.scale = branchIn.scale;
    end
end

function meta = build_meta_empty(branchIn, side, extOpts)
    meta = struct();
    meta.side = char(side);
    meta.seedIdx = [];
    meta.skippedIdx = [];
    meta.deltaMuSeed = NaN;
    meta.deltaMuTarget = NaN;
    meta.ds = NaN;
    meta.scale = [];
    meta.optsUsed = [];
    meta.segment = [];
    meta.nOld = numel(branchIn.y);
    meta.nSegment = 0;
    meta.nAdded = 0;
    meta.nNew = numel(branchIn.y);
    meta.note = 'No extension requested because extOpts.nSteps = 0.';
    meta.extOpts = extOpts;
end

function x = tocol(x)
    x = x(:);
end

function maybe_save(branchMaster, meta, saveFile)
    if nargin < 3 || isempty(saveFile)
        return
    end
    save(saveFile, 'branchMaster', 'meta', '-v7.3');
end
