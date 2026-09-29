function [famMaster, meta] = gmos_extend_fixC_family_adaptive_bidirectional( ...
    famMaster, N, mu, Ctarget, ode_opts, optsPALC, extOpts)
%GMOS_EXTEND_FIXC_FAMILY_ADAPTIVE_BIDIRECTIONAL
% Extend an existing fixed-C family from the head or tail using the
% adaptive PALC continuation routine, then merge automatically.
%
% extOpts fields
%   continuationSide = 'tail' or 'head'
%   sideSkip         = 1
%   dsScale          = 1.0
%   nSteps           = 100
%   saveFile         = ''
%   verbose          = true

    narginchk(7, 7);

    validate_family(famMaster, N);

    if nargin < 5 || isempty(ode_opts)
        ode_opts = [];
    end
    if nargin < 6 || isempty(optsPALC)
        optsPALC = struct();
    end
    if nargin < 7 || isempty(extOpts)
        extOpts = struct();
    end

    if ~isfield(extOpts, 'continuationSide') || isempty(extOpts.continuationSide)
        extOpts.continuationSide = 'tail';
    end
    if ~isfield(extOpts, 'sideSkip') || isempty(extOpts.sideSkip)
        if isfield(extOpts, 'tailSkip') && ~isempty(extOpts.tailSkip)
            extOpts.sideSkip = extOpts.tailSkip;
        else
            extOpts.sideSkip = 1;
        end
    end
    if ~isfield(extOpts, 'dsScale') || isempty(extOpts.dsScale)
        extOpts.dsScale = 1.0;
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

    side = lower(string(extOpts.continuationSide));
    if side ~= "tail" && side ~= "head"
        error('extOpts.continuationSide must be ''tail'' or ''head''.');
    end

    validateattributes(extOpts.sideSkip, {'numeric'}, ...
        {'scalar','integer','>=',0}, mfilename, 'extOpts.sideSkip');
    validateattributes(extOpts.dsScale, {'numeric'}, ...
        {'scalar','real','finite','nonzero'}, mfilename, 'extOpts.dsScale');
    validateattributes(extOpts.nSteps, {'numeric'}, ...
        {'scalar','integer','>=',0}, mfilename, 'extOpts.nSteps');

    nOld = numel(famMaster.z);
    if extOpts.nSteps == 0
        meta = build_meta_empty_fixC(famMaster, side, extOpts);
        maybe_save_fixC(famMaster, meta, extOpts.saveFile);
        return
    end

    [k0, k1, skippedIdx, goodIdx] = select_seed_indices_fixC(famMaster, extOpts.sideSkip, side);

    z0_sol = famMaster.z{k0};
    z1_sol = famMaster.z{k1};

    dsSeed = norm(z1_sol - z0_sol);
    if ~(isfinite(dsSeed) && dsSeed > 0)
        error('The fixed-C seed secant has zero or invalid norm.');
    end

    ds = extOpts.dsScale * dsSeed;

    if extOpts.verbose
        fprintf('Extending fixed-C family from the %s side\n', char(side));
        fprintf('Using seed indices k0 = %d, k1 = %d\n', k0, k1);
        fprintf('dsSeed   = %.6e\n', dsSeed);
        fprintf('ds       = %.6e\n', ds);
    end

    if side == "tail"
        famSeg = gmos_continue_fixC_family_adaptive( ...
            z0_sol, z1_sol, extOpts.nSteps, ds, N, mu, Ctarget, ode_opts, optsPALC);

        famMaster = merge_fixC_tail(famMaster, famSeg, k1);

    else
        % Reverse the seed order so the secant points outward from the head
        famSegRaw = gmos_continue_fixC_family_adaptive( ...
            z1_sol, z0_sol, extOpts.nSteps, ds, N, mu, Ctarget, ode_opts, optsPALC);

        famSeg = reverse_family_fixC(famSegRaw);

        famMaster = merge_fixC_head(famMaster, famSeg, k0);
    end

    meta = struct();
    meta.side = char(side);
    meta.seedIdx = [k0, k1];
    meta.skippedIdx = skippedIdx;
    meta.goodIdx = goodIdx;
    meta.dsSeed = dsSeed;
    meta.ds = ds;
    meta.nOld = nOld;
    meta.nSegment = numel(famSeg.z);
    meta.nAdded = max(numel(famSeg.z) - 2, 0);
    meta.nNew = numel(famMaster.z);

    maybe_save_fixC(famMaster, meta, extOpts.saveFile);
end


function validate_family(fam, N)
    if ~isstruct(fam)
        error('famMaster must be a struct.');
    end

    requiredFields = {'z','rho','T','Cavg','normInf','info','success'};
    for i = 1:numel(requiredFields)
        if ~isfield(fam, requiredFields{i})
            error('famMaster is missing required field "%s".', requiredFields{i});
        end
    end

    nFam = numel(fam.z);
    if nFam < 2
        error('famMaster must contain at least two family members.');
    end

    expectedLen = 6*N + 2;
    for k = 1:nFam
        if isempty(fam.z{k}) || numel(fam.z{k}) ~= expectedLen
            error('famMaster.z{%d} must have length 6*N + 2.', k);
        end
    end
end


function [k0, k1, skippedIdx, goodIdx] = select_seed_indices_fixC(famMaster, sideSkip, side)
    if isfield(famMaster, 'success') && ~isempty(famMaster.success)
        goodIdx = find(famMaster.success(:));
    else
        goodIdx = (1:numel(famMaster.z)).';
    end

    if numel(goodIdx) < sideSkip + 2
        error('Not enough converged family members to restart continuation.');
    end

    if side == "tail"
        k1 = goodIdx(end - sideSkip);
        k0 = goodIdx(end - sideSkip - 1);
        skippedIdx = goodIdx((end - sideSkip + 1):end);
    else
        k0 = goodIdx(1 + sideSkip);
        k1 = goodIdx(2 + sideSkip);
        skippedIdx = goodIdx(1:sideSkip);
    end
end


function famOut = merge_fixC_tail(famOld, famSeg, k1)
    famOut = struct();
    famOut.z       = [famOld.z(1:k1);       famSeg.z(3:end)];
    famOut.rho     = [famOld.rho(1:k1);     famSeg.rho(3:end)];
    famOut.T       = [famOld.T(1:k1);       famSeg.T(3:end)];
    famOut.Cavg    = [famOld.Cavg(1:k1);    famSeg.Cavg(3:end)];
    famOut.normInf = [famOld.normInf(1:k1); famSeg.normInf(3:end)];
    famOut.info    = [famOld.info(1:k1);    famSeg.info(3:end)];
    famOut.success = [famOld.success(1:k1); famSeg.success(3:end)];
    famOut.dsUsed  = [get_ds_fixC(famOld, k1); get_ds_fixC(famSeg, numel(famSeg.z), 3)];
end


function famOut = merge_fixC_head(famOld, famSeg, k0)
    famOut = struct();
    famOut.z       = [famSeg.z(1:end-2);       famOld.z(k0:end)];
    famOut.rho     = [famSeg.rho(1:end-2);     famOld.rho(k0:end)];
    famOut.T       = [famSeg.T(1:end-2);       famOld.T(k0:end)];
    famOut.Cavg    = [famSeg.Cavg(1:end-2);    famOld.Cavg(k0:end)];
    famOut.normInf = [famSeg.normInf(1:end-2); famOld.normInf(k0:end)];
    famOut.info    = [famSeg.info(1:end-2);    famOld.info(k0:end)];
    famOut.success = [famSeg.success(1:end-2); famOld.success(k0:end)];
    famOut.dsUsed  = [get_ds_fixC(famSeg, numel(famSeg.z)-2, 1); get_ds_fixC(famOld, numel(famOld.z), k0)];
end


function famOut = reverse_family_fixC(famIn)
    famOut = famIn;
    famOut.z       = famIn.z(end:-1:1);
    famOut.rho     = famIn.rho(end:-1:1);
    famOut.T       = famIn.T(end:-1:1);
    famOut.Cavg    = famIn.Cavg(end:-1:1);
    famOut.normInf = famIn.normInf(end:-1:1);
    famOut.info    = famIn.info(end:-1:1);
    famOut.success = famIn.success(end:-1:1);
    if isfield(famIn, 'dsUsed')
        famOut.dsUsed = famIn.dsUsed(end:-1:1);
    end
end


function ds = get_ds_fixC(fam, iEnd, iStart)
    if nargin < 3
        iStart = 1;
    end
    if isfield(fam, 'dsUsed') && ~isempty(fam.dsUsed)
        ds = fam.dsUsed(iStart:iEnd);
    else
        ds = nan(iEnd - iStart + 1, 1);
    end
end


function meta = build_meta_empty_fixC(famMaster, side, extOpts)
    meta = struct();
    meta.side = char(side);
    meta.seedIdx = [];
    meta.skippedIdx = [];
    meta.dsSeed = NaN;
    meta.ds = NaN;
    meta.nOld = numel(famMaster.z);
    meta.nSegment = 0;
    meta.nAdded = 0;
    meta.nNew = numel(famMaster.z);
    meta.note = 'No extension requested because extOpts.nSteps = 0.';
    meta.extOpts = extOpts;
end


function maybe_save_fixC(famMaster, meta, saveFile)
    if nargin < 3 || isempty(saveFile)
        return
    end
    save(saveFile, 'famMaster', 'meta', '-v7.3');
end