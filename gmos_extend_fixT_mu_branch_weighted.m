function [branchMaster, meta] = gmos_extend_fixT_mu_branch_weighted( ...
    branchMaster, N, Ttarget, ode_opts, optsPALC, extOpts)
%GMOS_EXTEND_FIXT_MU_BRANCH_WEIGHTED
% Extend an existing weighted mu-PALC branch in a systematic way.
%
% Inputs
%   branchMaster : existing accumulated branch struct
%   N            : number of curve nodes
%   Ttarget      : fixed stroboscopic period
%   ode_opts     : ODE options
%   optsPALC     : options for gmos_continue_fixT_mu_palc_weighted
%   extOpts      : extension options struct
%
% extOpts fields with defaults
%   .tailSkip        = 1
%   .deltaMuScale    = 0.25
%   .nSteps          = 100
%   .saveFile        = ''
%   .verbose         = true
%
% Output
%   branchMaster : updated accumulated branch
%   meta         : restart metadata

    if nargin < 6
        extOpts = struct();
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

    % --------------------------------------------
    % Pick restart seeds from successful branch members
    % --------------------------------------------
    if isfield(branchMaster, 'success') && ~isempty(branchMaster.success)
        goodIdx = find(branchMaster.success);
    else
        goodIdx = 1:numel(branchMaster.y);
    end

    if numel(goodIdx) < extOpts.tailSkip + 2
        error('Not enough converged branch members to restart continuation.');
    end

    k1 = goodIdx(end - extOpts.tailSkip);
    k0 = goodIdx(end - extOpts.tailSkip - 1);

    y0_sol = branchMaster.y{k0};
    y1_sol = branchMaster.y{k1};

    % --------------------------------------------
    % Rebuild weighted tangent and restart step
    % --------------------------------------------
    scale = gmos_build_palc_scaling_mu(y0_sol, y1_sol, N, struct());

    dy_seed_scaled = scale.w .* (y1_sol - y0_sol);
    t_hat0 = dy_seed_scaled / norm(dy_seed_scaled);

    deltaMu_seed = y1_sol(end) - y0_sol(end);
    deltaMu_target = extOpts.deltaMuScale * deltaMu_seed;

    ds = (scale.w(end) * deltaMu_target) / t_hat0(end);

    if extOpts.verbose
        fprintf('\n');
        fprintf('Restarting weighted mu-PALC from k0 = %d, k1 = %d\n', k0, k1);
        fprintf('mu0 = %.16f\n', y0_sol(end));
        fprintf('mu1 = %.16f\n', y1_sol(end));
        fprintf('deltaMu_seed   = %.6e\n', deltaMu_seed);
        fprintf('deltaMu_target = %.6e\n', deltaMu_target);
        fprintf('ds restart     = %.6e\n', ds);
    end

    % --------------------------------------------
    % Continue from those seeds
    % --------------------------------------------
    branchMore = gmos_continue_fixT_mu_palc_weighted( ...
        y0_sol, y1_sol, extOpts.nSteps, ds, N, Ttarget, ode_opts, optsPALC);

    % --------------------------------------------
    % Merge old branch up to k1 with new segment from 2:end
    % --------------------------------------------
    branchMaster = gmos_merge_fixT_mu_branch_segments(branchMaster, branchMore, k1);

    % --------------------------------------------
    % Save checkpoint if requested
    % --------------------------------------------
    if ~isempty(extOpts.saveFile)
        save(extOpts.saveFile, 'branchMaster', '-v7.3');
    end

    % --------------------------------------------
    % Metadata
    % --------------------------------------------
    meta = struct();
    meta.k0 = k0;
    meta.k1 = k1;
    meta.mu0 = y0_sol(end);
    meta.mu1 = y1_sol(end);
    meta.deltaMu_seed = deltaMu_seed;
    meta.deltaMu_target = deltaMu_target;
    meta.ds = ds;
    meta.nAdded = max(0, numel(branchMore.y) - 1);
end