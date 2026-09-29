%DEMO_GMOS_FIXT_MU_PALC_WEIGHTED_DRIVER Example driver for weighted mu continuation.
%
% Assumes you already built a naive branch in mu and want to restart with a
% robust weighted PALC corrector. Replace the placeholders below with your
% own data.

% Inputs expected in the workspace:
%   branchSimple : output of gmos_continue_fixT_mu_simple
%   N            : number of nodes
%   Ttarget      : fixed stroboscopic time
%   ode_opts     : ODE tolerances

% Example seed selection from a previously converged simple branch:
%   i0 = numel(branchSimple.y) - 2;
%   i1 = numel(branchSimple.y) - 1;
%   y0_sol = [branchSimple.z{i0}(:); branchSimple.T(i0); branchSimple.rho(i0); branchSimple.mu(i0)];
%   y1_sol = [branchSimple.z{i1}(:); branchSimple.T(i1); branchSimple.rho(i1); branchSimple.mu(i1)];
%
% In practice use gmos_pack_y_mu after unpacking z{i}.

scale = gmos_build_palc_scaling_mu(y0_sol, y1_sol, N, struct());
ds0 = norm(scale.w .* (y1_sol - y0_sol));

optsPALC = struct();
optsPALC.maxIter = 12;
optsPALC.tolInf = 1e-10;
optsPALC.verbose = true;
optsPALC.verboseStep = true;
optsPALC.computeEGMOS = true;
optsPALC.fdMu = 1e-8;
optsPALC.maxBacktrack = 6;
optsPALC.backtrackFactor = 0.5;
optsPALC.growFactor = 1.10;
optsPALC.enforceMuDirection = false;

nSteps = 20;
ds = 0.5 * ds0;

branchPALC = gmos_continue_fixT_mu_palc_weighted( ...
    y0_sol, y1_sol, nSteps, ds, N, Ttarget, ode_opts, optsPALC);

figure;
plot(branchPALC.mu, branchPALC.rho, '-o', 'LineWidth', 1.2);
grid on;
xlabel('$\mu$ [n.d.]', 'Interpreter', 'latex');
ylabel('$\rho$ [n.d.]', 'Interpreter', 'latex');
title('$\mu$ vs. $\rho$ along the weighted GMOS PALC branch', 'Interpreter', 'latex');
