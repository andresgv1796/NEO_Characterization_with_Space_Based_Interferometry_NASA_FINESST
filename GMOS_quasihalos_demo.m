
clear
clc
%% 

set(0,'DefaultFigureWindowStyle','normal');

set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
set(groot,'DefaultAxesFontSize',20);
set(groot,'DefaultTextFontSize',22);

% -------------------------
% Constants / seed orbit
% -------------------------
mu_EM_JPL = 1.215058560962404e-2;
mu = mu_EM_JPL;

---------------------
orbit_index = 200 ;

po = pick_po_from_JPL('L1_northern_halo_orbits_JPL_IC_Earth_Moon.csv', orbit_index);
xPO0 = po.x0(:);
TPO  = po.T;

 %% 
% ------------------------------------------------------------
% Choose target Jacobi constant
% ------------------------------------------------------------
% If you already have JC_Lyap in workspace, use that:
% Ctarget = JC_Lyap;

% Otherwise use the periodic orbit Jacobi constant:
Ctarget = Jacobi_Constant(xPO0.', mu);

fprintf('Target Jacobi constant = %.16f\n', Ctarget);
%% 

% -------------------------
% Quick trajectory plot of seed PO
% -------------------------
[t_eval, Y_eval] = IntegrateCR3BP_ODE89(xPO0, [0, TPO], mu, 1e3);

figure;
hold on; grid on; axis equal;
plot3(Y_eval(:,1), Y_eval(:,2), Y_eval(:,3), 'LineWidth', 1.2);
xlabel('x');
ylabel('y');
zlabel('z');
title('Seed periodic orbit propagation');

view([45,45,20])

% -------------------------
% ODE options and closure check
% -------------------------
ode_opts = odeset('RelTol',1e-12,'AbsTol',1e-12);

 
fprintf('PO closure error = %.3e\n', norm(Y_eval(end,1:6)' - xPO0));
%% 

% -------------------------
% GMOS setup
% -------------------------
N = 51;

% IMPORTANT:
% use 'vertical' if you want the vertical family
% use 'planar'   if you want the planar family
mode = 'vertical';

% Two amplitudes so PALC has two corrected starting members
K0 = 8e-3;
K1 = 1e-2;

demo0 = gmos_demo_invariance_only(mu, xPO0, TPO, N, K0, false, mode);
demo1 = gmos_demo_invariance_only(mu, xPO0, TPO, N, K1, false, mode);

init0 = demo0.init;
init1 = demo1.init;
%% 

% -------------------------
% First reference for phase conditions
% -------------------------
Xref   = init0.X0;
Tref   = init0.T0;
rhoref = init0.rho0;

% Initial unknown vectors z = [X0(:); T; rho]
z0 = gmos_pack_z(init0.X0, init0.T0, init0.rho0);
z1 = gmos_pack_z(init1.X0, init1.T0, init1.rho0);

optsN = struct();
optsN.maxIter = 16;
optsN.tolInf  = 1e-10;
optsN.verbose = true;

% -------------------------
% Correct first seed with fixed C
% -------------------------
%[z, info] = gmos_newton_fixC_analytic(z0, N, mu, Xref, Tref, rhoref, Ctarget, ode_opts, opts)

%function [z, info] = gmos_newton_fixT_analytic(z0, N, mu, Xref, Tref, rhoref, Ttarget, ode_opts, opts)

[z0_sol, info0] = gmos_newton_fixC_analytic( ...
    z0, N, mu, Xref, Tref, rhoref, Ctarget, ode_opts, optsN);

fprintf('seed 0 final ||F||inf = %.3e\n', info0.normInf(end));

[Xref1, Tref1, rhoref1] = gmos_unpack_z(z0_sol, N);

% -------------------------
% Correct second seed with fixed C
% -------------------------
[z1_sol, info1] = gmos_newton_fixC_analytic( ...
    z1, N, mu, Xref1, Tref1, rhoref1, Ctarget, ode_opts, optsN);

fprintf('seed 1 final ||F||inf = %.3e\n', info1.normInf(end));

% -------------------------
% Check corrected seeds
% -------------------------
[X0a, Ta, rhoa] = gmos_unpack_z(z0_sol, N);
[X0b, Tb, rhob] = gmos_unpack_z(z1_sol, N);

Ca = mean(Jacobi_Constant(X0a.', mu));
Cb = mean(Jacobi_Constant(X0b.', mu));

fprintf('seed 0: T = %.16f, rho = %.16f, mean(C) error = %.3e\n', ...
    Ta, rhoa, abs(Ca - Ctarget));
fprintf('seed 1: T = %.16f, rho = %.16f, mean(C) error = %.3e\n', ...
    Tb, rhob, abs(Cb - Ctarget));

% -------------------------
% Continue family with PALC
% -------------------------
ds     = 1e-2;
nSteps = 3;

if isfile('halo_fam_fixC.mat')
    load('halo_fam_fixC.mat','fam');
else
    fam = gmos_continue_fixC_family( ...
        z0_sol, z1_sol, nSteps, ds, N, mu, Ctarget, ode_opts, optsN);
   % save('halo_fam_fixC.mat','fam')
end

beep

fprintf('Computed %d family members.\n', numel(fam.z));
%% 

mu0=mu_EM_JPL
mu1=1.000020*mu_EM_JPL
z_start=fam.z{end}
Ttarget=fam.T(end)


mu1=1.020*mu_EM_JPL
%muVals =  linspace(mu0, mu1, 60);
muVals =  linspace(mu0, mu1, 30);
step_size=muVals(end)-muVals(end-1)

branch = gmos_continue_fixT_mu_simple( ...
    z_start, muVals, N, Ttarget, ode_opts, optsN);
%% 

figure;
plot(branch.mu, branch.rho, 'o-'); grid on;
xlabel('$\mu$'); ylabel('$\rho$');
title('$\rho$ along the $\mu$ branch');

figure;
semilogy(branch.mu, branch.EGMOS, 'o-'); grid on;
xlabel('$\mu$'); ylabel('$E_{GMOS}$');
title('GMOS accuracy along the $\mu$ branch');
%% 

k0 = numel(branch.z) - 1;
k1 = numel(branch.z);

% k0=1
% k0=2

z0_sol = branch.z{k0};
z1_sol = branch.z{k1};

mu0_sol = branch.mu(k0);
mu1_sol = branch.mu(k1);

y0_sol=[]
y1_sol=[]

y0_sol = [z0_sol(:); mu0_sol];
y1_sol = [z1_sol(:); mu1_sol];
% Rebuild weighted tangent and choose a more conservative restart step
scale = gmos_build_palc_scaling_mu(y0_sol, y1_sol, N, struct());

dy_seed_scaled = scale.w .* (y1_sol - y0_sol);
t_hat0 = dy_seed_scaled / norm(dy_seed_scaled);

deltaMu_seed = y1_sol(end) - y0_sol(end);

% Restart more conservatively than before
deltaMu_target = 0.25 * deltaMu_seed;
ds = (scale.w(end) * deltaMu_target) / t_hat0(end);

fprintf('Restart ds = %.6e\n', ds);
fprintf('Restart deltaMu_target = %.6e\n', deltaMu_target);
%% 

% Make restart options a bit safer

optsPALC = struct();
optsPALC.maxIter = 12;
optsPALC.tolInf = 1e-10;
optsPALC.verbose = true;
optsPALC.computeEGMOS = true;
optsPALC.fdMu = 1e-8;


optsPALC_restart = optsPALC;
optsPALC_restart.maxBacktrack = 14;
optsPALC_restart.backtrackFactor = 0.5;
optsPALC_restart.growFactor = 1.05;
optsPALC_restart.enforceMuDirection = true;

nSteps_restart = 3;
%% 

tic
branchPALC_more = gmos_continue_fixT_mu_palc_weighted( ...
    y0_sol, y1_sol, nSteps_restart, ds, N, Ttarget, ode_opts, optsPALC_restart);

%% 
 


%% 




extOpts = struct();
extOpts.tailSkip = 1;          % use 1, 2, or 3 if the tail looks suspicious
extOpts.deltaMuScale = 5.5;   % conservative restart
extOpts.nSteps = 2;
extOpts.saveFile = 'branchMaster_muPALC_Halos.mat';
extOpts.verbose = true;
extOpts.continuationSide = 'head';


tic
[branchMaster, meta1] = gmos_extend_fixT_mu_branch_weighted_bidirectional( ...
    branchMaster, N, Ttarget, ode_opts, optsPALC_restart, extOpts);



%% 



figure;
hold on;
grid on;

idx = 1:numel(branchMaster.mu);
scatter(branchMaster.mu, branchMaster.rho, 50, idx, 'filled');
plot(branchMaster.mu, branchMaster.rho, '-', 'LineWidth', 1.0);

xlabel('$\mu$ [n.d.]');
ylabel('$\rho$ [n.d.]');
title('$(\mu,\rho)$ colored by branch index');

cb = colorbar;
cb.Label.String = 'Branch index';
cb.Label.Interpreter = 'latex';


figure;
hold on;
grid on;

idx = 1:numel(branchMaster.mu);

scatter(branchMaster.mu, branchMaster.EGMOS, 50, idx, 'filled');
semilogy(branchMaster.mu, branchMaster.EGMOS, '-', 'LineWidth', 1.0);

set(gca, 'YScale', 'log');

xlabel('$\mu$ [n.d.]');
ylabel('$E_{GMOS}$');
title('$E_{GMOS}$ along the $\mu$ branch, colored by branch index');

cb = colorbar;
cb.Label.String = 'Branch index';
cb.Label.Interpreter = 'latex';
%% 


%% 


muAll  = branchMaster.mu(:);
rhoAll = branchMaster.rho(:);
nAll   = numel(muAll);

muMin = min(muAll);
muMax = max(muAll);
cmap  = parula(256);

% ------------------------------------------------------------
% Smart representative selection
% ------------------------------------------------------------
nShowBase = 10;   % number of roughly uniform-in-mu samples

muTargets = linspace(muMin, muMax, nShowBase);
idxMu = zeros(nShowBase,1);

for i = 1:nShowBase
    [~, idxMu(i)] = min(abs(muAll - muTargets(i)));
end

% Include endpoints
idxSpecial = [1; nAll];

% Include strongest bend in the (mu,rho) branch if possible
if nAll >= 4
    s1 = diff(rhoAll) ./ diff(muAll);
    bend = abs(diff(s1));
    [~, ib] = max(bend);
    idxSpecial(end+1) = ib + 1;
end

idxList = unique([idxSpecial; idxMu], 'stable');

%% % ------------------------------------------------------------
% Plot selected full tori
% ------------------------------------------------------------
figure;
ax = axes;
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');
view(ax, 3);

for a = 1:numel(idxList)
    k = idxList(a);

    z_k  = branchMaster.z{k};
    mu_k = branchMaster.mu(k);

    cIdx = 1 + round((size(cmap,1)-1) * (mu_k - muMin) / max(muMax - muMin, eps));
    cIdx = max(1, min(size(cmap,1), cIdx));
    thisColor = cmap(cIdx,:);

    traj_k = gmos_propagate_qpos_orbits( ...
        z_k, N, mu_k, 1:5:N, 2, 600, ode_opts);

    style_k = struct( ...
        'mode', 'rgb', ...
        'rgb', thisColor, ...
        'alpha', 0.50, ...
        'lineWidth', 1.1, ...
        'newFigure', false, ...
        'ax', ax);

    gmos_plot_qpos_orbits_3d(traj_k, style_k, 0);
end

xlabel('$x$ [n.d.]');
ylabel('$y$ [n.d.]');
zlabel('$z$ [n.d.]');
title('Representative full tori along the $\mu$ continuation branch, colored by $\mu$');

colormap(ax, cmap);
caxis(ax, [muMin muMax]);
cb = colorbar(ax);
cb.Label.String = '$\mu$';
cb.Label.Interpreter = 'latex';

%% % Pick one member of the continuation branch
k = 11;   % change as needed

z_k  = branchMaster.z{k};
mu_k = branchMaster.mu(k);

% Propagate enough curve samples to visualize one torus
traj_k = gmos_propagate_qpos_orbits( ...
    z_k, N, mu_k, 1:5:N, 17, 400, ode_opts);

% Plot that torus
figure;
ax = axes;
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');
view(ax, 3);

style_k = struct( ...
    'mode', 'rgb', ...
    'rgb', [0 0.4470 0.7410], ...
    'alpha', 0.9, ...
    'lineWidth', 1.1, ...
    'newFigure', false, ...
    'ax', ax);

gmos_plot_qpos_orbits_3d(traj_k, style_k, 1);

xlabel('$x$ [n.d.]','Interpreter','latex');
ylabel('$y$ [n.d.]','Interpreter','latex');
zlabel('$z$ [n.d.]','Interpreter','latex');
title(sprintf('One torus from continuation, k = %d, \\mu = %.8g', k, mu_k), ...
    'Interpreter', 'latex');
%% 

nB = numel(branchMaster.z);

xMin = nan(nB,1);
xMax = nan(nB,1);
yMax = nan(nB,1);
zMin = nan(nB,1);
zMax = nan(nB,1);

parfor k = 1:nB
    z_k  = branchMaster.z{k};
    mu_k = branchMaster.mu(k);

    traj_k = gmos_propagate_qpos_orbits( ...
        z_k, N, mu_k, 1:2:N, 1, 600, ode_opts);

    X = reshape(squeeze(traj_k.Y(1,:,:)), [], 1);
    Y = reshape(squeeze(traj_k.Y(2,:,:)), [], 1);
    Z = reshape(squeeze(traj_k.Y(3,:,:)), [], 1);

    xMin(k) = min(X);
    xMax(k) = max(X);
    yMax(k) = max(abs(Y));
    zMin(k) = min(abs(Z));
    zMax(k) = max(abs(Z));
end

zSpan = zMax - zMin;
%% 

muVals = branchMaster.mu(:);
idx    = (1:nB).';

xMid  = 0.5*(xMin + xMax);
xSpan = xMax - xMin;
zSpan = zMax - zMin;

figure;
% plot(muVals, xMin,  'o-', 'DisplayName', '$x_{\min}$'); hold on;
% plot(muVals, xMax,  'o-', 'DisplayName', '$x_{\max}$');
% plot(muVals, xMid,  'o-', 'DisplayName', '$x_{\mathrm{mid}}$');
%plot(muVals, xSpan, 'o-', 'DisplayName', '$x_{\mathrm{span}}$');
plot(muVals, zSpan, 'o-', 'DisplayName', '$z_{\mathrm{span}}$');
grid on;
xlabel('$\mu$ [n.d.]');
ylabel('$z$ size [n.d.]');
title('$z$ dimensions along the branch');
%legend('Interpreter','latex','Location','best');
%% 


figure;
% plot(muVals, xMin,  'o-', 'DisplayName', '$x_{\min}$'); hold on;
% plot(muVals, xMax,  'o-', 'DisplayName', '$x_{\max}$');
% plot(muVals, xMid,  'o-', 'DisplayName', '$x_{\mathrm{mid}}$');
plot(muVals, xSpan, 'o-', 'DisplayName', '$x_{\mathrm{span}}$');
%plot(muVals, zSpan, 'o-', 'DisplayName', '$z_{\mathrm{span}}$');
grid on;
xlabel('$\mu$ [n.d.]');
ylabel('$x$ size [n.d.]');
title('$x$ dimensions along the branch');
%legend('Interpreter','latex','Location','best');
% 
% 
 figure;
plot(muVals, yMax, 'o-', 'DisplayName', '$\max |y|$'); hold on;
plot(muVals, zMax, 'o-', 'DisplayName', '$\max |z|$');
grid on;
xlabel('$\mu$ [n.d.]');
ylabel('Extent [n.d.]');
title('Lateral and vertical extent along the branch');
legend('Interpreter','latex','Location','best');
%% 


%Testing with one of the largest halos as seed

%% 

 load('fam_fixT_halo.mat','fam')

% Plot one torus from the last family member
z = fam.z{100};

%z = fam.z{10};

jList = 1:1:N;
nPeriods = 2;
nSavePerPeriod = 400;

traj = gmos_propagate_qpos_orbits( ...
    z, N, ...
    mu_EM_JPL, jList, nPeriods, nSavePerPeriod, ode_opts);

hold on

gmos_plot_qpos_orbits_3d( ...
    traj, struct('mode','rgb','rgb',[0.7 0.0 0.1],'alpha',0.08,'lineWidth',1.2),1);

plot_Moon_nd(1e3,60,1);


%% 

% gmos_plot_qpos_orbits_3d( ...
%     traj, struct('mode','rgb','rgb',[0.7 0.0 0.1],'alpha',0.1,'lineWidth',1.2),0);
plot_Moon_nd(1e3,1,1);
%plot_Earth_nd(1e3,1,1);

title('Halo torus (GMOS)')
%% 

optsN = struct();
optsN.maxIter = 16;
optsN.tolInf  = 1e-10;
optsN.verbose = true;

mu0=mu_EM_JPL
mu1=1.000020*mu_EM_JPL
z_start=fam.z{end}
Ttarget=fam.T(end)


mu1=1.00020*mu_EM_JPL
%muVals =  linspace(mu0, mu1, 60);
muVals =  linspace(mu0, mu1, 2);
step_size=muVals(end)-muVals(end-1)

branch = gmos_continue_fixT_mu_simple( ...
    z_start, muVals, N, Ttarget, ode_opts, optsN);
%% 

figure;
plot(branch.mu, branch.rho, 'o-'); grid on;
xlabel('$\mu$'); ylabel('$\rho$');
title('$\rho$ along the $\mu$ branch');

figure;
semilogy(branch.mu, branch.EGMOS, 'o-'); grid on;
xlabel('$\mu$'); ylabel('$E_{GMOS}$');
title('GMOS accuracy along the $\mu$ branch');
%% 

k0 = numel(branch.z) - 1;
k1 = numel(branch.z);

% k0=1
% k0=2

z0_sol = branch.z{k0};
z1_sol = branch.z{k1};

mu0_sol = branch.mu(k0);
mu1_sol = branch.mu(k1);

y0_sol=[]
y1_sol=[]

y0_sol = [z0_sol(:); mu0_sol];
y1_sol = [z1_sol(:); mu1_sol];
% Rebuild weighted tangent and choose a more conservative restart step
scale = gmos_build_palc_scaling_mu(y0_sol, y1_sol, N, struct());

dy_seed_scaled = scale.w .* (y1_sol - y0_sol);
t_hat0 = dy_seed_scaled / norm(dy_seed_scaled);

deltaMu_seed = y1_sol(end) - y0_sol(end);

% Restart more conservatively than before
deltaMu_target = 0.25 * deltaMu_seed;
ds = (scale.w(end) * deltaMu_target) / t_hat0(end);

fprintf('Restart ds = %.6e\n', ds);
fprintf('Restart deltaMu_target = %.6e\n', deltaMu_target);
%% 

% Make restart options a bit safer

optsPALC = struct();
optsPALC.maxIter = 12;
optsPALC.tolInf = 1e-10;
optsPALC.verbose = true;
optsPALC.computeEGMOS = true;
optsPALC.fdMu = 1e-8;


optsPALC_restart = optsPALC;
optsPALC_restart.maxBacktrack = 14;
optsPALC_restart.backtrackFactor = 0.5;
optsPALC_restart.growFactor = 1.05;
optsPALC_restart.enforceMuDirection = true;

nSteps_restart = 3;
%% 

tic
branchPALC_large = gmos_continue_fixT_mu_palc_weighted( ...
    y0_sol, y1_sol, nSteps_restart, ds, N, Ttarget, ode_opts, optsPALC_restart);
%% 


load('branchMaster_muPALC_Halos_test')

 %   'tail' extends from the large-index end
%   'head' extends from the small-index end
extOpts = struct();
extOpts.tailSkip         = 1;   % skip suspicious tail points if needed
extOpts.deltaMuScale     = 1.0; % restart step scaling relative to local seed
extOpts.nSteps           = 15;   % short test extension
extOpts.saveFile         = 'branchMaster_muPALC_Halos_test.mat';
extOpts.verbose          = true;
extOpts.continuationSide = 'head';


[branchMaster_large, meta1] = gmos_extend_fixT_mu_branch_weighted_bidirectional( ...
    branchMaster_large, N, Ttarget0, ode_opts, optsPALC, extOpts);

beep

%% -------------------------
%  Choose which branch to plot
%  -------------------------
%
% Right now we plot the startup branch.
%branchToPlot = branchMaster0;

% If instead you want to plot the extended branch, use:
 branchToPlot = branchMaster_large;

% -------------------------
%  Plot rho versus mu, colored by branch index
%  -------------------------
figure;
hold on;
grid on;
box on;

idx = 1:numel(branchToPlot.mu);

scatter(branchToPlot.mu, branchToPlot.rho, 50, idx, 'filled');
plot(branchToPlot.mu, branchToPlot.rho, '-', 'LineWidth', 1.0);

xlabel('$\mu$ [n.d.]');
ylabel('$\rho$ [n.d.]');
title('$(\mu,\rho)$ along the continuation branch');

cb = colorbar;
cb.Label.String = 'Branch index';
cb.Label.Interpreter = 'latex';

% -------------------------
%  Plot EGMOS versus mu, colored by branch index
%  -------------------------
figure;
hold on;
grid on;
box on;

idx = 1:numel(branchToPlot.mu);

% Plot the points first so the color encodes branch order.
scatter(branchToPlot.mu, branchToPlot.EGMOS, 50, idx, 'filled');

% Use a standard line and then set the axis to log scale.
% This is cleaner than mixing scatter with semilogy.
plot(branchToPlot.mu, branchToPlot.EGMOS, '-', 'LineWidth', 1.0);
set(gca, 'YScale', 'log');

xlabel('$\mu$ [n.d.]');
ylabel('$E_{GMOS}$');
title('$E_{GMOS}$ along the $\mu$ continuation branch');

cb = colorbar;
cb.Label.String = 'Branch index';
cb.Label.Interpreter = 'latex';
muAll  = branchToPlot.mu(:);
rhoAll = branchToPlot.rho(:);
nAll   = numel(muAll);

muMin = min(muAll);
muMax = max(muAll);
cmap  = parula(256);
%% 

% ------------------------------------------------------------
% Smart representative selection
% ------------------------------------------------------------
nShowBase = 50;   % number of roughly uniform-in-mu samples

muTargets = linspace(muMin, muMax, nShowBase);
idxMu = zeros(nShowBase,1);

for i = 1:nShowBase
    [~, idxMu(i)] = min(abs(muAll - muTargets(i)));
end

% Include endpoints
idxSpecial = [1; nAll];

% Include strongest bend in the (mu,rho) branch if possible
if nAll >= 4
    s1 = diff(rhoAll) ./ diff(muAll);
    bend = abs(diff(s1));
    [~, ib] = max(bend);
    idxSpecial(end+1) = ib + 1;
end

idxList = unique([idxSpecial; idxMu], 'stable');

% ------------------------------------------------------------
% Plot selected full tori
% ------------------------------------------------------------
figure;
ax = axes;
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');
view(ax, 3);

for a = 1:numel(idxList)
    k = idxList(a);

    z_k  = branchToPlot.z{k};
    mu_k = branchToPlot.mu(k);

    cIdx = 1 + round((size(cmap,1)-1) * (mu_k - muMin) / max(muMax - muMin, eps));
    cIdx = max(1, min(size(cmap,1), cIdx));
    thisColor = cmap(cIdx,:);

    traj_k = gmos_propagate_qpos_orbits( ...
        z_k, N, mu_k, 1:5:N, 1, 600, ode_opts);

    style_k = struct( ...
        'mode', 'rgb', ...
        'rgb', thisColor, ...
        'alpha', 0.10, ...
        'lineWidth', 1.1, ...
        'newFigure', false, ...
        'ax', ax);

    gmos_plot_qpos_orbits_3d(traj_k, style_k, 0);
end

xlabel('$x$ [n.d.]');
ylabel('$y$ [n.d.]');
zlabel('$z$ [n.d.]');
title('Representative full tori along the $\mu$ continuation branch, colored by $\mu$');

colormap(ax, cmap);
caxis(ax, [muMin muMax]);
cb = colorbar(ax);
cb.Label.String = '$\mu$';
cb.Label.Interpreter = 'latex';

%% % Representative invariant curves from branchMaster
% Assumes:
%   [branchMaster, meta1] = gmos_extend_fixT_mu_branch_weighted(...)
% ============================================================

branchUse = branchToPlot;

% ------------------------------------------------------------
% Choose representative members
% ------------------------------------------------------------
nShow = 300;
idxList = unique(round(linspace(1, numel(branchUse.z), min(nShow, numel(branchUse.z)))));

muMin = min(branchUse.mu);
muMax = max(branchUse.mu);
cmap  = parula(256);

curveColors = zeros(numel(idxList), 3);
for a = 1:numel(idxList)
    mu_k = branchUse.mu(idxList(a));
    cIdx = 1 + round((size(cmap,1)-1) * (mu_k - muMin) / max(muMax - muMin, eps));
    cIdx = max(1, min(size(cmap,1), cIdx));
    curveColors(a,:) = cmap(cIdx,:);
end

% ------------------------------------------------------------
% x-y invariant curves
% ------------------------------------------------------------
figure;
ax1 = axes;
hold(ax1, 'on');
grid(ax1, 'on');
axis(ax1, 'equal');

for a = 1:numel(idxList)
    k = idxList(a);

    z_k = branchUse.z{k};
    [X0_k, ~, ~] = gmos_unpack_z(z_k, N);

    Xc = [X0_k, X0_k(:,1)];
    plot3(Xc(1,:), Xc(2,:), Xc(3,:),...
        'LineWidth', 1.6, ...
        'Color', curveColors(a,:));
end

xlabel('$x$ [n.d.]');
ylabel('$y$ [n.d.]');
title('Representative invariant curves in the $x$-$y$ plane, colored by $\mu$');

colormap(ax1, cmap);
caxis(ax1, [muMin muMax]);
cb1 = colorbar(ax1);
cb1.Label.String = '$\mu$';
cb1.Label.Interpreter = 'latex';

% ------------------------------------------------------------
% y-v_y invariant curves
% ------------------------------------------------------------
figure;
ax2 = axes;
hold(ax2, 'on');
grid(ax2, 'on');

for a = 1:numel(idxList)
    k = idxList(a);

    z_k = branchUse.z{k};
    [X0_k, ~, ~] = gmos_unpack_z(z_k, N);

    Xc = [X0_k, X0_k(:,1)];
    plot(Xc(2,:), Xc(5,:), ...
        'LineWidth', 1.6, ...
        'Color', curveColors(a,:));
end

xlabel('$y$ [n.d.]');
ylabel('$\dot{y}$ [n.d.]');
title('Representative invariant curves in the $y$-$\dot{y}$ plane, colored by $\mu$');

colormap(ax2, cmap);
caxis(ax2, [muMin muMax]);
cb2 = colorbar(ax2);
cb2.Label.String = '$\mu$';
cb2.Label.Interpreter = 'latex';

