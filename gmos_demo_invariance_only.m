function demo = gmos_demo_invariance_only(mu, xPO0, TPO, N, K, varargin)
%GMOS_DEMO_INVARIANCE_ONLY Minimal GMOS pipeline test: invariance residual only.
%
% demo = gmos_demo_invariance_only(mu, xPO0, TPO, N, K)
% demo = gmos_demo_invariance_only(mu, xPO0, TPO, N, K, mode)
% demo = gmos_demo_invariance_only(mu, xPO0, TPO, N, K, doPlot, mode)
%
% Inputs
%   mu    : CR3BP mass parameter
%   xPO0  : 6x1 periodic orbit initial condition
%   TPO   : scalar periodic orbit period
%   N     : number of curve points (theta1 samples)
%   K     : initial circle amplitude scale
%
% Optional inputs
%   doPlot: logical (default true)
%   mode  : 'planar' or 'vertical' (default 'planar')
%
% Output
%   demo  : struct with diagnostic fields

    % -----------------------
    % Parse optional inputs
    % -----------------------
    doPlot = true;
    mode = 'planar';

    if numel(varargin) == 1
        % If 6th arg is a string/char, treat as mode; else treat as doPlot
        if ischar(varargin{1}) || isstring(varargin{1})
            mode = char(varargin{1});
        else
            doPlot = logical(varargin{1});
        end
    elseif numel(varargin) >= 2
        doPlot = logical(varargin{1});
        mode = char(varargin{2});
    end

    % -----------------------
    % ODE tolerances
    % -----------------------
    ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

    % -----------------------
    % Center eigenpair
    % -----------------------
    tspan = [0, TPO];
    cen = gmos_center_eigpair_from_monodromy(tspan, xPO0, mu, mode);

    % -----------------------
    % Initialize curve
    % -----------------------
    init = gmos_init_curve_from_center_mode(xPO0, TPO, cen, N, K, doPlot);

    % -----------------------
    % Amplitude diagnostics
    % -----------------------
    amp = max(init.X0, [], 2) - min(init.X0, [], 2);
    [~, iMax] = max(amp);

    fprintf('Mode: %s\n', mode);
    if isfield(cen, 'proj')
        fprintf('Eigenvector projections: planar %.3e, vertical %.3e\n', ...
            cen.proj.planar, cen.proj.vertical);
    end
    fprintf('Init amplitude [x y z xd yd zd] = [%.2e %.2e %.2e %.2e %.2e %.2e]\n', amp);
    fprintf('Max init amplitude: component %d, amp = %.3e\n', iMax, amp(iMax));

    % -----------------------
    % Propagate curve
    % -----------------------
    X1 = gmos_flowmap_curve_CR3BP(init.X0, init.T0, mu, ode_opts);

    % fprintf('size(X0) = %dx%d, size(X1) = %dx%d\n', ...
    % size(init.X0,1), size(init.X0,2), size(X1,1), size(X1,2));

    % -----------------------
    % Pick sign of rho
    % -----------------------


%    X0 = init.X0;

% X0p = gmos_shift_curve_fft(X0,  init.rho0);
% X0m = gmos_shift_curve_fft(X0, -init.rho0);
% 
% fprintf('Shift self test:\n');
% fprintf('  ||X0p - X0||inf = %.3e\n', norm(X0p(:) - X0(:), inf));
% fprintf('  ||X0m - X0||inf = %.3e\n', norm(X0m(:) - X0(:), inf));
% fprintf('  ||X0p - X0m||inf = %.3e\n', norm(X0p(:) - X0m(:), inf));

% Get monodromy matrix itself from your monodromy routine
M_eigs = monodromy_CR3BP([0,TPO], xPO0, mu);
V = M_eigs(:,1:6);
D = M_eigs(:,7:12);

% Reconstruct M = V*D/V
M = V * D / V;

% Pick one sample j (or loop)
j = 1;
dx0 = init.X0(:,j) - xPO0;
dx1 = X1(:,j) - xPO0;
dx1_lin = M * dx0;

fprintf('Linearization test at j=%d:\n', j);
fprintf('  ||dx0||      = %.3e\n', norm(dx0));
fprintf('  ||dx1||      = %.3e\n', norm(dx1));
fprintf('  ||dx1_lin||  = %.3e\n', norm(dx1_lin));
fprintf('  ||dx1-dx1lin|| = %.3e\n', norm(dx1 - dx1_lin));



    X1shift_p = gmos_shift_curve_fft(X1,  init.rho0);
    X1shift_m = gmos_shift_curve_fft(X1, -init.rho0);

    Finv_p = X1shift_p - init.X0;
    Finv_m = X1shift_m - init.X0;

    np = norm(Finv_p(:), inf);
    nm = norm(Finv_m(:), inf);

    if nm < np
        init.rho0 = -init.rho0;
        X1shift = X1shift_m;
        Finv = Finv_m;
        signUsed = -1;
    else
        X1shift = X1shift_p;
        Finv = Finv_p;
        signUsed = +1;
    end

    fprintf('Inf residual: +rho %.3e, -rho %.3e, using sign %d\n', np, nm, signUsed);

    % -----------------------
    % Norms
    % -----------------------
    norm2 = norm(Finv(:), 2);
    normInf = norm(Finv(:), inf);

    fprintf('||Finv||2   = %.3e\n', norm2);
    fprintf('||Finv||inf = %.3e\n', normInf);

    % -----------------------
    % Package outputs
    % -----------------------
    demo = struct();
    demo.init = init;
    demo.cen = cen;
    demo.X1 = X1;
    demo.X1shift = X1shift;
    demo.Finv = Finv;
    demo.amp = amp;
    demo.iMax = iMax;
    demo.signUsed = signUsed;
    demo.norm2 = norm2;
    demo.normInf = normInf;
    demo.mode = mode;

    % -----------------------
    % Plots
    % -----------------------
    if doPlot
        figure;
        plot(1:N, init.X0(iMax,:), 'o-'); hold on;  
        plot(1:N, X1shift(iMax,:), '.-');
        grid on;
        xlabel('index j');
        ylabel(sprintf('state component %d', iMax));
        legend('X0', 'R_{-rho}(Phi^T(X0))');
        title('GMOS invariance demo (max-amplitude component)');

        figure;
        plot(1:N, init.X0(1,:), 'o-'); hold on;  
        plot(1:N, X1shift(1,:), '.-');
        grid on;
        xlabel('index j');
        ylabel('X component');
        legend('X0', 'R_{-rho}(Phi^T(X0))');
        title('GMOS invariance demo (X component)');
    end
end