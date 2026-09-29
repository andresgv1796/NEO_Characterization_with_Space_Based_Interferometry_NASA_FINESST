function init = gmos_init_curve_from_center_mode(xPO0, TPO, cen, N, K, doPlot)
%GMOS_INIT_CURVE_FROM_CENTER_MODE Initialize a GMOS invariant-curve guess from a center eigenmode.
%
% init = gmos_init_curve_from_center_mode(xPO0, TPO, cen, N, K, doPlot)
%
% Inputs
%   xPO0    : 6x1 periodic orbit initial condition
%   TPO     : scalar orbit period (used as initial stroboscopic time)
%   cen     : struct with fields
%               w     : 6x1 complex center eigenvector
%               rho0  : scalar, angle of center eigenvalue
%   N       : number of curve points (theta1 samples)
%   K       : small amplitude scaling for initial circle
%   doPlot  : (optional) logical, if true makes a 2x2 diagnostic plot (default false)
%
% Output
%   init    : struct with fields
%               X0     : 6xN initial curve samples
%               T0     : scalar, initial stroboscopic time (TPO)
%               rho0   : scalar, initial rotation number
%               theta  : 1xN theta1 grid

    if nargin < 6 || isempty(doPlot)
        doPlot = false;
    end

    % Theta1 grid
    theta = (0:N-1) * (2*pi/N);

    % Center eigenvector real/imag parts define the ellipse plane
    wr = real(cen.w);
    wi = imag(cen.w);

    % Build initial invariant circle samples around the periodic orbit IC
    X0 = zeros(6, N);
    for j = 1:N
        X0(:, j) = xPO0 + K * (cos(theta(j)) * wr - sin(theta(j)) * wi);
    end

    % Pack outputs
    init = struct();
    init.X0 = X0;
    init.T0 = TPO;
    init.rho0 = cen.rho0;
    init.theta = theta;

    % Optional diagnostic plots (2x2)
    if doPlot
        % Panel definitions: [xIndex, yIndex, xlabel, ylabel, title]
        panels = { ...
            1, 2, 'x',    'y',    'x-y'; ...
            1, 4, 'x',    'xdot', 'x-xdot'; ...
            2, 5, 'y',    'ydot', 'y-ydot'; ...
            1, 5, 'x',    'ydot', 'x-ydot' ...
        };

        figure;
        for p = 1:4
            subplot(2,2,p);

            ix = panels{p,1};
            iy = panels{p,2};

            % Center marker at the periodic-orbit IC
            plot(xPO0(ix), xPO0(iy), 'kx', 'MarkerSize', 9, 'LineWidth', 2); hold on;

            % Curve samples
            plot(X0(ix,:), X0(iy,:), 'o-');

            grid on;
            axis equal;

            xlabel(panels{p,3});
            ylabel(panels{p,4});
            title(['GMOS init: ', panels{p,5}]);
        end
    end
end