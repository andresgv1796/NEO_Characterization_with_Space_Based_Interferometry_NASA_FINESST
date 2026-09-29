function gmos_plot_qpos_orbits_3d(traj, style, t_i)
%GMOS_PLOT_QPOS_ORBITS_3D Plot 3D trajectories with selectable styling.
%
% Allowed calls:
%   gmos_plot_qpos_orbits_3d(traj)
%   gmos_plot_qpos_orbits_3d(traj, style)
%   gmos_plot_qpos_orbits_3d(traj, t_i)
%   gmos_plot_qpos_orbits_3d(traj, style, t_i)
%
% Optional style fields:
%   style.mode            = 'index' | 'gray' | 'rgb'
%   style.lineWidth       = scalar
%   style.alpha           = scalar in [0,1]
%   style.gray            = scalar gray level
%   style.rgb             = 1x3 RGB color
%   style.newFigure       = true/false
%   style.ax              = axes handle
%
% Optional ring styling fields (used only if t_i is provided):
%   style.showRingLine    = true/false
%   style.showRingScatter = true/false
%   style.ringLineWidth   = scalar
%   style.ringLineColor   = 1x3 RGB color
%   style.ringMarkerSize  = scalar
%   style.ringMarkerAlpha = scalar in [0,1]

    % -------------------------
    % Flexible input parsing
    % -------------------------
    if nargin < 2
        style = struct();
        t_i = [];
    elseif nargin < 3
        if isstruct(style) || isempty(style)
            t_i = [];
        else
            t_i = style;
            style = struct();
        end
    end

    if isempty(style)
        style = struct();
    end
    if ~isstruct(style)
        error('Second input must be a style struct or a scalar time t_i.');
    end

    % -------------------------
    % Defaults
    % -------------------------
    if ~isfield(style, 'mode') || isempty(style.mode)
        style.mode = 'index';
    end
    if ~isfield(style, 'lineWidth') || isempty(style.lineWidth)
        style.lineWidth = 1.0;
    end
    if ~isfield(style, 'alpha') || isempty(style.alpha)
        style.alpha = 0.25;
    end
    if ~isfield(style, 'newFigure') || isempty(style.newFigure)
        style.newFigure = true;
    end
    if ~isfield(style, 'ax')
        style.ax = [];
    end

    % Ring defaults
    if ~isfield(style, 'showRingLine') || isempty(style.showRingLine)
        style.showRingLine = false;
    end
    if ~isfield(style, 'showRingScatter') || isempty(style.showRingScatter)
        style.showRingScatter = true;
    end
    if ~isfield(style, 'ringLineWidth') || isempty(style.ringLineWidth)
        style.ringLineWidth = 1.5;
    end
    if ~isfield(style, 'ringMarkerSize') || isempty(style.ringMarkerSize)
        style.ringMarkerSize = 36;
    end
    if ~isfield(style, 'ringMarkerAlpha') || isempty(style.ringMarkerAlpha)
        style.ringMarkerAlpha = 1.0;
    end

    mode = lower(string(style.mode));
    [~, ~, M] = size(traj.Y);

    if mode == "index"
        cols = lines(max(M,1));
    elseif mode == "gray"
        if ~isfield(style, 'gray') || isempty(style.gray)
            style.gray = 0.25;
        end
        cols = repmat([style.gray style.gray style.gray], M, 1);
    elseif mode == "rgb"
        if ~isfield(style, 'rgb') || isempty(style.rgb)
            style.rgb = [0 0 0];
        end
        cols = repmat(style.rgb(:).', M, 1);
    else
        error("style.mode must be 'index', 'gray', or 'rgb'.");
    end

    % -------------------------
    % Axes / figure handling
    % -------------------------
    if ~isempty(style.ax)
        ax = style.ax;
        if ~ishandle(ax) || ~strcmp(get(ax, 'Type'), 'axes')
            error('style.ax must be a valid axes handle.');
        end
        axes(ax);
    elseif style.newFigure
        figure;
        ax = gca;
    else
        ax = gca;
    end

    hold(ax, 'on');
    grid(ax, 'on');
    axis(ax, 'equal');
    view(ax, 3);
    xlabel(ax, '$x$ [n.d.]');
    ylabel(ax, '$y$ [n.d.]');
    zlabel(ax, '$z$ [n.d.]');

    % -------------------------
    % Draw trajectories
    % -------------------------
    for m = 1:M
        r = squeeze(traj.Y(1:3,:,m));

        patch(ax, ...
            'XData', [r(1,:) nan], ...
            'YData', [r(2,:) nan], ...
            'ZData', [r(3,:) nan], ...
            'FaceColor', 'none', ...
            'EdgeColor', cols(m,:), ...
            'EdgeAlpha', style.alpha, ...
            'LineWidth', style.lineWidth);
    end

    % -------------------------
    % Optional ring at fixed time t_i
    % -------------------------
    if ~isempty(t_i)
        if ~isscalar(t_i) || ~isnumeric(t_i) || ~isfinite(t_i)
            error('t_i must be a finite numeric scalar.');
        end

        if isfield(traj, 't')
            tvec = traj.t;
        elseif isfield(traj, 'tout')
            tvec = traj.tout;
        else
            error('traj must contain a time vector field, such as traj.t or traj.tout.');
        end

        [~, k_i] = min(abs(tvec - t_i));

        x_ring = squeeze(traj.Y(1,k_i,:));
        y_ring = squeeze(traj.Y(2,k_i,:));
        z_ring = squeeze(traj.Y(3,k_i,:));

        % Ring line
        if style.showRingLine
            if isfield(style, 'ringLineColor') && ~isempty(style.ringLineColor)
                ringLineColor = style.ringLineColor;
            else
                % If all orbit colors are the same, use that color.
                % Otherwise use the mean color as a neutral representative line.
                ringLineColor = mean(cols, 1);
            end

            plot3(ax, x_ring, y_ring, z_ring, '-', ...
                'LineWidth', style.ringLineWidth, ...
                'Color', ringLineColor);
        end

        % Ring scatter points, each matching its orbit color
        if style.showRingScatter
            hsc = scatter3(ax, x_ring, y_ring, z_ring, ...
                style.ringMarkerSize, cols, 'filled');

            if isprop(hsc, 'MarkerFaceAlpha')
                hsc.MarkerFaceAlpha = style.ringMarkerAlpha;
            end
            if isprop(hsc, 'MarkerEdgeAlpha')
                hsc.MarkerEdgeAlpha = style.ringMarkerAlpha;
            end
        end
    end
end