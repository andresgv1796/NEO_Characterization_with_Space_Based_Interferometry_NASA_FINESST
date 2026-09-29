function figs = plot_accessibility_and_fixed_target_uv(geom, access, uv)
%PLOT_ACCESSIBILITY_AND_FIXED_TARGET_UV Make diagnostic plots.
%
% This function keeps the user's global plot defaults, but creates larger
% figures, short titles, readable legends, and scaled UV units.

    if isfield(geom, 'dim')
        xTime = geom.dim.timePlot;
        xLabel = geom.dim.timePlotLabel;
        L12 = geom.dim.baselines.L12_km;
        L13 = geom.dim.baselines.L13_km;
        L23 = geom.dim.baselines.L23_km;
        offsetTotal = geom.dim.centroidOffset_km;
        inPlaneOffset = geom.dim.inPlaneOffset_km;
        absPlaneOffset = geom.dim.absPlaneOffset_km;
        areaVals = geom.dim.area_km2;
        baselineLabel = 'Baseline length [km]';
        offsetLabel = 'Offset [km]';
        areaLabel = 'Area [km$^2$]';
    else
        xTime = geom.tau;
        xLabel = '$t/T_C$';
        L12 = geom.baselines.L12;
        L13 = geom.baselines.L13;
        L23 = geom.baselines.L23;
        offsetTotal = geom.centroidOffset;
        inPlaneOffset = geom.inPlaneOffset;
        absPlaneOffset = geom.absPlaneOffset;
        areaVals = geom.area;
        baselineLabel = 'Baseline length [nondim.]';
        offsetLabel = 'Offset [nondim.]';
        areaLabel = 'Area [nondim.$^2$]';
    end

    r = geom.r;
    rC = geom.rC;
    acc = access.accessible;

    [uvScale, uvUnitLabel] = choose_uv_plot_scale(uv);
    uv12 = uv.uvw12(:, 1:2)/uvScale;
    uv13 = uv.uvw13(:, 1:2)/uvScale;
    uv23 = uv.uvw23(:, 1:2)/uvScale;
    rho12 = uv.rho12/uvScale;
    rho13 = uv.rho13/uvScale;
    rho23 = uv.rho23/uvScale;
    w12 = uv.w12/uvScale;
    w13 = uv.w13/uvScale;
    w23 = uv.w23/uvScale;

    figs = gobjects(0);

    % 3D trajectory geometry in synodic frame.
    figs(end+1) = figure('Name', 'QPO collectors and periodic combiner', ...
        'Color', 'w', 'Position', [80, 80, 1250, 800]);
    hold on; grid on; axis equal;
    plot3(rC(:,1), rC(:,2), rC(:,3), 'k-', 'LineWidth', 1.8, ...
        'DisplayName', 'Combiner PO');
    for i = 1:3
        plot3(r(:,1,i), r(:,2,i), r(:,3,i), 'LineWidth', 1.4, ...
            'DisplayName', sprintf('Collector %d', i));
    end
    plot_triangle_snapshot(r, 1, 'Initial triangle');
    plot_triangle_snapshot(r, size(r,1), 'Final triangle');
    xlabel('$x_S$'); ylabel('$y_S$'); zlabel('$z_S$');
    title('Combiner and QPO collectors');
    legend('Location', 'eastoutside');
    view(45, 25);

    % Baselines.
    figs(end+1) = figure('Name', 'Collector-collector baselines', ...
        'Color', 'w', 'Position', [100, 100, 1200, 650]);
    hold on; grid on;
    plot(xTime, L12, 'DisplayName', '$\|b_{12}\|$');
    plot(xTime, L13, 'DisplayName', '$\|b_{13}\|$');
    plot(xTime, L23, 'DisplayName', '$\|b_{23}\|$');
    xlabel(xLabel);
    ylabel(baselineLabel);
    title('Collector--collector baselines');
    legend('Location', 'eastoutside');

    % Combiner-centroid offset.
    figs(end+1) = figure('Name', 'Combiner-centroid offset', ...
        'Color', 'w', 'Position', [100, 100, 1200, 800]);
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile; hold on; grid on;
    plot(xTime, offsetTotal, 'DisplayName', '$\|\bar r-r_C\|$');
    plot(xTime, inPlaneOffset, 'DisplayName', 'In-plane');
    plot(xTime, absPlaneOffset, 'DisplayName', 'Plane-normal');
    xlabel(xLabel);
    ylabel(offsetLabel);
    title('Combiner offset');
    legend('Location', 'eastoutside');
    nexttile; hold on; grid on;
    plot(xTime, geom.normalizedCentroidOffset, 'DisplayName', ...
        '$\|\bar r-r_C\|/R_{\rm rms}$');
    plot(xTime, geom.normalizedInPlaneOffset, 'DisplayName', ...
        'In-plane$/R_{\rm rms}$');
    plot(xTime, geom.normalizedAbsPlaneOffset, 'DisplayName', ...
        'Plane-normal$/R_{\rm rms}$');
    xlabel(xLabel);
    ylabel('Normalized offset');
    title('Normalized combiner offset');
    legend('Location', 'eastoutside');

    % Triangle area and conditioning.
    figs(end+1) = figure('Name', 'Triangle area and conditioning', ...
        'Color', 'w', 'Position', [100, 100, 1200, 800]);
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile; hold on; grid on;
    plot(xTime, areaVals, 'DisplayName', '$A_\triangle$');
    xlabel(xLabel);
    ylabel(areaLabel);
    title('Collector triangle area');
    legend('Location', 'eastoutside');
    nexttile; hold on; grid on;
    plot(xTime, geom.sinAngle12_13, 'DisplayName', '$\sin\angle(b_{12},b_{13})$');
    xlabel(xLabel);
    ylabel('Conditioning');
    title('Normal conditioning');
    ylim([0, 1.05]);
    legend('Location', 'eastoutside');

    % Natural normal in synodic frame.
    figs(end+1) = figure('Name', 'Natural normal path synodic', ...
        'Color', 'w', 'Position', [80, 80, 1100, 800]);
    hold on; grid on; axis equal;
    [xs, ys, zs] = sphere(40);
    surf(xs, ys, zs, 'FaceAlpha', 0.06, 'EdgeAlpha', 0.08, ...
        'HandleVisibility', 'off');
    plot3(geom.normalSynodic(:,1), geom.normalSynodic(:,2), geom.normalSynodic(:,3), ...
        'LineWidth', 1.8, 'DisplayName', '$\hat n_S(t)$');
    plot3(access.sTargetSynodic(:,1), access.sTargetSynodic(:,2), access.sTargetSynodic(:,3), ...
        '--', 'LineWidth', 1.5, 'DisplayName', '$\hat s_S(t)$');
    xlabel('$x_S$'); ylabel('$y_S$'); zlabel('$z_S$');
    title('Natural normal, synodic frame');
    legend('Location', 'eastoutside');
    view(45, 25);

    % Natural normal in inertial frame.
    figs(end+1) = figure('Name', 'Natural normal path inertial', ...
        'Color', 'w', 'Position', [80, 80, 1100, 800]);
    hold on; grid on; axis equal;
    [xi, yi, zi] = sphere(40);
    surf(xi, yi, zi, 'FaceAlpha', 0.06, 'EdgeAlpha', 0.08, ...
        'HandleVisibility', 'off');
    plot3(geom.normalInertial(:,1), geom.normalInertial(:,2), geom.normalInertial(:,3), ...
        'LineWidth', 1.8, 'DisplayName', '$\hat n_I(t)$');
    scatter3(access.sTargetInertial(:,1), access.sTargetInertial(:,2), access.sTargetInertial(:,3), ...
        '*', 'LineWidth', 1.5, 'DisplayName', '$\hat s_I(t)$');
    xlabel('$x_I$'); ylabel('$y_I$'); zlabel('$z_I$');
    title('Natural normal, inertial frame');
    legend('Location', 'eastoutside');
    view(45, 25);

    % Off-normal/accessibility.
    figs(end+1) = figure('Name', 'Target accessibility angle', ...
        'Color', 'w', 'Position', [100, 100, 1200, 650]);
    hold on; grid on;
    plot(xTime, access.thetaOffDeg, 'DisplayName', '$\theta_{\rm off}(t)$');
    yline(access.thetaMaxDeg, '--', 'DisplayName', '$\theta_{\max}$');
    scatter(xTime(acc), access.thetaOffDeg(acc), 12, 'filled', ...
        'DisplayName', 'Accessible samples');
    xlabel(xLabel);
    ylabel('Off-normal angle [deg]');
    title('Fixed-target accessibility');
    legend('Location', 'eastoutside');
    ylim([0, 90]);

        % Fixed target UV tracks.
    figs(end+1) = figure('Name', 'Fixed-target UV coverage', ...
        'Color', 'w', 'Position', [100, 100, 950, 850]);

    ax = axes('Parent', figs(end));
    hold(ax, 'on');
    grid(ax, 'on');
    axis(ax, 'equal');

    baselineColors = lines(3);

    hUv = gobjects(3, 1);

    hUv(1) = plot_uv_track_with_accessibility( ...
        ax, uv12, acc, '$\mathbf{B}_{12}$', baselineColors(1,:));

    hUv(2) = plot_uv_track_with_accessibility( ...
        ax, uv13, acc, '$\mathbf{B}_{13}$', baselineColors(2,:));

    hUv(3) = plot_uv_track_with_accessibility( ...
        ax, uv23, acc, '$\mathbf{B}_{23}$', baselineColors(3,:));

    xlabel(ax, sprintf('$u$ [%s]', uvUnitLabel));
    ylabel(ax, sprintf('$v$ [%s]', uvUnitLabel));
    title(ax, 'Fixed-target UV coverage');

    legend(ax, hUv, ...
        {'$\mathbf{B}_{12}$', '$\mathbf{B}_{13}$', '$\mathbf{B}_{23}$'}, ...
        'Location', 'eastoutside');

    set_symmetric_uv_limits(ax, [uv12; uv13; uv23; -uv12; -uv13; -uv23]);

    % UV radius and w component.
    figs(end+1) = figure('Name', 'Fixed-target uvw versus time', ...
        'Color', 'w', 'Position', [100, 100, 1200, 800]);
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile; hold on; grid on;
    plot(xTime, rho12, 'DisplayName', '$\rho_{uv,12}$');
    plot(xTime, rho13, 'DisplayName', '$\rho_{uv,13}$');
    plot(xTime, rho23, 'DisplayName', '$\rho_{uv,23}$');
    xlabel(xLabel); ylabel(sprintf('$\rho_{uv}$ [%s]', uvUnitLabel));
    title('Projected baseline radius');
    legend('Location', 'eastoutside');
    nexttile; hold on; grid on;
    plot(xTime, w12, 'DisplayName', '$w_{12}$');
    plot(xTime, w13, 'DisplayName', '$w_{13}$');
    plot(xTime, w23, 'DisplayName', '$w_{23}$');
    xlabel(xLabel); ylabel(sprintf('$w$ [%s]', uvUnitLabel));
    title('Line-of-sight baseline');
    legend('Location', 'eastoutside');
end

function [scale, label] = choose_uv_plot_scale(uv)
%CHOOSE_UV_PLOT_SCALE Pick a readable unit for UV plots.

    vals = [uv.uvw12(:); uv.uvw13(:); uv.uvw23(:)];
    maxAbs = max(abs(vals), [], 'omitnan');

    if isempty(maxAbs) || ~isfinite(maxAbs) || maxAbs == 0
        scale = 1.0;
        label = 'wavelengths';
    elseif maxAbs >= 1e9
        scale = 1e9;
        label = 'G$\lambda$';
    elseif maxAbs >= 1e6
        scale = 1e6;
        label = 'M$\lambda$';
    elseif maxAbs >= 1e3
        scale = 1e3;
        label = 'k$\lambda$';
    else
        scale = 1.0;
        label = 'wavelengths';
    end
end

function set_symmetric_uv_limits(ax, uv2)
%SET_SYMMETRIC_UV_LIMITS Use symmetric limits for an equal-aspect UV plot.

    lim = max(abs(uv2(:)), [], 'omitnan');
    if isempty(lim) || ~isfinite(lim) || lim == 0
        lim = 1;
    end
    lim = 1.05*lim;
    xlim(ax, [-lim, lim]);
    ylim(ax, [-lim, lim]);
end


function hMain = plot_uv_track_with_accessibility(ax, uv2, acc, labelText, colorVal)
%PLOT_UV_TRACK_WITH_ACCESSIBILITY Plot one baseline and its Hermitian copy.
%
% solid line  : sampled baseline track (u,v)
% dashed line : Hermitian counterpart (-u,-v)
% markers     : accessible samples
%
% Each physical baseline keeps one color.

    acc = acc(:);

    hMain = plot(ax, uv2(:,1), uv2(:,2), '-', ...
        'Color', colorVal, ...
        'LineWidth', 1.6, ...
        'DisplayName', labelText);

    plot(ax, -uv2(:,1), -uv2(:,2), '--', ...
        'Color', colorVal, ...
        'LineWidth', 1.2, ...
        'HandleVisibility', 'off');

    if any(acc)
        scatter(ax, uv2(acc,1), uv2(acc,2), 18, ...
            'MarkerFaceColor', colorVal, ...
            'MarkerEdgeColor', colorVal, ...
            'HandleVisibility', 'off');

        scatter(ax, -uv2(acc,1), -uv2(acc,2), 18, ...
            'MarkerFaceColor', colorVal, ...
            'MarkerEdgeColor', colorVal, ...
            'HandleVisibility', 'off');
    end
end