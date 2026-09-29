function figsDirty = plot_dirty_image_results(simulated, dirty)
%PLOT_DIRTY_IMAGE_RESULTS
%
% Plot the simulated target, dirty beam, and dirty image.
%
% The simulated target is nonnegative and is shown with the same hot
% relative-intensity scale used for the simulated-target gallery.
%
% The dirty beam and dirty image are signed quantities, so they are shown
% with a centered diverging colormap.

    figsDirty = struct();

    figsDirty.main = figure( ...
        'Name', 'Dirty-image reconstruction', ...
        'Color', 'w', ...
        'Position', [100 100 1450 450]);

    tl = tiledlayout(1, 3, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    %% Simulated target

    ax1 = nexttile(tl);

    Ishow = simulated.I/max(simulated.I(:));

    imagesc(ax1, simulated.l_mas, simulated.m_mas, Ishow);
    axis(ax1, 'image');
    set(ax1, 'YDir', 'normal');

    colormap(ax1, hot(256));
    clim(ax1, [0 1]);

    format_dark_image_axes(ax1);

    cb1 = colorbar(ax1);
    cb1.Color = 'k';
    cb1.Label.String = 'Relative intensity';
    cb1.Label.Interpreter = 'latex';
    cb1.Label.Color = 'k';

    xlabel(ax1, '$l$ [mas]', 'Interpreter', 'latex', 'Color', 'k');
    ylabel(ax1, '$m$ [mas]', 'Interpreter', 'latex', 'Color', 'k');
    title(ax1, 'Simulated target', 'Interpreter', 'latex', 'Color', 'k');


    %% Dirty beam

    ax2 = nexttile(tl);

    beamShow = dirty.beam;
    beamLim = max(abs(beamShow(:)));

    imagesc(ax2, dirty.l_mas, dirty.m_mas, beamShow);
    axis(ax2, 'image');
    set(ax2, 'YDir', 'normal');

    colormap(ax2, red_white_blue_colormap(256));
    clim(ax2, [-beamLim beamLim]);

    format_dark_image_axes(ax2);

    cb2 = colorbar(ax2);
    cb2.Color = 'k';
    cb2.Label.String = 'Dirty beam';
    cb2.Label.Interpreter = 'latex';
    cb2.Label.Color = 'k';

    xlabel(ax2, '$l$ [mas]', 'Interpreter', 'latex', 'Color', 'k');
    ylabel(ax2, '$m$ [mas]', 'Interpreter', 'latex', 'Color', 'k');
    title(ax2, 'Dirty beam', 'Interpreter', 'latex', 'Color', 'k');


    %% Dirty image

    ax3 = nexttile(tl);

    imageShow = dirty.image;
    imageLim = max(abs(imageShow(:)));

    imagesc(ax3, dirty.l_mas, dirty.m_mas, imageShow);
    axis(ax3, 'image');
    set(ax3, 'YDir', 'normal');

    colormap(ax3, red_white_blue_colormap(256));
    clim(ax3, [-imageLim imageLim]);

    format_dark_image_axes(ax3);

    cb3 = colorbar(ax3);
    cb3.Color = 'k';
    cb3.Label.String = 'Dirty image';
    cb3.Label.Interpreter = 'latex';
    cb3.Label.Color = 'k';

    xlabel(ax3, '$l$ [mas]', 'Interpreter', 'latex', 'Color', 'k');
    ylabel(ax3, '$m$ [mas]', 'Interpreter', 'latex', 'Color', 'k');
    title(ax3, 'Dirty image', 'Interpreter', 'latex', 'Color', 'k');


    %% Main title

    title(tl, sprintf('Dirty imaging: %s', ...
        strrep(simulated.type, '_', ' ')), ...
        'Interpreter', 'latex', 'Color', 'k');
end


function format_dark_image_axes(ax)
%FORMAT_DARK_IMAGE_AXES Apply common image-axis formatting.

    ax.Color = 'k';
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.FontSize = 11;
    ax.LineWidth = 0.9;
    ax.Box = 'on';
end


function cmap = red_white_blue_colormap(n)
%RED_WHITE_BLUE_COLORMAP Simple centered diverging colormap.

    x = linspace(-1, 1, n).';

    r = ones(n,1);
    g = 1 - abs(x);
    b = ones(n,1);

    r(x < 0) = 1 + x(x < 0);
    b(x > 0) = 1 - x(x > 0);

    cmap = [r, g, b];
end