function figsClean = plot_clean_image_results(simulated, dirty, clean)
%PLOT_CLEAN_IMAGE_RESULTS
%
% Plot the complete Hogbom CLEAN sequence.

    figsClean = struct();

    figsClean.main = figure( ...
        'Name', 'Hogbom CLEAN', ...
        'Color', 'w', ...
        'Position', [40 70 1700 850]);

    tl = tiledlayout(2, 4, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');


    %% Simulated target

    ax1 = nexttile(tl);

    targetDisplay = simulated.I/max(simulated.I(:));

    imagesc( ...
        ax1, ...
        simulated.l_mas, ...
        simulated.m_mas, ...
        targetDisplay);

    format_image_axis(ax1);

    colormap(ax1, hot(256));
    clim(ax1, [0 1]);

    cb1 = colorbar(ax1);
    cb1.Label.String = 'Relative intensity';
    cb1.Label.Interpreter = 'latex';

    title(ax1, 'Simulated target', 'Interpreter', 'latex');


    %% Dirty beam

    ax2 = nexttile(tl);

    imagesc( ...
        ax2, ...
        dirty.l_mas, ...
        dirty.m_mas, ...
        real(dirty.beam));

    format_image_axis(ax2);

    colormap(ax2, red_white_blue_colormap(256));
    clim(ax2, signed_clim(real(dirty.beam)));

    cb2 = colorbar(ax2);
    cb2.Label.String = 'Dirty beam';
    cb2.Label.Interpreter = 'latex';

    title(ax2, 'Dirty beam', 'Interpreter', 'latex');


    %% Fitted clean beam

    ax3 = nexttile(tl);

    imagesc( ...
        ax3, ...
        clean.l_mas, ...
        clean.m_mas, ...
        clean.cleanBeam);

    format_image_axis(ax3);

    colormap(ax3, hot(256));
    clim(ax3, [0 1]);

    cb3 = colorbar(ax3);
    cb3.Label.String = 'Clean beam';
    cb3.Label.Interpreter = 'latex';

    title(ax3, sprintf( ...
        ['Fitted clean beam\n', ...
         '$\\Theta_{\\mathrm{maj}}=%.3f$ mas, ', ...
         '$\\Theta_{\\mathrm{min}}=%.3f$ mas'], ...
        clean.cleanBeamFit.fwhmMajor_mas, ...
        clean.cleanBeamFit.fwhmMinor_mas), ...
        'Interpreter', 'latex');


    %% Dirty image

    ax4 = nexttile(tl);

    imagesc( ...
        ax4, ...
        dirty.l_mas, ...
        dirty.m_mas, ...
        real(dirty.image));

    format_image_axis(ax4);

    colormap(ax4, red_white_blue_colormap(256));
    clim(ax4, signed_clim(real(dirty.image)));

    cb4 = colorbar(ax4);
    cb4.Label.String = 'Dirty image';
    cb4.Label.Interpreter = 'latex';

    title(ax4, 'Initial dirty map', 'Interpreter', 'latex');


    %% Hogbom delta-component model

    ax5 = nexttile(tl);

    imagesc( ...
        ax5, ...
        clean.l_mas, ...
        clean.m_mas, ...
        clean.model);

    format_image_axis(ax5);

    colormap(ax5, red_white_blue_colormap(256));
    clim(ax5, signed_clim(clean.model));

    cb5 = colorbar(ax5);
    cb5.Label.String = 'Signed CLEAN components';
    cb5.Label.Interpreter = 'latex';

    title(ax5, 'Delta-component model', 'Interpreter', 'latex');


    %% Restored components

    ax6 = nexttile(tl);

    imagesc( ...
        ax6, ...
        clean.l_mas, ...
        clean.m_mas, ...
        clean.restoredComponents);

    format_image_axis(ax6);

    colormap(ax6, red_white_blue_colormap(256));
    clim(ax6, signed_clim(clean.restoredComponents));

    cb6 = colorbar(ax6);
    cb6.Label.String = 'Restored components';
    cb6.Label.Interpreter = 'latex';

    title(ax6, '$M_{\mathrm{CLEAN}}*B_C$', ...
        'Interpreter', 'latex');


    %% Final residual

    ax7 = nexttile(tl);

    imagesc( ...
        ax7, ...
        clean.l_mas, ...
        clean.m_mas, ...
        clean.residual);

    format_image_axis(ax7);

    colormap(ax7, red_white_blue_colormap(256));
    clim(ax7, signed_clim(clean.residual));

    cb7 = colorbar(ax7);
    cb7.Label.String = 'Residual';
    cb7.Label.Interpreter = 'latex';

    title(ax7, 'Final remaining map', 'Interpreter', 'latex');


    %% Final restored CLEAN image

    ax8 = nexttile(tl);

    imagesc( ...
        ax8, ...
        clean.l_mas, ...
        clean.m_mas, ...
        clean.image);

    format_image_axis(ax8);

    colormap(ax8, red_white_blue_colormap(256));
    clim(ax8, signed_clim(clean.image));

    cb8 = colorbar(ax8);
    cb8.Label.String = 'Restored CLEAN image';
    cb8.Label.Interpreter = 'latex';

    title(ax8, '$M_{\mathrm{CLEAN}}*B_C+R$', ...
        'Interpreter', 'latex');


    title(tl, sprintf( ...
        'Hogbom CLEAN: %s, %d iterations', ...
        strrep(simulated.type, '_', ' '), ...
        clean.nIter), ...
        'Interpreter', 'latex');


    %% Convergence history

    figsClean.history = figure( ...
        'Name', 'Hogbom CLEAN convergence', ...
        'Color', 'w', ...
        'Position', [100 100 1250 720]);

    tlHistory = tiledlayout(2, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');


    axH1 = nexttile(tlHistory);

    semilogy( ...
        axH1, ...
        clean.history.iteration, ...
        clean.history.absPeak, ...
        'o-', ...
        'LineWidth', 1.2, ...
        'MarkerSize', 4);

    hold(axH1, 'on');

    yline( ...
        axH1, ...
        clean.params.stopThreshold, ...
        '--', ...
        'Stopping threshold', ...
        'Interpreter', 'latex');

    xlabel(axH1, 'Iteration', 'Interpreter', 'latex');
    ylabel(axH1, '$\max |R|$', 'Interpreter', 'latex');
    title(axH1, 'Residual peak', 'Interpreter', 'latex');

    grid(axH1, 'on');
    box(axH1, 'on');


    axH2 = nexttile(tlHistory);

    semilogy( ...
        axH2, ...
        clean.history.iteration, ...
        clean.history.residualRms, ...
        'o-', ...
        'LineWidth', 1.2, ...
        'MarkerSize', 4);

    xlabel(axH2, 'Iteration', 'Interpreter', 'latex');
    ylabel(axH2, 'Residual RMS', 'Interpreter', 'latex');
    title(axH2, 'Residual RMS', 'Interpreter', 'latex');

    grid(axH2, 'on');
    box(axH2, 'on');


    axH3 = nexttile(tlHistory);

    plot( ...
        axH3, ...
        clean.history.iteration, ...
        clean.history.modelFlux, ...
        'o-', ...
        'LineWidth', 1.2, ...
        'MarkerSize', 4);

    xlabel(axH3, 'Iteration', 'Interpreter', 'latex');
    ylabel(axH3, '$\sum M_{\mathrm{CLEAN}}$', ...
        'Interpreter', 'latex');

    title(axH3, 'Accumulated component sum', ...
        'Interpreter', 'latex');

    grid(axH3, 'on');
    box(axH3, 'on');


    axH4 = nexttile(tlHistory);

    scatter( ...
        axH4, ...
        clean.history.l_mas, ...
        clean.history.m_mas, ...
        28, ...
        clean.history.iteration, ...
        'filled');

    axis(axH4, 'image');

    xlabel(axH4, '$l$ [mas]', 'Interpreter', 'latex');
    ylabel(axH4, '$m$ [mas]', 'Interpreter', 'latex');

    title(axH4, 'Selected component positions', ...
        'Interpreter', 'latex');

    cbH4 = colorbar(axH4);
    cbH4.Label.String = 'Iteration';
    cbH4.Label.Interpreter = 'latex';

    grid(axH4, 'on');
    box(axH4, 'on');
end


function format_image_axis(ax)

    axis(ax, 'image');
    set(ax, 'YDir', 'normal');

    xlabel(ax, '$l$ [mas]', 'Interpreter', 'latex');
    ylabel(ax, '$m$ [mas]', 'Interpreter', 'latex');

    box(ax, 'on');
end


function limits = signed_clim(arrayIn)

    maximumAbsoluteValue = max(abs(arrayIn(:)));

    if maximumAbsoluteValue == 0
        maximumAbsoluteValue = 1;
    end

    limits = [
        -maximumAbsoluteValue, ...
         maximumAbsoluteValue
    ];
end


function colormapOut = red_white_blue_colormap(nColors)

    coordinate = linspace(-1, 1, nColors).';

    red = ones(nColors,1);
    green = 1 - abs(coordinate);
    blue = ones(nColors,1);

    red(coordinate < 0) = ...
        1 + coordinate(coordinate < 0);

    blue(coordinate > 0) = ...
        1 - coordinate(coordinate > 0);

    colormapOut = [red, green, blue];
end