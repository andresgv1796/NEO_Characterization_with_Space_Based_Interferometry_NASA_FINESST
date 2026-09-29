function save_all_figures(figs, figDir, filePrefix)
%SAVE_ALL_FIGURES Save figure handles contained in a struct or array.
%
% save_all_figures(figs, figDir)
% save_all_figures(figs, figDir, filePrefix)

    if nargin < 3 || isempty(filePrefix)
        filePrefix = '';
    end

    if ~exist(figDir, 'dir')
        mkdir(figDir);
    end

    if isstruct(figs)
        names = fieldnames(figs);

        for k = 1:numel(names)
            figName = names{k};
            hFig = figs.(figName);

            if isempty(hFig) || ~isvalid(hFig)
                continue;
            end

            save_one_figure(hFig, figDir, filePrefix, figName);
        end

    elseif isa(figs, 'matlab.ui.Figure')
        for k = 1:numel(figs)
            figName = sprintf('figure_%02d', k);
            save_one_figure(figs(k), figDir, filePrefix, figName);
        end

    else
        error('figs must be a struct of figure handles or an array of figure handles.');
    end
end

function save_one_figure(hFig, figDir, filePrefix, figName)

    figName = regexprep(figName, '[^a-zA-Z0-9_+-]', '_');

    if isempty(filePrefix)
        baseName = figName;
    else
        baseName = sprintf('%s_%s', filePrefix, figName);
    end

    pngFile = fullfile(figDir, [baseName, '.png']);
    figFile = fullfile(figDir, [baseName, '.fig']);

    drawnow;

    exportgraphics(hFig, pngFile, 'Resolution', 300);
    savefig(hFig, figFile);

    fprintf('Saved figure: %s\n', pngFile);
end