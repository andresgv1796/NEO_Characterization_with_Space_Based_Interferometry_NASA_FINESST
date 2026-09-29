function figs = plot_target_suitability_scan(scan)
%PLOT_TARGET_SUITABILITY_SCAN Plot sampled-target suitability metrics.
%
% figs = plot_target_suitability_scan(scan)
%
% Required scan fields:
%   raDeg, decDeg, score, accessibleFraction, uvAnisotropy
%
% Optional scan fields:
%   thetaMinDeg, uvMax_Glambda, bestIndex, maxOffNormalDeg

    requiredFields = {'raDeg', 'decDeg', 'score', ...
        'accessibleFraction', 'uvAnisotropy'};

    for k = 1:numel(requiredFields)
        if ~isfield(scan, requiredFields{k})
            error('scan.%s is required by plot_target_suitability_scan.', ...
                requiredFields{k});
        end
    end

    raDeg = scan.raDeg(:);
    decDeg = scan.decDeg(:);
    score = scan.score(:);
    fAcc = scan.accessibleFraction(:);
    etaUv = scan.uvAnisotropy(:);

    if isfield(scan, 'thetaMinDeg')
        thetaMinDeg = scan.thetaMinDeg(:);
    else
        thetaMinDeg = NaN(size(score));
    end

    if isfield(scan, 'uvMax_Glambda')
        uvMax_Glambda = scan.uvMax_Glambda(:);
    else
        uvMax_Glambda = NaN(size(score));
    end

    if isfield(scan, 'bestIndex') && ~isempty(scan.bestIndex)
        bestIdx = scan.bestIndex;
    else
        [~, bestIdx] = max(score);
    end

    figs = struct();

    figs.scoreMap = figure('Name', 'target suitability score', ...
        'Color', 'w', 'Position', [100, 100, 1150, 650]);

    scatter(raDeg, decDeg, 34, score, 'filled');
    hold on;
    plot(raDeg(bestIdx), decDeg(bestIdx), 'kp', ...
        'MarkerSize', 16, 'MarkerFaceColor', 'y');
    hold off;

    grid on;
    box on;
    xlabel('RA-like angle [deg]', 'Interpreter', 'latex');
    ylabel('Dec-like angle [deg]', 'Interpreter', 'latex');
    title('Target suitability score', 'Interpreter', 'latex');
    cb = colorbar;
    cb.Label.String = '$f_{\rm acc}\eta_{uv}$';
    cb.Label.Interpreter = 'latex';
    xlim([0, 360]);
    ylim([-90, 90]);

    figs.accessibilityMap = figure('Name', 'target accessibility fraction', ...
        'Color', 'w', 'Position', [130, 130, 1150, 650]);

    scatter(raDeg, decDeg, 34, fAcc, 'filled');
    hold on;
    plot(raDeg(bestIdx), decDeg(bestIdx), 'kp', ...
        'MarkerSize', 16, 'MarkerFaceColor', 'y');
    hold off;

    grid on;
    box on;
    xlabel('RA-like angle [deg]', 'Interpreter', 'latex');
    ylabel('Dec-like angle [deg]', 'Interpreter', 'latex');
    title('Accessibility fraction', 'Interpreter', 'latex');
    cb = colorbar;
    cb.Label.String = '$f_{\rm acc}$';
    cb.Label.Interpreter = 'latex';
    xlim([0, 360]);
    ylim([-90, 90]);

    figs.uvAnisotropyMap = figure('Name', 'target uv anisotropy', ...
        'Color', 'w', 'Position', [160, 160, 1150, 650]);

    scatter(raDeg, decDeg, 34, etaUv, 'filled');
    hold on;
    plot(raDeg(bestIdx), decDeg(bestIdx), 'kp', ...
        'MarkerSize', 16, 'MarkerFaceColor', 'y');
    hold off;

    grid on;
    box on;
    xlabel('RA-like angle [deg]', 'Interpreter', 'latex');
    ylabel('Dec-like angle [deg]', 'Interpreter', 'latex');
    title('UV anisotropy', 'Interpreter', 'latex');
    cb = colorbar;
    cb.Label.String = '$\eta_{uv}$';
    cb.Label.Interpreter = 'latex';
    xlim([0, 360]);
    ylim([-90, 90]);

    figs.scoreTradeoff = figure('Name', 'target suitability tradeoff', ...
        'Color', 'w', 'Position', [190, 190, 950, 760]);

    scatter(fAcc, etaUv, 44, score, 'filled');
    hold on;
    plot(fAcc(bestIdx), etaUv(bestIdx), 'kp', ...
        'MarkerSize', 16, 'MarkerFaceColor', 'y');
    hold off;

    grid on;
    box on;
    xlabel('$f_{\rm acc}$', 'Interpreter', 'latex');
    ylabel('$\eta_{uv}$', 'Interpreter', 'latex');
    title('Accessibility versus UV anisotropy', 'Interpreter', 'latex');
    cb = colorbar;
    cb.Label.String = '$f_{\rm acc}\eta_{uv}$';
    cb.Label.Interpreter = 'latex';
    xlim([0, 1]);
    ylim([0, 1]);

    figs.summary = figure('Name', 'target suitability summary', ...
        'Color', 'w', 'Position', [220, 220, 1050, 720]);

    [sortedScore, idxSort] = sort(score, 'descend');
    nShow = min(30, numel(sortedScore));
    topIdx = idxSort(1:nShow);

    bar(1:nShow, sortedScore(1:nShow));
    grid on;
    box on;
    xlabel('Target rank', 'Interpreter', 'latex');
    ylabel('$f_{\rm acc}\eta_{uv}$', 'Interpreter', 'latex');
    title('Top sampled targets', 'Interpreter', 'latex');

    labels = compose('%.0f/%.0f', raDeg(topIdx), decDeg(topIdx));
    xticks(1:nShow);
    xticklabels(labels);
    xtickangle(45);

    fprintf('\nTop sampled targets:\n');
    fprintf('  rank      RA       Dec      score      f_acc     eta_uv   rho_max[Glam]  theta_min[deg]\n');

    nPrint = min(10, numel(topIdx));
    for k = 1:nPrint
        ii = topIdx(k);
        fprintf('  %4d  %8.3f  %8.3f  %8.4f  %8.4f  %8.4f  %12.4g  %12.4f\n', ...
            k, raDeg(ii), decDeg(ii), score(ii), fAcc(ii), etaUv(ii), ...
            uvMax_Glambda(ii), thetaMinDeg(ii));
    end
end
