function plot_uv_track_with_accessibility(uv2, accessible, displayName)
%PLOT_UV_TRACK_WITH_ACCESSIBILITY Plot UV track and highlight accessible samples.

    plot(uv2(:,1), uv2(:,2), '-', 'DisplayName', displayName);
    scatter(uv2(accessible,1), uv2(accessible,2), 10, 'filled', ...
        'HandleVisibility', 'off');
end
