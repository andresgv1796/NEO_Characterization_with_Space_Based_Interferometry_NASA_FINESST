function plot_triangle_snapshot(r, k, displayName)
%PLOT_TRIANGLE_SNAPSHOT Plot closed collector triangle at one time index.

    tri = squeeze(r(k, :, :)).';       % 3 x 3, rows are collectors
    tri = [tri; tri(1, :)];            % close the loop
    plot3(tri(:,1), tri(:,2), tri(:,3), '--', 'LineWidth', 1.3, ...
        'DisplayName', displayName);
end
