function gmos_plot_qpos_points(qpos, proj)
%GMOS_PLOT_QPOS_POINTS Plot a cloud of qpos points.
%
% qpos : 3xM
% proj : 'xy','xz','yz','xyz'

    if nargin < 2 || isempty(proj)
        proj = 'xy';
    end

    figure; hold on; grid on;
    if strcmpi(proj, 'xyz')
        plot3(qpos(1,:), qpos(2,:), qpos(3,:), '.');
        view(3);
        xlabel('x'); ylabel('y'); zlabel('z');
    else
        axis equal;
        switch lower(proj)
            case 'xy'
                plot(qpos(1,:), qpos(2,:), '.');
                xlabel('x'); ylabel('y');
            case 'xz'
                plot(qpos(1,:), qpos(3,:), '.');
                xlabel('x'); ylabel('z');
            case 'yz'
                plot(qpos(2,:), qpos(3,:), '.');
                xlabel('y'); ylabel('z');
            otherwise
                error('Unknown proj. Use xy, xz, yz, or xyz.');
        end
    end

    title('GMOS qpos: stroboscopic iterates');
end