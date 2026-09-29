function demo_itokawa_rotation(savePath,itokawaFile)
%DEMO_ITOKAWA_ROTATION Animate the existing Gaskell Q=64 Itokawa model.
%
%   demo_itokawa_rotation
%   demo_itokawa_rotation('itokawa_rotation.mp4')
%   demo_itokawa_rotation('',fullfile(...,'ver64q.tab'))
%
% Requires the existing FINESST function:
%   load_gaskell_itokawa_shape.m
%
% The Gaskell body frame is retained. Its +z_B axis is the rotation pole.
% R_IB0 below is only the chosen t=0 reference attitude for this relative
% rotation demo; it is not an absolute UTC/SPICE attitude solution.

    if nargin < 1
        savePath = '';
    end

    if nargin < 2 || isempty(itokawaFile)
        itokawaFile = fullfile('ShapeModels','Itokawa','ver64q.tab');
    end

    if ~isfile(itokawaFile)
        error('Could not find Gaskell Itokawa model:\n%s',itokawaFile);
    end

    shape = load_gaskell_itokawa_shape(itokawaFile);

    % Gaskell/IAU model magnitude: 712.14376110 deg/day.
    % This corresponds to 12.132381791 h per rotation.
    rot.period_s   = 12.132381791 * 3600;
    rot.spinAxis_B = [0;0;1];
    rot.R_IB0      = eye(3);
    rot.t0_s       = 0;

    % Relative phase convention for the present FINESST experiment.
    % Keep this explicit because an absolute IAU/SPICE attitude layer can
    % later determine the physical sign and phase at a real UTC epoch.
    rot.spinSense  = +1;

    % Two cheap sanity checks replace a separate test package.
    R0 = asteroid_attitude_at_time(rot,rot.t0_s);
    RT = asteroid_attitude_at_time(rot,rot.t0_s + rot.period_s);

    fprintf('\nGaskell Itokawa rotating-target demo\n');
    fprintf('  vertices            = %d\n',size(shape.vertices,1));
    fprintf('  facets              = %d\n',size(shape.faces,1));
    fprintf('  rotation period     = %.9f h\n',rot.period_s/3600);
    fprintf('  one-period closure  = %.3e\n',norm(RT-R0,'fro'));
    fprintf('  det(R0)             = %.15f\n\n',det(R0));

    opts.nFrames = 180;
    opts.fps = 30;
    opts.savePath = savePath;
    opts.view = [35 20];
    opts.showSpinAxis = true;

    animate_rotating_asteroid(shape,rot,opts);
end
