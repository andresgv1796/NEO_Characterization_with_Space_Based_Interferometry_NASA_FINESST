clear;
clc;
close all;


%% ============================================================
%  ITOKAWA ROTATION ANIMATION DEMO
%
%  Uses:
%    - Gaskell Q=64 Itokawa shape model
%    - native Gaskell body-fixed frame
%    - +z_B as rotation pole
%    - constant principal-axis rotation
%  ============================================================


%% ============================================================
%  1. Load Gaskell Itokawa shape model
%  ============================================================

itokawaFile = fullfile( ...
    'ShapeModels', ...
    'Itokawa', ...
    'ver64q.tab');


if ~isfile(itokawaFile)

    error( ...
        ['Could not find the Gaskell Itokawa model:\n\n%s\n\n', ...
         'Current MATLAB folder:\n%s'], ...
        itokawaFile, ...
        pwd);

end


shapeItokawa = ...
    load_gaskell_itokawa_shape( ...
        itokawaFile);


fprintf('\n');
fprintf('============================================================\n');
fprintf('GASKELL ITOKAWA MODEL\n');
fprintf('============================================================\n');
fprintf('Vertices : %d\n',size(shapeItokawa.vertices,1));
fprintf('Facets   : %d\n',size(shapeItokawa.faces,1));
fprintf('============================================================\n');


%% ============================================================
%  2. Define rotation model
%  ============================================================
%
%  Gaskell body-frame convention:
%
%       +z_B = rotation pole
%       +x_B = zero longitude
%
%  r_I = R_IB(t) r_B
%
%  For this demonstration:
%
%       R_IB(t0) = I
%
%  so initially the body frame and inertial plotting frame coincide.
%
%  This is NOT intended to represent the absolute inertial attitude
%  of Itokawa at a particular UTC epoch.
%  ============================================================

rot = struct();


% Rotation period
rot.period_s = ...
    12.132381791 * 3600;


% Physical spin axis in Gaskell body coordinates
rot.spinAxis_B = ...
    [0;0;1];


% Initial body-to-inertial attitude
rot.R_IB0 = ...
    eye(3);


% Reference epoch
rot.t0_s = ...
    0.0;


% +1 = positive right-handed rotation about spinAxis_B
rot.spinSense = ...
    +1;


fprintf('\n');
fprintf('ROTATION MODEL\n');
fprintf('------------------------------------------------------------\n');
fprintf('Period       = %.9f h\n',rot.period_s/3600);
fprintf('Spin axis B  = [%+.3f %+.3f %+.3f]\n', ...
    rot.spinAxis_B);
fprintf('Spin sense   = %+d\n',rot.spinSense);
fprintf('------------------------------------------------------------\n');


%% ============================================================
%  3. Animation options
%  ============================================================

opts = struct();


% Number of frames over one full rotation
opts.nFrames = ...
    240;


% Playback / video frame rate
opts.fps = ...
    30;


% Leave empty if you do not want to save video
opts.savePath = ...
    'itokawa_rotation.mp4';


% Fixed inertial camera.
%
% Important:
% Do NOT look directly down the spin axis, otherwise the motion
% appears nearly two-dimensional.
opts.view = ...
    [42 24];


% Visualization aids
opts.showSpinAxis    = true;
opts.showOrigin      = true;
opts.showBodyAxes    = true;
opts.showInertialAxes = true;


% Surface appearance
opts.faceAlpha = ...
    1.0;


%% ============================================================
%  4. Animate
%  ============================================================

out = ...
    animate_rotating_asteroid( ...
        shapeItokawa, ...
        rot, ...
        opts);