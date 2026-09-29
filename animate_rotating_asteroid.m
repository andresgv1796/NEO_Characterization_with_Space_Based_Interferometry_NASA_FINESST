function out = animate_rotating_asteroid(shape,rot,opts)
%ANIMATE_ROTATING_ASTEROID
% Animate a triangular asteroid mesh under a prescribed attitude R_IB(t).
%
% The body-fixed mesh is never modified:
%
%       r_I(t) = R_IB(t) r_B
%
% EXPECTED SHAPE FIELDS
%
%   shape.vertices       Nx3 body-frame vertices
%   shape.faces          Mx3 triangular facets
%
% EXPECTED ROTATION FIELDS
%
%   rot.period_s
%   rot.spinAxis_B
%   rot.R_IB0
%   rot.t0_s
%   rot.spinSense
%
% OPTIONAL opts FIELDS
%
%   nFrames           default 240
%   fps               default 30
%   savePath          default ''
%   view              default [42 24]
%
%   showSpinAxis      default true
%   showOrigin        default true
%   showBodyAxes      default true
%   showInertialAxes  default true
%
%   faceAlpha         default 1
%
% OUTPUT
%
%   out.figure
%   out.axes
%   out.patch
%   out.times_s


    %% ============================================================
    %  Options
    %  ============================================================

    if nargin < 3 || isempty(opts)

        opts = struct();

    end


    if ~isfield(opts,'nFrames')
        opts.nFrames = 240;
    end

    if ~isfield(opts,'fps')
        opts.fps = 30;
    end

    if ~isfield(opts,'savePath')
        opts.savePath = '';
    end

    if ~isfield(opts,'view')
        opts.view = [42 24];
    end

    if ~isfield(opts,'showSpinAxis')
        opts.showSpinAxis = true;
    end

    if ~isfield(opts,'showOrigin')
        opts.showOrigin = true;
    end

    if ~isfield(opts,'showBodyAxes')
        opts.showBodyAxes = true;
    end

    if ~isfield(opts,'showInertialAxes')
        opts.showInertialAxes = true;
    end

    if ~isfield(opts,'faceAlpha')
        opts.faceAlpha = 1.0;
    end


    %% ============================================================
    %  Mesh validation
    %  ============================================================

    if ~isfield(shape,'vertices')

        error( ...
            'shape.vertices is missing.');

    end


    if ~isfield(shape,'faces')

        error( ...
            'shape.faces is missing.');

    end


    V_B = ...
        shape.vertices;

    F = ...
        shape.faces;


    if size(V_B,2) ~= 3

        error( ...
            'shape.vertices must be Nx3.');

    end


    if size(F,2) ~= 3

        error( ...
            'shape.faces must be Mx3 triangular facets.');

    end


    %% ============================================================
    %  Time samples
    %
    %  Do not duplicate t=0 as the final movie frame.
    %  ============================================================

    fraction = ...
        (0:opts.nFrames-1).' / ...
        opts.nFrames;


    t = ...
        rot.t0_s + ...
        fraction*rot.period_s;


    %% ============================================================
    %  Initial attitude
    %  ============================================================

    [R_IB,phi] = ...
        asteroid_attitude_at_time( ...
            rot, ...
            t(1));


    V_I = ...
        (R_IB * V_B.').';


    %% ============================================================
    %  Characteristic physical scale
    %
    %  We do NOT recenter the mesh.
    %
    %  The Gaskell reference origin remains the rotation origin.
    %  ============================================================

    coordinateMax = ...
        max(abs(V_B(:)));


    if coordinateMax <= 0

        error( ...
            'Asteroid mesh has zero spatial extent.');

    end


    plotLimit = ...
        1.35*coordinateMax;


    axisLength = ...
        1.15*coordinateMax;


    bodyAxisLength = ...
        0.85*coordinateMax;


    inertialAxisLength = ...
        0.65*coordinateMax;


    %% ============================================================
    %  Figure
    %  ============================================================

    fig = figure( ...
        'Color','w', ...
        'Name','Itokawa rotation', ...
        'Position',[100 80 1050 850]);


    ax = ...
        axes(fig);


    hold(ax,'on');


    %% ============================================================
    %  Asteroid surface
    %  ============================================================

    hp = patch( ...
        ax, ...
        'Vertices',V_I, ...
        'Faces',F, ...
        'FaceColor',[0.70 0.70 0.70], ...
        'EdgeColor','none', ...
        'FaceAlpha',opts.faceAlpha, ...
        'FaceLighting','gouraud');


    %% ============================================================
    %  Physical 3-D plotting environment
    %  ============================================================

    xlim(ax, ...
        [-plotLimit plotLimit]);

    ylim(ax, ...
        [-plotLimit plotLimit]);

    zlim(ax, ...
        [-plotLimit plotLimit]);


    axis(ax,'manual');

    axis(ax,'vis3d');

    daspect(ax, ...
        [1 1 1]);

    pbaspect(ax, ...
        [1 1 1]);


    grid(ax,'on');

    box(ax,'on');


    view(ax, ...
        opts.view);

camproj(ax, ...
    'orthographic');


    camup(ax, ...
        [0 0 1]);


    xlabel(ax, ...
        '$x_I$ [km]', ...
        'Interpreter','latex');

    ylabel(ax, ...
        '$y_I$ [km]', ...
        'Interpreter','latex');

    zlabel(ax, ...
        '$z_I$ [km]', ...
        'Interpreter','latex');


    %% ============================================================
    %  Lighting
    %  ============================================================

    camlight(ax, ...
        'headlight');

    lighting(ax, ...
        'gouraud');

    material(ax, ...
        'dull');


    %% ============================================================
    %  Physical rotation origin
    %  ============================================================

    if opts.showOrigin

        hOrigin = ...
            plot3( ...
                ax, ...
                0,0,0, ...
                'ko', ...
                'MarkerFaceColor','k', ...
                'MarkerSize',6);

    else

        hOrigin = [];

    end


    %% ============================================================
    %  Fixed inertial coordinate triad
    %
    %  These NEVER rotate.
    %  ============================================================

    if opts.showInertialAxes

        hIX = plot3( ...
            ax, ...
            [0 inertialAxisLength], ...
            [0 0], ...
            [0 0], ...
            ':', ...
            'LineWidth',1.2);


        hIY = plot3( ...
            ax, ...
            [0 0], ...
            [0 inertialAxisLength], ...
            [0 0], ...
            ':', ...
            'LineWidth',1.2);


        hIZ = plot3( ...
            ax, ...
            [0 0], ...
            [0 0], ...
            [0 inertialAxisLength], ...
            ':', ...
            'LineWidth',1.2);


        text( ...
            ax, ...
            1.08*inertialAxisLength, ...
            0, ...
            0, ...
            '$x_I$', ...
            'Interpreter','latex');


        text( ...
            ax, ...
            0, ...
            1.08*inertialAxisLength, ...
            0, ...
            '$y_I$', ...
            'Interpreter','latex');


        text( ...
            ax, ...
            0, ...
            0, ...
            1.08*inertialAxisLength, ...
            '$z_I$', ...
            'Interpreter','latex');

    else

        hIX = [];
        hIY = [];
        hIZ = [];

    end


    %% ============================================================
    %  Spin axis
    %  ============================================================

    if opts.showSpinAxis

        sB = ...
            rot.spinAxis_B(:);

        sB = ...
            sB/norm(sB);


        sI = ...
            R_IB*sB;


        hSpin = ...
            plot3( ...
                ax, ...
                axisLength*[-sI(1) sI(1)], ...
                axisLength*[-sI(2) sI(2)], ...
                axisLength*[-sI(3) sI(3)], ...
                'k--', ...
                'LineWidth',1.8);

    else

        hSpin = [];

    end


    %% ============================================================
    %  Rotating BODY triad
    %
    %  These rotate with Itokawa.
    %  ============================================================

    eXB = ...
        [1;0;0];

    eYB = ...
        [0;1;0];

    eZB = ...
        [0;0;1];


    if opts.showBodyAxes

        eXI = ...
            R_IB*eXB;

        eYI = ...
            R_IB*eYB;

        eZI = ...
            R_IB*eZB;


        hBX = ...
            plot3( ...
                ax, ...
                [0 bodyAxisLength*eXI(1)], ...
                [0 bodyAxisLength*eXI(2)], ...
                [0 bodyAxisLength*eXI(3)], ...
                '-', ...
                'LineWidth',2.0);


        hBY = ...
            plot3( ...
                ax, ...
                [0 bodyAxisLength*eYI(1)], ...
                [0 bodyAxisLength*eYI(2)], ...
                [0 bodyAxisLength*eYI(3)], ...
                '-', ...
                'LineWidth',2.0);


        hBZ = ...
            plot3( ...
                ax, ...
                [0 bodyAxisLength*eZI(1)], ...
                [0 bodyAxisLength*eZI(2)], ...
                [0 bodyAxisLength*eZI(3)], ...
                '-', ...
                'LineWidth',2.0);


        hTX = ...
            text( ...
                ax, ...
                1.08*bodyAxisLength*eXI(1), ...
                1.08*bodyAxisLength*eXI(2), ...
                1.08*bodyAxisLength*eXI(3), ...
                '$x_B$', ...
                'Interpreter','latex', ...
                'FontWeight','bold');


        hTY = ...
            text( ...
                ax, ...
                1.08*bodyAxisLength*eYI(1), ...
                1.08*bodyAxisLength*eYI(2), ...
                1.08*bodyAxisLength*eYI(3), ...
                '$y_B$', ...
                'Interpreter','latex', ...
                'FontWeight','bold');


        hTZ = ...
            text( ...
                ax, ...
                1.08*bodyAxisLength*eZI(1), ...
                1.08*bodyAxisLength*eZI(2), ...
                1.08*bodyAxisLength*eZI(3), ...
                '$z_B$', ...
                'Interpreter','latex', ...
                'FontWeight','bold');

    else

        hBX = [];
        hBY = [];
        hBZ = [];

        hTX = [];
        hTY = [];
        hTZ = [];

    end


    %% ============================================================
    %  Title
    %  ============================================================

    hTitle = ...
        title( ...
            ax, ...
            sprintf( ...
                'Itokawa:  t/T = %.3f, phase = %.1f deg', ...
                0, ...
                rad2deg(phi)), ...
            'Interpreter','none');


    %% ============================================================
    %  Video
    %  ============================================================

    writer = [];


    if ~isempty(opts.savePath)

        writer = ...
            VideoWriter( ...
                opts.savePath, ...
                'MPEG-4');


        writer.FrameRate = ...
            opts.fps;


        open(writer);

    end


    cleanup = ...
        onCleanup( ...
            @() local_close_writer(writer));  


    %% ============================================================
    %  Animation loop
    %  ============================================================

    for k = 1:numel(t)


        %% --------------------------------------------------------
        %  Current attitude
        %  --------------------------------------------------------

        [R_IB,phi] = ...
            asteroid_attitude_at_time( ...
                rot, ...
                t(k));


        %% --------------------------------------------------------
        %  Rotate asteroid
        %  --------------------------------------------------------

        V_I = ...
            (R_IB*V_B.').';


        set( ...
            hp, ...
            'Vertices',V_I);


        %% --------------------------------------------------------
        %  Spin axis
        %  --------------------------------------------------------

        if ~isempty(hSpin)

            sI = ...
                R_IB*sB;


            set( ...
                hSpin, ...
                'XData',axisLength*[-sI(1) sI(1)], ...
                'YData',axisLength*[-sI(2) sI(2)], ...
                'ZData',axisLength*[-sI(3) sI(3)]);

        end


        %% --------------------------------------------------------
        %  Rotating body axes
        %  --------------------------------------------------------

        if opts.showBodyAxes

            eXI = ...
                R_IB*eXB;

            eYI = ...
                R_IB*eYB;

            eZI = ...
                R_IB*eZB;


            set( ...
                hBX, ...
                'XData',[0 bodyAxisLength*eXI(1)], ...
                'YData',[0 bodyAxisLength*eXI(2)], ...
                'ZData',[0 bodyAxisLength*eXI(3)]);


            set( ...
                hBY, ...
                'XData',[0 bodyAxisLength*eYI(1)], ...
                'YData',[0 bodyAxisLength*eYI(2)], ...
                'ZData',[0 bodyAxisLength*eYI(3)]);


            set( ...
                hBZ, ...
                'XData',[0 bodyAxisLength*eZI(1)], ...
                'YData',[0 bodyAxisLength*eZI(2)], ...
                'ZData',[0 bodyAxisLength*eZI(3)]);


            set( ...
                hTX, ...
                'Position', ...
                1.08*bodyAxisLength*eXI.');


            set( ...
                hTY, ...
                'Position', ...
                1.08*bodyAxisLength*eYI.');


            set( ...
                hTZ, ...
                'Position', ...
                1.08*bodyAxisLength*eZI.');

        end


        %% --------------------------------------------------------
        %  Title
        %  --------------------------------------------------------

        set( ...
            hTitle, ...
            'String', ...
            sprintf( ...
                'Itokawa:  t/T = %.3f, phase = %.1f deg', ...
                fraction(k), ...
                rad2deg(phi)));


        %% --------------------------------------------------------
        %  Render
        %  --------------------------------------------------------

        drawnow;


        %% --------------------------------------------------------
        %  Save / playback
        %  --------------------------------------------------------

        if ~isempty(writer)

            writeVideo( ...
                writer, ...
                getframe(fig));

        elseif k < numel(t)

            pause( ...
                1/opts.fps);

        end

    end


    %% ============================================================
    %  Output
    %  ============================================================

    out = struct();

    out.figure = ...
        fig;

    out.axes = ...
        ax;

    out.patch = ...
        hp;

    out.origin = ...
        hOrigin;

    out.spinAxis = ...
        hSpin;

    out.bodyAxes = ...
        [hBX hBY hBZ];

    out.inertialAxes = ...
        [hIX hIY hIZ];

    out.times_s = ...
        t;

end


function local_close_writer(writer)
%LOCAL_CLOSE_WRITER Close the video if one was opened.

    if ~isempty(writer)

        try

            close(writer);

        catch

        end

    end

end