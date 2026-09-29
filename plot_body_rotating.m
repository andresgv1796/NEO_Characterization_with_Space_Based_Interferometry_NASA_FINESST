function h = plot_body_rotating(ax,bodyName,center,opts)
%PLOT_BODY_ROTATING General textured/schematic body renderer.
%
%   h = plot_body_rotating(ax,bodyName,center, ...)
%
% Plots a body centered at `center` in the current coordinate system using a
% simple rotation about the z-axis, similarly to the older plot_Earth/plot_Moon/
% plot_Sun helpers.  The body radius is taken from a small internal catalog.
%
% Inputs
% ------
% ax                 target axes (use gca if empty)
% bodyName           e.g. "Sun", "Earth", "Moon", "Mars"
% center             [x y z] center in current coordinates
%
% Name-value options
% ------------------
% RotationDeg        z-axis rotation angle applied to the texture grid [deg]
% ScaleFactor        multiplicative display scale on the physical radius
% LengthUnit_km      length unit of the current coordinates [km]. For CR3BP,
%                    pass the system L* in km so the body size becomes
%                    radius_km / Lstar_km.
% Resolution         surface resolution
% TextureFile        custom texture file (optional)
% EquatorialRadius_km, PolarRadius_km   override internal radii
% FaceColor          fallback solid color when no texture file is found
%
% Notes
% -----
% The function tries to load a texture file if one is available either in the
% current folder or on the MATLAB path. Otherwise it falls back to a solid-color
% sphere/ellipsoid so the visualization still works.

arguments
    ax = []
    bodyName {mustBeTextScalar} = "Earth"
    center (1,3) double = [0 0 0]
    opts.RotationDeg (1,1) double = 0
    opts.ScaleFactor (1,1) double {mustBePositive} = 1
    opts.LengthUnit_km (1,1) double {mustBePositive} = 1
    opts.Resolution (1,1) double {mustBeInteger,mustBePositive} = 80
    opts.TextureFile string = ""
    opts.EquatorialRadius_km (1,1) double = NaN
    opts.PolarRadius_km (1,1) double = NaN
    opts.FaceColor (1,3) double = [0.65 0.65 0.65]
end

if isempty(ax)
    ax = gca;
end

body = lookup_body_defaults(bodyName);
if isfinite(opts.EquatorialRadius_km)
    body.req_km = opts.EquatorialRadius_km;
end
if isfinite(opts.PolarRadius_km)
    body.rpol_km = opts.PolarRadius_km;
end
if strlength(opts.TextureFile) > 0
    body.texture = opts.TextureFile;
end
if any(isnan([body.req_km body.rpol_km]))
    error('No radius data are available for body "%s".',string(bodyName));
end

req = opts.ScaleFactor * body.req_km / opts.LengthUnit_km;
rpol = opts.ScaleFactor * body.rpol_km / opts.LengthUnit_km;
[X,Y,Z] = ellipsoid(0,0,0,req,req,rpol,opts.Resolution);

C = [cosd(opts.RotationDeg), -sind(opts.RotationDeg); ...
     sind(opts.RotationDeg),  cosd(opts.RotationDeg)];
rotXY = C*[X(:).'; Y(:).'];
X = reshape(rotXY(1,:),size(X)) + center(1);
Y = reshape(rotXY(2,:),size(Y)) + center(2);
Z = -Z + center(3);

texPath = resolve_texture_file(body.texture);
if strlength(texPath) > 0
    cdata = imread(texPath);
    h = surf(ax,X,Y,Z, ...
        'FaceColor','texturemap', ...
        'CData',cdata, ...
        'EdgeColor','none');
else
    h = surf(ax,X,Y,Z, ...
        'FaceColor',body.color, ...
        'EdgeColor','none');
end

hide_handle_from_legend(h);
end


function body = lookup_body_defaults(bodyName)
name = lower(regexprep(char(string(bodyName)),'[^a-zA-Z0-9]',''));
body = struct('req_km',NaN,'rpol_km',NaN,'texture',"",'color',[0.65 0.65 0.65]);

switch name
    case 'sun'
        body.req_km = 696342;
        body.rpol_km = 695700;
        body.texture = "Sun_texturemap.jpg";
        body.color = [0.98 0.72 0.16];

    case 'earth'
        body.req_km = 6378.137;
        body.rpol_km = 6356.7521736;
        body.texture = "Earth_texturemap.jpg";
        body.color = [0.20 0.45 0.85];

    case 'moon'
        body.req_km = 1738.1;
        body.rpol_km = 1736.0;
        body.texture = "Moon_texturemap.jpg";
        body.color = [0.72 0.72 0.72];

    case 'mars'
        body.req_km = 3396.19;
        body.rpol_km = 3376.20;
        body.texture = "Mars_texturemap.jpg";
        body.color = [0.80 0.42 0.22];

    case 'venus'
        body.req_km = 6051.8;
        body.rpol_km = 6051.8;
        body.texture = "Venus_texturemap.jpg";
        body.color = [0.82 0.72 0.50];

    case 'mercury'
        body.req_km = 2439.7;
        body.rpol_km = 2439.7;
        body.texture = "Mercury_texturemap.jpg";
        body.color = [0.72 0.67 0.62];

    case 'jupiter'
        body.req_km = 71492;
        body.rpol_km = 66854;
        body.texture = "Jupiter_texturemap.jpg";
        body.color = [0.84 0.70 0.52];

    case 'saturn'
        body.req_km = 60268;
        body.rpol_km = 54364;
        body.texture = "Saturn_texturemap.jpg";
        body.color = [0.90 0.84 0.60];

    case 'uranus'
        body.req_km = 25559;
        body.rpol_km = 24973;
        body.texture = "Uranus_texturemap.jpg";
        body.color = [0.63 0.86 0.90];

    case 'neptune'
        body.req_km = 24764;
        body.rpol_km = 24341;
        body.texture = "Neptune_texturemap.jpg";
        body.color = [0.28 0.47 0.87];

    otherwise
        error('Unsupported bodyName "%s". Add it to lookup_body_defaults.',string(bodyName));
end
end


function texPath = resolve_texture_file(textureFile)
texPath = "";
if strlength(textureFile) == 0
    return
end
if isfile(textureFile)
    texPath = textureFile;
    return
end
w = which(char(textureFile));
if ~isempty(w)
    texPath = string(w);
end
end


function hide_handle_from_legend(h)
if isempty(h) || ~isgraphics(h), return, end
try
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
catch
    set(h,'HandleVisibility','off');
end
end
