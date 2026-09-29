function simulated = make_simulated_target_from_raster( ...
    img, range_km, targetType)
%MAKE_SIMULATED_TARGET_FROM_RASTER Convert physical raster to angular image.
%
%   simulated = make_simulated_target_from_raster( ...
%       img, range_km, targetType)
%
% Inputs
% ------
% img :
%       Output structure from rasterize_visible_facets.
%
% range_km :
%       Observer-to-target-center distance [km].
%
% targetType :
%       Descriptive target name used in plots and diagnostics.
%       Example:
%
%           'triaxial ellipsoid'
%           'irregular ellipsoid'
%           'contact binary'
%
% Output
% ------
% simulated :
%       Target structure compatible with the interferometric-observation
%       pipeline.
%
% Small-angle convention
% ----------------------
%
%       l = x_sky / D
%       m = y_sky / D
%
% where x_sky and y_sky are aligned with the interferometric +eU and +eV
% directions respectively.

arguments
    img struct
    range_km (1,1) double {mustBePositive}
    targetType char = 'synthetic asteroid'
end

%% Unit conversion

masPerRad = mas_per_rad();


%% Physical sky-plane coordinates -> angular coordinates

l_rad = img.x / range_km;
m_rad = img.y / range_km;

l_mas = l_rad * masPerRad;
m_mas = m_rad * masPerRad;

[L_rad, M_rad] = meshgrid(l_rad, m_rad);
[L_mas, M_mas] = meshgrid(l_mas, m_mas);


%% Normalize image for visibility calculations

I = img.I;

imageSum = sum(I(:));

if imageSum <= 0
    error('Rasterized image has zero total intensity.');
end

I = I / imageSum;


%% Package simulated target

simulated.type = targetType;

simulated.I = I;

simulated.l_rad = l_rad;
simulated.m_rad = m_rad;

simulated.L_rad = L_rad;
simulated.M_rad = M_rad;

simulated.l_mas = l_mas;
simulated.m_mas = m_mas;

simulated.L_mas = L_mas;
simulated.M_mas = M_mas;


%% Physical metadata

simulated.range_km = range_km;

simulated.x_km = img.x;
simulated.y_km = img.y;

simulated.X_km = img.X;
simulated.Y_km = img.Y;

simulated.nPix = size(I,1);

% Angular half-field of view.
simulated.fov_mas = max(abs([l_mas(:); m_mas(:)]));


%% Preserve pre-normalization information

if isfield(img, 'Iraw')
    simulated.Iraw = img.Iraw;
end

if isfield(img, 'totalFlux')
    simulated.totalRelativeFlux = img.totalFlux;
end

if isfield(img, 'pixelArea')
    simulated.pixelArea_km2 = img.pixelArea;
end

end