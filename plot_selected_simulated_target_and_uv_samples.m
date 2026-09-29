function fig = plot_selected_simulated_target_and_uv_samples( ...
    simulated, ...
    uvFixed, ...
    access, ...
    baselineFields, ...
    baselineNames)
%PLOT_SELECTED_SIMULATED_TARGET_AND_UV_SAMPLES
%
% Plot the selected simulated target and the corresponding accessible
% physical UV samples.
%
% Image-grid sampling diagnostics are intentionally NOT displayed here.
% Those checks are performed separately by
%
%       validate_image_uv_sampling
%
% and printed using
%
%       print_image_uv_sampling_summary.
%
% This figure therefore displays only the physical quantities of interest:
%
%       1. simulated sky brightness I(l,m)
%       2. accessible interferometric UV coverage
%
% INPUTS
%   simulated      simulated target structure
%   uvFixed        fixed-target UVW structure
%   access         accessibility structure
%   baselineFields baseline field names in uvFixed
%   baselineNames  display names for each baseline
%
% OUTPUT
%   fig            figure handle


%% Validate baseline inputs

if isstring(baselineFields)
    baselineFields = cellstr(baselineFields);
end

if isstring(baselineNames)
    baselineNames = cellstr(baselineNames);
end

if numel(baselineFields) ~= numel(baselineNames)

    error([ ...
        'baselineFields and baselineNames must contain ', ...
        'the same number of entries.']);

end

nBaselines = numel(baselineFields);


%% Determine number of time samples

uvwFirst = local_extract_uvw( ...
    uvFixed, ...
    baselineFields{1});

nTime = size(uvwFirst,1);


%% Determine accessible time indices

accessibleIndex = local_accessible_index( ...
    access, ...
    nTime);

if isempty(accessibleIndex)

    error('No accessible UV samples are available to plot.');

end


%% Extract accessible UV samples

uCell = cell(nBaselines,1);
vCell = cell(nBaselines,1);

for b = 1:nBaselines

    uvw = local_extract_uvw( ...
        uvFixed, ...
        baselineFields{b});

    if size(uvw,1) ~= nTime

        error( ...
            'Baseline %s has an inconsistent number of time samples.', ...
            baselineFields{b});

    end

    uCell{b} = ...
        uvw(accessibleIndex,1)/1e9;

    vCell{b} = ...
        uvw(accessibleIndex,2)/1e9;

end


%% Image axes

masPerRad = mas_per_rad();

if isfield(simulated,'l_mas')

    l_mas = simulated.l_mas(:);

else

    l_mas = simulated.l_rad(:)*masPerRad;

end

if isfield(simulated,'m_mas')

    m_mas = simulated.m_mas(:);

else

    m_mas = simulated.m_rad(:)*masPerRad;

end


%% Normalize image for display

Ishow = real(simulated.I);

Imax = max(Ishow(:));

if ~isfinite(Imax) || Imax <= 0

    error('simulated.I contains no positive finite intensity.');

end

Ishow = Ishow/Imax;


%% Create figure

fig = figure( ...
    'Name', 'Selected simulated target and UV samples', ...
    'Color', 'w', ...
    'Position', [100 100 1320 570]);

tl = tiledlayout( ...
    fig, ...
    1, ...
    2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');


%% Simulated target

axImage = nexttile(tl,1);

imagesc( ...
    axImage, ...
    l_mas, ...
    m_mas, ...
    Ishow);

set(axImage, ...
    'YDir', 'normal');

axis(axImage,'image');

colormap(axImage,hot(256));
clim(axImage,[0 1]);

xlabel( ...
    axImage, ...
    '$l$ [mas]', ...
    'Interpreter','latex');

ylabel( ...
    axImage, ...
    '$m$ [mas]', ...
    'Interpreter','latex');

if isfield(simulated,'type')

    targetType = strrep( ...
        simulated.type, ...
        '_', ...
        ' ');

else

    targetType = 'simulated target';

end

title( ...
    axImage, ...
    sprintf('Simulated target: %s',targetType), ...
    'Interpreter','none');

cb = colorbar(axImage);

cb.Label.String = ...
    'Relative intensity';

axImage.Box = 'on';
axImage.Layer = 'top';


%% Accessible UV coverage

axUV = nexttile(tl,2);

hold(axUV,'on');
grid(axUV,'on');
box(axUV,'on');

baselineColors = lines(nBaselines);

for b = 1:nBaselines

    uPlot = uCell{b};
    vPlot = vCell{b};

    % Physical measured samples.
    plot( ...
        axUV, ...
        uPlot, ...
        vPlot, ...
        '.', ...
        'MarkerSize',8, ...
        'Color',baselineColors(b,:), ...
        'DisplayName',baselineNames{b});

    % Hermitian-conjugate locations.
    plot( ...
        axUV, ...
        -uPlot, ...
        -vPlot, ...
        'o', ...
        'MarkerSize',4, ...
        'Color',baselineColors(b,:), ...
        'HandleVisibility','off');

end
%% Symmetric UV plotting limits

uAll = vertcat(uCell{:});
vAll = vertcat(vCell{:});

uvExtent = max(abs([uAll; vAll]));

if ~isfinite(uvExtent) || uvExtent <= 0
    uvExtent = 1;
end

plotExtent = 1.10*uvExtent;

axis(axUV,'equal');

xlim( ...
    axUV, ...
    [-plotExtent plotExtent]);

ylim( ...
    axUV, ...
    [-plotExtent plotExtent]);


%% Overall title

title( ...
    tl, ...
    'Selected simulated target and accessible physical UV coverage', ...
    'Interpreter','latex');

end


%% ========================================================================
function uvw = local_extract_uvw(uvFixed, fieldName)
%LOCAL_EXTRACT_UVW Extract baseline UVW as an N x 3 array.


if ~isfield(uvFixed,fieldName)

    error('uvFixed.%s does not exist.',fieldName);

end

data = uvFixed.(fieldName);


if isnumeric(data)

    if size(data,2) == 3

        uvw = data;

    elseif size(data,1) == 3

        uvw = data.';

    else

        error( ...
            'uvFixed.%s must be N x 3 or 3 x N.', ...
            fieldName);

    end


elseif isstruct(data)

    if isfield(data,'u') && ...
       isfield(data,'v')

        u = data.u(:);
        v = data.v(:);

        if isfield(data,'w')

            w = data.w(:);

        else

            w = zeros(size(u));

        end

        if numel(u) ~= numel(v) || ...
           numel(u) ~= numel(w)

            error( ...
                'u, v, and w in uvFixed.%s must have equal lengths.', ...
                fieldName);

        end

        uvw = [u v w];


    elseif isfield(data,'uvw')

        uvw = data.uvw;

        if size(uvw,2) ~= 3 && ...
           size(uvw,1) == 3

            uvw = uvw.';

        end


    else

        error( ...
            'Cannot interpret uvFixed.%s as UVW data.', ...
            fieldName);

    end


else

    error( ...
        'Unsupported type for uvFixed.%s.', ...
        fieldName);

end


if size(uvw,2) ~= 3

    error( ...
        'Extracted uvFixed.%s must contain [u v w].', ...
        fieldName);

end

end


%% ========================================================================
function idx = local_accessible_index(access,nTime)
%LOCAL_ACCESSIBLE_INDEX Resolve accessible time indices.


if islogical(access) || ...
   (isnumeric(access) && isvector(access))

    values = access(:);

    if numel(values) == nTime

        idx = find(values ~= 0);

    else

        idx = values;

    end


elseif isstruct(access)

    idx = [];


    indexFields = { ...
        'accessibleIndex', ...
        'idxAccessible', ...
        'accessibleIdx', ...
        'indexAccessible'};

    for k = 1:numel(indexFields)

        fieldName = indexFields{k};

        if isfield(access,fieldName) && ...
           ~isempty(access.(fieldName))

            idx = access.(fieldName)(:);
            break

        end

    end


    if isempty(idx)

        maskFields = { ...
            'isAccessible', ...
            'accessible', ...
            'accessibleMask', ...
            'mask', ...
            'Iacc'};

        for k = 1:numel(maskFields)

            fieldName = maskFields{k};

            if isfield(access,fieldName)

                mask = access.(fieldName);

                if numel(mask) == nTime

                    idx = find(mask(:) ~= 0);
                    break

                end

            end

        end

    end


    if isempty(idx) && ...
       isfield(access,'thetaOffDeg') && ...
       isfield(access,'maxOffNormalDeg')

        thetaOffDeg = access.thetaOffDeg(:);

        if numel(thetaOffDeg) == nTime

            idx = find( ...
                thetaOffDeg <= ...
                access.maxOffNormalDeg);

        end

    end


    if isempty(idx)

        error( ...
            'Could not determine accessible samples from access.');

    end


else

    error('Unsupported accessibility input.');

end


idx = double(idx(:));


if any(~isfinite(idx)) || ...
   any(idx < 1) || ...
   any(idx > nTime) || ...
   any(abs(idx-round(idx)) > 0)

    error('Accessible indices are invalid.');

end


idx = unique(round(idx),'stable');

end