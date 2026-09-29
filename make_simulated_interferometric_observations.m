function obs = make_simulated_interferometric_observations( ...
    simulated, ...
    uvFixed, ...
    access, ...
    baselineFields, ...
    baselineNames, ...
    blockSize)
%MAKE_SIMULATED_INTERFEROMETRIC_OBSERVATIONS
%
% Generate ideal physical complex-visibility observations at all accessible
% baseline samples.
%
% Before evaluating the visibility model, this function checks whether the
% image grid supports the complete requested Cartesian UV extent using
%
%       |u| <= 1/(2*Delta_l)
%       |v| <= 1/(2*Delta_m).
%
% The resulting diagnostics are stored in
%
%       obs.gridCheck
%
% INPUTS
%   simulated      simulated target structure
%   uvFixed        fixed-target UVW structure
%   access         accessibility structure
%   baselineFields cell array of fields in uvFixed
%   baselineNames  display names for each baseline
%   blockSize      optional visibility evaluation block size
%
% OUTPUT
%   obs            physical interferometric observation structure


%% Defaults

if nargin < 6 || isempty(blockSize)

    blockSize = 64;

end


if ~isscalar(blockSize) || ...
   ~isfinite(blockSize) || ...
   blockSize < 1

    error('blockSize must be a positive finite scalar.');

end

blockSize = floor(blockSize);


%% Validate baseline inputs

if isstring(baselineFields)

    baselineFields = cellstr(baselineFields);

end


if isstring(baselineNames)

    baselineNames = cellstr(baselineNames);

end


if ~iscell(baselineFields) || isempty(baselineFields)

    error('baselineFields must be a nonempty cell array.');

end


if ~iscell(baselineNames)

    error('baselineNames must be a cell array.');

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


nTime = size(uvwFirst, 1);


if nTime < 1

    error('UV data contain no time samples.');

end


%% Accessibility indices

accessibleIndex = local_accessible_index( ...
    access, ...
    nTime);


if isempty(accessibleIndex)

    error('There are no accessible samples for this target.');

end


nAccessible = numel(accessibleIndex);


%% Extract all physical accessible UVW samples

uCell = cell(nBaselines, 1);
vCell = cell(nBaselines, 1);
wCell = cell(nBaselines, 1);

sampleCount = zeros(nBaselines, 1);


for b = 1:nBaselines

    fieldName = baselineFields{b};

    uvw = local_extract_uvw( ...
        uvFixed, ...
        fieldName);


    if size(uvw, 1) ~= nTime

        error([ ...
            'Baseline %s has %d time samples but expected %d.'], ...
            fieldName, ...
            size(uvw, 1), ...
            nTime);

    end


    uvwAccessible = ...
        uvw(accessibleIndex, :);


    uCell{b} = uvwAccessible(:, 1);
    vCell{b} = uvwAccessible(:, 2);
    wCell{b} = uvwAccessible(:, 3);

    sampleCount(b) = ...
        size(uvwAccessible, 1);

end


uAll = vertcat(uCell{:});
vAll = vertcat(vCell{:});
wAll = vertcat(wCell{:});


%% Mandatory image-grid support check

gridCheck = validate_image_uv_sampling( ...
    simulated, ...
    uAll, ...
    vAll);


%% Evaluate all physical complex visibilities

VAll = local_sample_visibility( ...
    simulated, ...
    uAll, ...
    vAll, ...
    blockSize);


VAll = VAll(:);


if numel(VAll) ~= numel(uAll)

    error([ ...
        'Visibility sampler returned %d values for %d UV samples.'], ...
        numel(VAll), ...
        numel(uAll));

end


%% Derived flattened quantities

rhoAll = hypot(uAll, vAll);

amplitudeAll = abs(VAll);
phaseRadAll = angle(VAll);
phaseDegAll = rad2deg(phaseRadAll);


uAll_Glambda = uAll/1e9;
vAll_Glambda = vAll/1e9;
wAll_Glambda = wAll/1e9;
rhoAll_Glambda = rhoAll/1e9;


%% Build per-baseline structures

baselineTemplate = struct( ...
    'field', '', ...
    'name', '', ...
    'accessibleIndex', [], ...
    'u', [], ...
    'v', [], ...
    'w', [], ...
    'rho', [], ...
    'u_Glambda', [], ...
    'v_Glambda', [], ...
    'w_Glambda', [], ...
    'rho_Glambda', [], ...
    'V', [], ...
    'amplitude', [], ...
    'phase_rad', [], ...
    'phase_deg', []);


baselines = repmat( ...
    baselineTemplate, ...
    nBaselines, ...
    1);


baselineIndexAll = ...
    zeros(numel(uAll), 1);

baselineNameAll = ...
    strings(numel(uAll), 1);


offset = 0;


for b = 1:nBaselines

    nThis = sampleCount(b);

    idx = ...
        offset + (1:nThis);


    baselines(b).field = ...
        baselineFields{b};

    baselines(b).name = ...
        baselineNames{b};

    baselines(b).accessibleIndex = ...
        accessibleIndex;

    baselines(b).u = ...
        uAll(idx);

    baselines(b).v = ...
        vAll(idx);

    baselines(b).w = ...
        wAll(idx);

    baselines(b).rho = ...
        rhoAll(idx);

    baselines(b).u_Glambda = ...
        uAll_Glambda(idx);

    baselines(b).v_Glambda = ...
        vAll_Glambda(idx);

    baselines(b).w_Glambda = ...
        wAll_Glambda(idx);

    baselines(b).rho_Glambda = ...
        rhoAll_Glambda(idx);

    baselines(b).V = ...
        VAll(idx);

    baselines(b).amplitude = ...
        amplitudeAll(idx);

    baselines(b).phase_rad = ...
        phaseRadAll(idx);

    baselines(b).phase_deg = ...
        phaseDegAll(idx);


    baselineIndexAll(idx) = b;

    baselineNameAll(idx) = ...
        string(baselineNames{b});


    offset = offset + nThis;

end


%% Assemble output

obs = struct();

if isfield(simulated, 'type')

    obs.targetType = simulated.type;

else

    obs.targetType = '';

end


obs.accessibleIndex = ...
    accessibleIndex;

obs.nAccessibleTimeSamples = ...
    nAccessible;

obs.nBaselines = ...
    nBaselines;

obs.baselineFields = ...
    baselineFields;

obs.baselineNames = ...
    baselineNames;

obs.baselines = ...
    baselines;


% Flattened physical observation set.

obs.u = uAll;
obs.v = vAll;
obs.w = wAll;

obs.rho = rhoAll;

obs.u_Glambda = uAll_Glambda;
obs.v_Glambda = vAll_Glambda;
obs.w_Glambda = wAll_Glambda;
obs.rho_Glambda = rhoAll_Glambda;

obs.V = VAll;

obs.amplitude = amplitudeAll;
obs.phase_rad = phaseRadAll;
obs.phase_deg = phaseDegAll;

obs.baselineIndex = ...
    baselineIndexAll;

obs.baselineName = ...
    baselineNameAll;

obs.nPhysicalSamples = ...
    numel(VAll);

obs.blockSize = ...
    blockSize;


% Image-grid support diagnostics.

obs.gridCheck = ...
    gridCheck;

end


%% ========================================================================
function V = local_sample_visibility( ...
    simulated, ...
    u, ...
    v, ...
    blockSize)
%LOCAL_SAMPLE_VISIBILITY
%
% Memory-bounded direct nonuniform discrete Fourier transform:
%
%   V_k = sum_p I_p exp[-2*pi*i*(u_k*l_p + v_k*m_p)].


I = double(simulated.I);

l = double(simulated.l_rad(:));
m = double(simulated.m_rad(:));


if size(I, 2) ~= numel(l) || ...
   size(I, 1) ~= numel(m)

    error('Image dimensions are inconsistent with angular axes.');

end


if any(~isfinite(I(:)))

    error('simulated.I contains nonfinite values.');

end


totalFlux = sum(I(:));


if ~isfinite(totalFlux) || ...
   abs(totalFlux) <= eps

    error('simulated.I has zero or invalid total flux.');

end


% The current pipeline uses normalized morphology. Normalize defensively
% here so V(0,0) = 1.

I = I/totalFlux;


[L, M] = meshgrid(l, m);


lPix = L(:);
mPix = M(:);
IPix = I(:);


% Zero pixels do not contribute and can be removed exactly.

active = IPix ~= 0;

lPix = lPix(active);
mPix = mPix(active);
IPix = IPix(active);


u = u(:);
v = v(:);


nVis = numel(u);

V = complex(zeros(nVis, 1));


% Limit temporary phase matrices to a few million entries.

maxPhaseElements = 2e6;


for k0 = 1:blockSize:nVis

    k1 = min( ...
        k0 + blockSize - 1, ...
        nVis);

    idxVis = k0:k1;

    uBlock = u(idxVis);
    vBlock = v(idxVis);

    nBlock = numel(idxVis);


    pixelChunk = max( ...
        1, ...
        floor(maxPhaseElements/nBlock));


    VBlock = complex( ...
        zeros(nBlock, 1));


    for p0 = 1:pixelChunk:numel(IPix)

        p1 = min( ...
            p0 + pixelChunk - 1, ...
            numel(IPix));

        idxPix = p0:p1;


        phase = ...
            uBlock*lPix(idxPix).' + ...
            vBlock*mPix(idxPix).';


        VBlock = ...
            VBlock + ...
            exp(-2*pi*1i*phase)*IPix(idxPix);

    end


    V(idxVis) = VBlock;

end

end


%% ========================================================================
function uvw = local_extract_uvw(uvFixed, fieldName)
%LOCAL_EXTRACT_UVW
%
% Extract one baseline as an N x 3 [u v w] array.


if ~isfield(uvFixed, fieldName)

    error('uvFixed.%s does not exist.', fieldName);

end


data = uvFixed.(fieldName);


if isnumeric(data)

    if size(data, 2) == 3

        uvw = data;

    elseif size(data, 1) == 3

        uvw = data.';

    else

        error([ ...
            'uvFixed.%s must be N x 3 or 3 x N.'], ...
            fieldName);

    end


elseif isstruct(data)

    if isfield(data, 'u') && ...
       isfield(data, 'v')

        u = data.u(:);
        v = data.v(:);


        if isfield(data, 'w')

            w = data.w(:);

        else

            w = zeros(size(u));

        end


        if numel(u) ~= numel(v) || ...
           numel(u) ~= numel(w)

            error([ ...
                'u, v, and w fields in uvFixed.%s ', ...
                'must have equal lengths.'], ...
                fieldName);

        end


        uvw = [u, v, w];


    elseif isfield(data, 'uvw')

        uvw = data.uvw;


        if size(uvw, 2) ~= 3 && ...
           size(uvw, 1) == 3

            uvw = uvw.';

        end


    else

        error([ ...
            'Cannot interpret uvFixed.%s as UVW data.'], ...
            fieldName);

    end


else

    error([ ...
        'Unsupported type for uvFixed.%s.'], ...
        fieldName);

end


if size(uvw, 2) ~= 3

    error([ ...
        'Extracted uvFixed.%s must contain three columns [u v w].'], ...
        fieldName);

end


if any(~isfinite(uvw(:)))

    error('uvFixed.%s contains nonfinite values.', fieldName);

end

end


%% ========================================================================
function idx = local_accessible_index(access, nTime)
%LOCAL_ACCESSIBLE_INDEX
%
% Resolve the accessible time-sample indices from the access structure.


if islogical(access) || ...
   (isnumeric(access) && isvector(access))

    values = access(:);


    if numel(values) == nTime

        idx = find(values ~= 0);

    else

        idx = values;

    end


elseif isstruct(access)

    indexFields = { ...
        'accessibleIndex', ...
        'idxAccessible', ...
        'accessibleIdx', ...
        'indexAccessible'};


    idx = [];


    for k = 1:numel(indexFields)

        fieldName = indexFields{k};


        if isfield(access, fieldName) && ...
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


            if isfield(access, fieldName)

                mask = access.(fieldName);


                if numel(mask) ~= nTime

                    continue

                end


                idx = find(mask(:) ~= 0);

                break

            end

        end

    end


    if isempty(idx) && ...
       isfield(access, 'thetaOffDeg') && ...
       isfield(access, 'maxOffNormalDeg')

        thetaOffDeg = ...
            access.thetaOffDeg(:);


        if numel(thetaOffDeg) == nTime

            idx = find( ...
                thetaOffDeg <= ...
                access.maxOffNormalDeg);

        end

    end


    if isempty(idx)

        error([ ...
            'Could not determine accessible samples from ', ...
            'the access structure.']);

    end


else

    error('Unsupported accessibility input.');

end


idx = double(idx(:));


if any(~isfinite(idx)) || ...
   any(idx < 1) || ...
   any(idx > nTime) || ...
   any(abs(idx - round(idx)) > 0)

    error('Accessible indices are invalid.');

end


idx = unique(round(idx), 'stable');

end