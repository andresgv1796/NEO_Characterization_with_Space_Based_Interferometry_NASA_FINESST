function check = validate_image_uv_sampling(simulated, u, v)
%VALIDATE_IMAGE_UV_SAMPLING
%
% Check whether a regularly sampled angular image grid can support the
% requested interferometric spatial-frequency samples.
%
% For angular grid spacings
%
%       Delta_l = spacing in l [rad]
%       Delta_m = spacing in m [rad]
%
% the Cartesian Nyquist limits are
%
%       u_Nyq = 1/(2*Delta_l)
%       v_Nyq = 1/(2*Delta_m).
%
% Define
%
%       q_u = max|u| / u_Nyq = 2*Delta_l*max|u|
%       q_v = max|v| / v_Nyq = 2*Delta_m*max|v|.
%
% Formal sampling validity requires
%
%       q_u <= 1
%       q_v <= 1.
%
% The image-grid Fourier support is therefore rectangular:
%
%       |u| <= u_Nyq
%       |v| <= v_Nyq.
%
% The number of image samples across the finest Fourier fringe is
%
%       Nfringe_u = 1/(Delta_l*max|u|) = 2/q_u
%       Nfringe_v = 1/(Delta_m*max|v|) = 2/q_v.
%
% A working recommendation of q <= 0.5 corresponds to at least four image
% samples across the finest fringe. This is a numerical margin rather than
% a hard physical requirement.
%
% INPUTS
%   simulated : simulated target structure containing
%               I, l_rad, m_rad
%   u         : sampled u coordinates [wavelengths]
%   v         : sampled v coordinates [wavelengths]
%
% OUTPUT
%   check     : structure containing all sampling diagnostics


%% Validate inputs

requiredFields = {'I', 'l_rad', 'm_rad'};

for k = 1:numel(requiredFields)

    if ~isfield(simulated, requiredFields{k})

        error('simulated.%s is required.', ...
            requiredFields{k});

    end

end


l = simulated.l_rad(:);
m = simulated.m_rad(:);

u = u(:);
v = v(:);


if numel(u) ~= numel(v)

    error('u and v must contain the same number of samples.');

end


if isempty(u)

    error('At least one UV sample is required.');

end


if any(~isfinite(u)) || any(~isfinite(v))

    error('u and v must contain only finite values.');

end


if numel(l) < 2 || numel(m) < 2

    error('Image axes must contain at least two samples.');

end


if size(simulated.I, 2) ~= numel(l) || ...
   size(simulated.I, 1) ~= numel(m)

    error([ ...
        'simulated.I dimensions are inconsistent with ', ...
        'simulated.l_rad and simulated.m_rad.']);

end


%% Check image-grid regularity

dl = diff(l);
dm = diff(m);


if any(dl == 0) || any(dm == 0)

    error('Image axes must be strictly monotonic.');

end


if any(sign(dl) ~= sign(dl(1)))

    error('simulated.l_rad must be monotonic.');

end


if any(sign(dm) ~= sign(dm(1)))

    error('simulated.m_rad must be monotonic.');

end


deltaL = mean(abs(dl));
deltaM = mean(abs(dm));


relativeNonuniformityL = ...
    max(abs(abs(dl) - deltaL))/deltaL;

relativeNonuniformityM = ...
    max(abs(abs(dm) - deltaM))/deltaM;


uniformityTolerance = 1e-9;


if relativeNonuniformityL > uniformityTolerance

    error('simulated.l_rad is not uniformly sampled.');

end


if relativeNonuniformityM > uniformityTolerance

    error('simulated.m_rad is not uniformly sampled.');

end


%% Cartesian Nyquist limits

uNyq = 1/(2*deltaL);
vNyq = 1/(2*deltaM);


uMax = max(abs(u));
vMax = max(abs(v));


qU = uMax/uNyq;
qV = vMax/vNyq;

qMax = max(qU, qV);


%% Samples across finest Fourier fringe

if uMax > 0

    pixelsPerFringeU = 1/(uMax*deltaL);

else

    pixelsPerFringeU = Inf;

end


if vMax > 0

    pixelsPerFringeV = 1/(vMax*deltaM);

else

    pixelsPerFringeV = Inf;

end


pixelsPerFringeMin = min( ...
    pixelsPerFringeU, ...
    pixelsPerFringeV);


%% Image dimensions and field of view

nL = numel(l);
nM = numel(m);

fullWidthL_rad = abs(l(end) - l(1));
fullWidthM_rad = abs(m(end) - m(1));


masPerRad = ...
    180*3600*1000/pi;


deltaL_mas = deltaL*masPerRad;
deltaM_mas = deltaM*masPerRad;

fullWidthL_mas = fullWidthL_rad*masPerRad;
fullWidthM_mas = fullWidthM_rad*masPerRad;


%% Minimum grid for the current field of view

% Formal Nyquist:
%
%       2*Delta_l*uMax <= 1
%
% with
%
%       Delta_l = fullWidthL/(N_l - 1)

nLMinNyquist = ceil( ...
    1 + 2*fullWidthL_rad*uMax);

nMMinNyquist = ceil( ...
    1 + 2*fullWidthM_rad*vMax);


nLMinNyquist = max(nLMinNyquist, 3);
nMMinNyquist = max(nMMinNyquist, 3);


% Use odd recommendations so that the image origin can lie exactly on a
% pixel for symmetric grids.

if mod(nLMinNyquist, 2) == 0

    nLMinNyquist = nLMinNyquist + 1;

end


if mod(nMMinNyquist, 2) == 0

    nMMinNyquist = nMMinNyquist + 1;

end


%% Recommended working grid

qRecommended = 0.5;


nLRecommended = ceil( ...
    1 + ...
    2*fullWidthL_rad*uMax/qRecommended);

nMRecommended = ceil( ...
    1 + ...
    2*fullWidthM_rad*vMax/qRecommended);


nLRecommended = max(nLRecommended, 3);
nMRecommended = max(nMRecommended, 3);


if mod(nLRecommended, 2) == 0

    nLRecommended = nLRecommended + 1;

end


if mod(nMRecommended, 2) == 0

    nMRecommended = nMRecommended + 1;

end


%% Classify

numericalTolerance = 1e-12;


isNyquistValid = ...
    qU <= 1 + numericalTolerance && ...
    qV <= 1 + numericalTolerance;


isWellSampled = ...
    qU <= qRecommended && ...
    qV <= qRecommended;


if ~isNyquistValid

    status = 'INVALID';

elseif isWellSampled

    status = 'WELL SAMPLED';

else

    status = 'NYQUIST VALID - LIMITED MARGIN';

end


%% Store diagnostics

check = struct();

check.deltaL_rad = deltaL;
check.deltaM_rad = deltaM;

check.deltaL_mas = deltaL_mas;
check.deltaM_mas = deltaM_mas;

check.fullWidthL_rad = fullWidthL_rad;
check.fullWidthM_rad = fullWidthM_rad;

check.fullWidthL_mas = fullWidthL_mas;
check.fullWidthM_mas = fullWidthM_mas;

check.nL = nL;
check.nM = nM;

check.uNyq = uNyq;
check.vNyq = vNyq;

check.uNyq_Glambda = uNyq/1e9;
check.vNyq_Glambda = vNyq/1e9;

check.uMax = uMax;
check.vMax = vMax;

check.uMax_Glambda = uMax/1e9;
check.vMax_Glambda = vMax/1e9;

check.qU = qU;
check.qV = qV;
check.qMax = qMax;

check.pixelsPerFringeU = pixelsPerFringeU;
check.pixelsPerFringeV = pixelsPerFringeV;
check.pixelsPerFringeMin = pixelsPerFringeMin;

check.qRecommended = qRecommended;

check.nLMinNyquist = nLMinNyquist;
check.nMMinNyquist = nMMinNyquist;

check.nLRecommended = nLRecommended;
check.nMRecommended = nMRecommended;

check.isNyquistValid = isNyquistValid;
check.isWellSampled = isWellSampled;

check.status = status;


%% Hard failure only for an actual Cartesian Nyquist violation

if ~isNyquistValid

    error([ ...
        'Image grid does not support the requested UV samples.\n', ...
        '  q_u = max|u|/u_Nyq = %.6f\n', ...
        '  q_v = max|v|/v_Nyq = %.6f\n', ...
        '  current grid         = %d x %d\n', ...
        '  minimum Nyquist grid = %d x %d\n', ...
        'Increase image resolution, reduce angular FOV, ', ...
        'or reduce the sampled UV extent.'], ...
        qU, ...
        qV, ...
        nL, ...
        nM, ...
        nLMinNyquist, ...
        nMMinNyquist);

end


%% Warn when formally valid but weakly oversampled

if isNyquistValid && ~isWellSampled

    warning([ ...
        'Image grid is Nyquist-valid but only weakly oversampled.\n', ...
        '  q_u = %.4f -> %.2f pixels/fringe\n', ...
        '  q_v = %.4f -> %.2f pixels/fringe\n', ...
        '  recommended q <= %.2f\n', ...
        '  recommended grid for current FOV = %d x %d\n', ...
        'Visibility convergence should be checked for production results.'], ...
        qU, ...
        pixelsPerFringeU, ...
        qV, ...
        pixelsPerFringeV, ...
        qRecommended, ...
        nLRecommended, ...
        nMRecommended);

end

end