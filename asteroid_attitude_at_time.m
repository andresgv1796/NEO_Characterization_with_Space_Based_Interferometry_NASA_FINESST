function [R_IB,phaseWrapped,phaseUnwrapped] = ...
    asteroid_attitude_at_time(rot,t_s)
%ASTEROID_ATTITUDE_AT_TIME
% Body-to-inertial attitude of a constant-rate rotating asteroid.
%
% Convention:
%
%       r_I = R_IB(t) r_B
%
% and
%
%       R_IB(t)
%       =
%       R_IB0 * Rot(spinAxis_B,phi)
%
% where
%
%       phi
%       =
%       spinSense * 2*pi/period * (t-t0)
%
% INPUTS
%   rot.period_s
%   rot.spinAxis_B
%   rot.R_IB0
%   rot.t0_s
%   rot.spinSense
%
%   t_s
%       scalar time [s]
%
% OUTPUTS
%   R_IB
%       3x3 body-to-inertial attitude matrix
%
%   phaseWrapped
%       phase in [0,2*pi)
%
%   phaseUnwrapped
%       continuous phase


    %% Validate required fields

    requiredFields = { ...
        'period_s', ...
        'spinAxis_B', ...
        'R_IB0', ...
        't0_s', ...
        'spinSense'};

    for k = 1:numel(requiredFields)

        if ~isfield(rot,requiredFields{k})

            error( ...
                'rot.%s is missing.', ...
                requiredFields{k});

        end

    end


    %% Spin axis

    sB = ...
        rot.spinAxis_B(:);

    sNorm = ...
        norm(sB);

    if sNorm < 1e-14

        error( ...
            'rot.spinAxis_B must be nonzero.');

    end

    sB = ...
        sB/sNorm;


    %% Rotation phase

    omega = ...
        rot.spinSense * ...
        2*pi / rot.period_s;


    phaseUnwrapped = ...
        omega * ...
        (t_s - rot.t0_s);


    phaseWrapped = ...
        mod(phaseUnwrapped,2*pi);


    %% Rodrigues rotation matrix

    sx = sB(1);
    sy = sB(2);
    sz = sB(3);


    K = [ ...
         0, -sz,  sy;
        sz,   0, -sx;
       -sy,  sx,   0];


    c = ...
        cos(phaseUnwrapped);

    s = ...
        sin(phaseUnwrapped);


    Rspin = ...
        eye(3) + ...
        s*K + ...
        (1-c)*(K*K);


    %% Body -> inertial attitude

    R_IB = ...
        rot.R_IB0 * Rspin;

end