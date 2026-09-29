function s0 = gmos_constraint_fix_T(T, Ttarget)
%GMOS_CONSTRAINT_FIX_T Fix the stroboscopic time (paper s0(T)).
%
% Inputs
%   T       : current stroboscopic time
%   Ttarget : desired stroboscopic time
%
% Output
%   s0      : T - Ttarget

    s0 = T - Ttarget;
end