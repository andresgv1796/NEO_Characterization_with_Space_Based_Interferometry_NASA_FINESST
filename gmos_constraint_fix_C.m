function s0 = gmos_constraint_fix_C(X0, mu, Ctarget)
%GMOS_CONSTRAINT_FIX_C Fix the average Jacobi constant of the invariant curve.
%
% Paper-style parametrizing equation for the GMOS solver when the Jacobi
% constant is held fixed and the stroboscopic time T is free:
%
%   s0(X0) = (1/N) * sum_j C(X0(:,j)) - Ctarget
%
% Inputs
%   X0      : 6xN curve samples
%   mu      : CR3BP mass parameter
%   Ctarget : desired Jacobi constant
%
% Output
%   s0      : scalar residual

    if size(X0, 1) ~= 6
        error('X0 must be 6xN.');
    end
    if ~isscalar(mu) || ~isscalar(Ctarget)
        error('mu and Ctarget must be scalars.');
    end

    Cvals = Jacobi_Constant(X0.', mu);
    s0 = mean(Cvals) - Ctarget;
end
