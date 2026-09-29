function [F, cache] = gmos_residual_full_fixC(z, N, mu, Xref, Tref, rhoref, Ctarget, ode_opts)
%GMOS_RESIDUAL_FULL_FIXC Residual with invariance, phase, and fixed-C constraint.
%
% F = [Finv(:); p0; p1; s0], where:
%   Finv = R_{-rho}(Phi^T(X0)) - X0
%   p0, p1 are the paper phase conditions using (Xref, Tref, rhoref)
%   s0 = mean_j C(X0(:,j)) - Ctarget
%
% Inputs
%   z       : packed unknowns [X0(:); T; rho]
%   N       : number of curve points
%   mu      : CR3BP mass parameter
%   Xref    : 6xN reference curve (fixed during Newton)
%   Tref    : scalar reference stroboscopic time
%   rhoref  : scalar reference rotation number
%   Ctarget : scalar target Jacobi constant
%   ode_opts: (optional) ODE options
%
% Outputs
%   F       : residual vector
%   cache   : struct with useful intermediate values

    if nargin < 8
        ode_opts = [];
    end

    [X0, T, rho] = gmos_unpack_z(z, N);

    X1 = gmos_flowmap_curve_CR3BP(X0, T, mu, ode_opts);
    X1shift = gmos_shift_curve_fft(X1, rho);
    Finv = X1shift - X0;

    [p0, p1, tang0, tang1] = gmos_phase_conditions(X0, Xref, mu, Tref, rhoref);

    Cvals = Jacobi_Constant(X0.', mu);
    s0 = mean(Cvals) - Ctarget;

    F = [Finv(:); p0; p1; s0];

    cache = struct();
    cache.X0 = X0;
    cache.T = T;
    cache.rho = rho;
    cache.X1 = X1;
    cache.X1shift = X1shift;
    cache.Finv = Finv;
    cache.p0 = p0;
    cache.p1 = p1;
    cache.s0 = s0;
    cache.Cvals = Cvals;
    cache.tang0 = tang0;
    cache.tang1 = tang1;
end
