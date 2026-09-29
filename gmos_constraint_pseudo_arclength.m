function [s1, t_hat] = gmos_constraint_pseudo_arclength(z, z_prev, z_prevprev, ds)
%GMOS_CONSTRAINT_PSEUDO_ARCLENGTH Pseudo arclength continuation constraint.
%
% Inputs
%   z          : current unknown vector
%   z_prev     : previous converged solution
%   z_prevprev : solution before z_prev
%   ds         : step size along the continuation direction
%
% Outputs
%   s1    : scalar constraint value t_hat'*(z - z_prev) - ds
%   t_hat : normalized secant direction

    t = z_prev - z_prevprev;
    nt = norm(t);

    if nt == 0
        error('Pseudo arclength: previous solutions are identical.');
    end

    t_hat = t / nt;

    s1 = t_hat.' * (z - z_prev) - ds;
end