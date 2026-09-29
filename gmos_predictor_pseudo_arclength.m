function z_pred = gmos_predictor_pseudo_arclength(z_prev, z_prevprev, ds)
%GMOS_PREDICTOR_PSEUDO_ARCLENGTH Secant predictor for continuation.
%
% Inputs
%   z_prev     : previous converged solution
%   z_prevprev : solution before z_prev
%   ds         : step size
%
% Output
%   z_pred     : predicted next solution

    t = z_prev - z_prevprev;
    nt = norm(t);

    if nt == 0
        error('Predictor: previous solutions are identical.');
    end

    t_hat = t / nt;
    z_pred = z_prev + ds * t_hat;
end