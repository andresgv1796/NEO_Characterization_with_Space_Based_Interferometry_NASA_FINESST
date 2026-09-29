function [y_pred, t_hat_s] = gmos_predictor_pseudo_arclength_weighted(y_prev, y_prevprev, ds, w)
%GMOS_PREDICTOR_PSEUDO_ARCLENGTH_WEIGHTED Weighted secant predictor in scaled space.
%
% Inputs
%   y_prev, y_prevprev : previous two converged solutions [X0(:); T; rho; mu]
%   ds                 : pseudo arclength step in the scaled space
%   w                  : diagonal weight vector, so ys = w .* y
%
% Outputs
%   y_pred  : predicted next solution in physical coordinates
%   t_hat_s : normalized secant in the scaled space

    [~, Jrow, t_hat_s] = gmos_constraint_pseudo_arclength_weighted(y_prev, y_prev, y_prevprev, 0.0, w);

    %#ok<NASGU> Jrow is not used, but the call above consistently builds t_hat_s.
    y_pred = y_prev + ds * (t_hat_s ./ w);
end
