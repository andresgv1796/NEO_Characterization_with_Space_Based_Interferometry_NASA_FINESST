function [s1, Jrow, t_hat_s] = gmos_constraint_pseudo_arclength_weighted(y, y_prev, y_prevprev, ds, w)
%GMOS_CONSTRAINT_PSEUDO_ARCLENGTH_WEIGHTED Weighted PALC constraint in scaled space.
%
% We define the scaled variable ys = D*y with diagonal weight vector w,
% meaning ys = w .* y. The pseudo arclength constraint is
%
%   s1(y) = t_hat_s' * (ys - ys_prev) - ds = 0,
%
% where
%
%   t_hat_s = (ys_prev - ys_prevprev) / ||ys_prev - ys_prevprev||_2.
%
% Inputs
%   y, y_prev, y_prevprev : extended unknowns [X0(:); T; rho; mu]
%   ds                    : pseudo arclength step in the scaled space
%   w                     : diagonal scaling weights, same length as y
%
% Outputs
%   s1      : scalar PALC residual
%   Jrow    : 1 x numel(y) Jacobian row ds1/dy
%   t_hat_s : normalized secant in the scaled space

    y = y(:);
    y_prev = y_prev(:);
    y_prevprev = y_prevprev(:);
    w = w(:);

    if ~(numel(y) == numel(y_prev) && numel(y) == numel(y_prevprev) && numel(y) == numel(w))
        error('y, y_prev, y_prevprev, and w must all have the same length.');
    end

    dy_prev_s = w .* (y_prev - y_prevprev);
    nd = norm(dy_prev_s);
    if ~isfinite(nd) || nd == 0
        error('Weighted PALC: previous scaled secant has zero or invalid norm.');
    end

    t_hat_s = dy_prev_s / nd;

    s1 = t_hat_s.' * (w .* (y - y_prev)) - ds;
    Jrow = (t_hat_s .* w).';
end
