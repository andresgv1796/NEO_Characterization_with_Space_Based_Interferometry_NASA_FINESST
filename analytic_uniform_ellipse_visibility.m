function [V, q] = analytic_uniform_ellipse_visibility(u, v, a_rad, b_rad)
%ANALYTIC_UNIFORM_ELLIPSE_VISIBILITY
% Normalized visibility of a centered, axis-aligned uniform ellipse.
%
% Inputs:
%   u, v   - spatial frequencies [wavelengths]
%   a_rad  - angular semi-axis along l [rad]
%   b_rad  - angular semi-axis along m [rad]
%
% Outputs:
%   V      - normalized real visibility
%   q      - elliptical spatial-frequency coordinate

u = u(:);
v = v(:);

q = sqrt( ...
    (a_rad*u).^2 + ...
    (b_rad*v).^2);

x = 2*pi*q;

V = ones(size(x));

idx = abs(x) > 1e-12;

V(idx) = 2*besselj(1, x(idx)) ./ x(idx);

end