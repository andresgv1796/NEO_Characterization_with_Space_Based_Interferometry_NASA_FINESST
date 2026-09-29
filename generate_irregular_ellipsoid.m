function shape = generate_irregular_ellipsoid(a, b, c, N, irregularity)
%GENERATE_IRREGULAR_ELLIPSOID
% Generate a smooth irregular asteroid from a perturbed triaxial ellipsoid.
%
% Inputs:
%   a,b,c          - nominal semi-axes [km]
%   N              - sphere mesh resolution
%   irregularity   - dimensionless perturbation amplitude
%                    e.g. 0.10 = mild, 0.25 = pronounced
%
% Output:
%   shape.vertices
%   shape.faces
%   shape.axes
%   shape.irregularity
%   shape.N

arguments
    a (1,1) double {mustBePositive}
    b (1,1) double {mustBePositive}
    c (1,1) double {mustBePositive}
    N (1,1) double {mustBeInteger,mustBePositive} = 60
    irregularity (1,1) double {mustBeNonnegative} = 0.15
end

% Unit-sphere parameterization
[Xs, Ys, Zs] = sphere(N);

theta = acos(max(-1,min(1,Zs)));
phi   = atan2(Ys,Xs);

% Smooth low-order shape perturbation
f = ...
      0.50*sin(2*theta).*cos(3*phi) ...
    + 0.30*sin(3*theta).*sin(2*phi) ...
    + 0.20*cos(4*theta);

% Normalize perturbation so irregularity directly controls amplitude
f = f / max(abs(f(:)));

R = 1 + irregularity*f;

% Perturbed triaxial ellipsoid
X = a * R .* Xs;
Y = b * R .* Ys;
Z = c * R .* Zs;

% Triangulate
mesh = surf2patch(X,Y,Z,'triangles');

V = mesh.vertices;
F = mesh.faces;

% Remove degenerate faces
v1 = V(F(:,1),:);
v2 = V(F(:,2),:);
v3 = V(F(:,3),:);

doubleArea = vecnorm( ...
    cross(v2-v1,v3-v1,2), ...
    2,2);

tol = 1e-12*max(doubleArea);

F = F(doubleArea > tol,:);

shape.vertices     = V;
shape.faces        = F;
shape.axes         = [a,b,c];
shape.irregularity = irregularity;
shape.N            = N;

end