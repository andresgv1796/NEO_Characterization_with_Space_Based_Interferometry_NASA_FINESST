function shape = generate_triaxial_ellipsoid(a, b, c, N)
%GENERATE_TRIAXIAL_ELLIPSOID Generate a triangulated triaxial ellipsoid.

arguments
    a (1,1) double {mustBePositive}
    b (1,1) double {mustBePositive}
    c (1,1) double {mustBePositive}
    N (1,1) double {mustBeInteger, mustBePositive} = 40
end

% Generate ellipsoid surface
[X, Y, Z] = sphere(N);

X = a * X;
Y = b * Y;
Z = c * Z;

% Convert to triangular mesh
mesh = surf2patch(X, Y, Z, 'triangles');

V = mesh.vertices;
F = mesh.faces;

% Remove degenerate zero-area triangles
v1 = V(F(:,1), :);
v2 = V(F(:,2), :);
v3 = V(F(:,3), :);

doubleArea = vecnorm( ...
    cross(v2 - v1, v3 - v1, 2), ...
    2, 2);

tol = 1e-12 * max(doubleArea);

keep = doubleArea > tol;

F = F(keep, :);

% Store
shape.vertices = V;
shape.faces    = F;
shape.axes     = [a, b, c];
shape.N        = N;

end