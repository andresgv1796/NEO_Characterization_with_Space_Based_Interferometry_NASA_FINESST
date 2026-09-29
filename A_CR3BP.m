function A = A_CR3BP(Y, mu)
    % Computes the Jacobian matrix A for the CR3BP at state Y.
    % Y = [x, y, z, vx, vy, vz] and mu is the mass parameter.
    %
    % The matrix A is defined as:
    %   A = [zeros(3), eye(3);
    %        U_xx  U_xy  U_xz  0   2   0;
    %        U_xy  U_yy  U_yz -2   0   0;
    %        U_xz  U_yz  U_zz  0   0   0]
    % where U_ij are the second partial derivatives of the effective potential.
    
    % Extract only the necessary state variables
    x = Y(1);
    y = Y(2); 
    z = Y(3);
    
    % Distance calculations from the two primaries
    r13 = sqrt((x + mu)^2 + y^2 + z^2);
    r23 = sqrt((x - 1 + mu)^2 + y^2 + z^2);
    
    % Compute second derivatives of the effective potential U
    U_xx = 1 - (1 - mu) / r13^3 - mu / r23^3 ...
           + 3 * (1 - mu) * (x + mu)^2 / r13^5 ...
           + 3 * mu * (x - 1 + mu)^2 / r23^5;
       
    U_xy = 3 * (1 - mu) * (x + mu) * y / r13^5 ...
           + 3 * mu * (x - 1 + mu) * y / r23^5;
       
    U_xz = 3 * (1 - mu) * (x + mu) * z / r13^5 ...
           + 3 * mu * (x - 1 + mu) * z / r23^5;
       
    U_yy = 1 - (1 - mu) / r13^3 - mu / r23^3 ...
           + 3 * (1 - mu) * y^2 / r13^5 ...
           + 3 * mu * y^2 / r23^5;
       
    U_yz = 3 * (1 - mu) * y * z / r13^5 ...
           + 3 * mu * y * z / r23^5;
       
    U_zz = -(1 - mu) / r13^3 - mu / r23^3 ...
           + 3 * (1 - mu) * z^2 / r13^5 ...
           + 3 * mu * z^2 / r23^5;
    
    % Assemble the A matrix (Jacobian of the CR3BP equations)
    A = [0,   0,   0,   1,  0,  0;
         0,   0,   0,   0,  1,  0;
         0,   0,   0,   0,  0,  1;
         U_xx, U_xy, U_xz, 0,  2,  0;
         U_xy, U_yy, U_yz, -2, 0,  0;
         U_xz, U_yz, U_zz, 0,  0,  0];
end
