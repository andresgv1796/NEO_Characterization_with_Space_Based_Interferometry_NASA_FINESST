function dYdt = CR3BP_mu(t, Y, mu)
% CR3BP_mu  Equations of motion for the Circular Restricted Three-Body Problem (CR3BP)
%
%   dYdt = CR3BP_mu(t, Y, mu) computes the time derivatives of the state vector
%   Y for the CR3BP system in a rotating frame. The state vector Y is defined as:
%       Y = [x; y; z; vx; vy; vz]
%
%   INPUTS:
%       t  - Time variable (scalar). Not used explicitly here since the system is autonomous.
%       Y  - State vector [x; y; z; vx; vy; vz] (6x1)
%       mu - Mass ratio of the two primaries (scalar), typically the Earth-Moon mass ratio.
%
%   OUTPUT:
%       dYdt - Time derivative of the state vector, [vx; vy; vz; ax; ay; az] (6x1)
%
%   The equations of motion include:
%       - The Coriolis force (2*n*vy and -2*n*vx, where n is the angular rate)
%       - The centrifugal force (n^2*x and n^2*y)
%       - The gravitational attractions of the two primaries.
%
%   For the CR3BP, the primaries are located at (-mu, 0, 0) and (1-mu, 0, 0). 
%   The distances from the third body to the primaries are computed as:
%       r13 = sqrt((x+mu)^2 + y^2 + z^2)
%       r23 = sqrt((x-1+mu)^2 + y^2 + z^2)
%
%   The accelerations in the rotating frame are given by:
%       xddot = 2*n*vy + n^2*x - ((1-mu)*(x+mu))/(r13^3) - (mu*(x-1+mu))/(r23^3)
%       yddot = -2*n*vx + n^2*y - ((1-mu)*y)/(r13^3) - (mu*y)/(r23^3)
%       zddot = -((1-mu)*z)/(r13^3) - (mu*z)/(r23^3)


    % Angular rate in the rotating frame (n) is often normalized to 1.
    n = 1;

    % Unpack the state vector.
    x = Y(1);
    y = Y(2);
    z = Y(3);
    vx = Y(4);
    vy = Y(5);
    vz = Y(6);

    % Compute the distances from the third body to the two primaries.
    r13 = sqrt((x + mu)^2 + y^2 + z^2);      % Distance to primary 1 (located at (-mu, 0, 0))
    r23 = sqrt((x - 1 + mu)^2 + y^2 + z^2);    % Distance to primary 2 (located at (1-mu, 0, 0))

    % Compute the accelerations in the rotating frame.
    xddot = 2 * n * vy + n^2 * x - ((1 - mu) * (x + mu)) / (r13^3) - (mu * (x - 1 + mu)) / (r23^3);
    yddot = -2 * n * vx + n^2 * y - ((1 - mu) * y) / (r13^3) - (mu * y) / (r23^3);
    zddot = - ((1 - mu) * z) / (r13^3) - (mu * z) / (r23^3);

    % Pack the derivatives into the output vector.
    dYdt = [vx; vy; vz; xddot; yddot; zddot];
end
