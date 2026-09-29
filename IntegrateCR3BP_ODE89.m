function [t_eval, Y_eval] = IntegrateCR3BP_ODE89(Y0, tspan, mu, N)
%IntegrateCR3BP_ODE89  Integrate Circular Restricted 3 body problem in Hill frame
%   [t_eval, Y_eval] = IntegrateCR3BP_ODE89(Y0, tspan, mu, N) returns the state
%   history Y_eval evaluated at evenly spaced times between tspan(1) and
%   tspan(2) using the ode89 solver with tight tolerances.
%
%   Inputs
%     Y0     6×1 initial state vector [position; velocity]
%     tspan  [t0, tf] integration interval
%     mu     Mass ratio of the two primaries (scalar)
%     N      number of output points
%   Outputs
%     t_eval   N×1 vector of times where solution is evaluated
%     Y_eval   N×6 matrix of state vectors at each t_eval

   

    % build evaluation grid
    t_eval = linspace(tspan(1), tspan(2), N)';

    % set tight relative and absolute tolerances
    options = odeset( ...
        'RelTol', 1e-12, ...
        'AbsTol', 1e-12 ...
    );

    % pack right hand side
    odefun = @(t, Y) CR3BP_mu(t, Y, mu);

    % perform numeric integration
    sol = ode89(odefun, [tspan(1), tspan(2)], Y0, options);

    % sample solution on our grid
    Y_eval = deval(sol, t_eval)';

end
