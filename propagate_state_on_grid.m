function [Y, sol] = propagate_state_on_grid(x0, tGrid, mu, odeOpts)
%PROPAGATE_STATE_ON_GRID Integrate one CR3BP state and sample on tGrid.
%
% Legacy usage:
%
%   Y = propagate_state_on_grid(x0,tGrid,mu,odeOpts)
%
% preserves the existing project behavior and preferentially uses
% IntegrateCR3BP_ODE89 when available.
%
% Dense-solution usage:
%
%   [Y,sol] = propagate_state_on_grid(x0,tGrid,mu,odeOpts)
%
% uses a direct MATLAB ODE integration and returns the dense solution
% structure sol, which can subsequently be evaluated using deval.
%
% No integration failure is silently caught. If an available integrator
% fails, the error propagates to the caller.

    x0 = x0(:);
    tGrid = tGrid(:);

    if numel(x0) ~= 6
        error('x0 must contain exactly 6 CR3BP state components.');
    end

    if any(~isfinite(x0))
        error('x0 contains nonfinite values.');
    end

    if numel(tGrid) < 2
        error('tGrid must contain at least two time points.');
    end

    if any(~isfinite(tGrid))
        error('tGrid contains nonfinite values.');
    end

    if any(diff(tGrid) <= 0)
        error('tGrid must be strictly increasing.');
    end

    if nargin < 4 || isempty(odeOpts)
        odeOpts = odeset( ...
            'RelTol',1e-12, ...
            'AbsTol',1e-12);
    end

    nGrid = numel(tGrid);
    tspan = [tGrid(1),tGrid(end)];

    wantDenseSolution = nargout >= 2;


    %% Legacy one-output path
    %
    % Preserve the existing project integration convention exactly when
    % only Y is requested.

    if ~wantDenseSolution && exist('IntegrateCR3BP_ODE89','file') == 2

        [tEval,Y] = ...
            IntegrateCR3BP_ODE89( ...
                x0, ...
                tspan, ...
                mu, ...
                nGrid);

        tEval = tEval(:);

        if size(Y,1) ~= numel(tEval) || size(Y,2) ~= 6
            error([ ...
                'IntegrateCR3BP_ODE89 returned inconsistent output ', ...
                'dimensions. Expected nTime x 6 states.']);
        end

        if any(~isfinite(Y(:)))
            error('IntegrateCR3BP_ODE89 returned nonfinite state values.');
        end

        if any(diff(tEval) <= 0)
            error('IntegrateCR3BP_ODE89 returned a nonmonotonic time vector.');
        end

        timeTolerance = ...
            1e-10*max(1,max(abs(tGrid)));

        gridMismatch = ...
            numel(tEval) ~= nGrid;

        if ~gridMismatch
            gridMismatch = ...
                max(abs(tEval-tGrid)) > timeTolerance;
        end

        if gridMismatch

            if tGrid(1) < tEval(1) || tGrid(end) > tEval(end)
                error([ ...
                    'Requested tGrid extends outside the time interval ', ...
                    'returned by IntegrateCR3BP_ODE89.']);
            end

            warning( ...
                'propagate_state_on_grid:GridMismatch', ...
                ['IntegrateCR3BP_ODE89 did not return the requested ', ...
                 'tGrid. Applying the legacy PCHIP resampling step.']);

            Y = interp1( ...
                tEval, ...
                Y, ...
                tGrid, ...
                'pchip');

            if any(~isfinite(Y(:)))
                error('PCHIP resampling produced nonfinite state values.');
            end

        end

        sol = [];
        return

    end


    %% Dense-solution path
    %
    % This path is used whenever a second output is requested.
    % It intentionally bypasses IntegrateCR3BP_ODE89 because that wrapper
    % does not expose a dense solution object suitable for deval.

    rhs = @(t,y) CR3BP_mu(t,y,mu);

    if exist('ode89','file') == 2

        sol = ode89( ...
            rhs, ...
            tspan, ...
            x0, ...
            odeOpts);

    else

        warning( ...
            'propagate_state_on_grid:NoODE89', ...
            'ode89 not found. Using ode113 for the dense solution.');

        sol = ode113( ...
            rhs, ...
            tspan, ...
            x0, ...
            odeOpts);

    end


    %% Evaluate dense solution on requested grid

    Y = deval(sol,tGrid).';

    if size(Y,1) ~= nGrid || size(Y,2) ~= 6
        error('Dense ODE evaluation returned an unexpected state-array size.');
    end

    if any(~isfinite(Y(:)))
        error('Dense ODE evaluation returned nonfinite state values.');
    end

end