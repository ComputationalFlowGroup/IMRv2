% file f_odesolve.m
% brief contains function f_odesolve

% brief This function computes time marching for the ODE and PDE system of
% equations computed as part of the IMR solver. The function features known
% integration functions: ODE15, ODE23tb (most stable), and ODE45.
function [t,X] = f_odesolve(bubble, init, method, divisions, Tol, tspan, maxwalltime)
    if nargin < 7 || isempty(maxwalltime)
        maxwalltime = 0;
    end
    
    if divisions == 0
        options = odeset('RelTol',Tol(1),'AbsTol',Tol(2));
    else
        options = odeset('MaxStep',tspan(end)/divisions,'RelTol',Tol(1), ...
            'AbsTol',Tol(2));
    end
    if isfinite(maxwalltime) && maxwalltime > 0
        startTime = tic;
        options = odeset(options, 'OutputFcn', ...
            @(t, y, flag) stopOnWallTime(t, y, flag, startTime, maxwalltime));
    end
    
    if method == 15
        [t,X] = ode15s(bubble,tspan,init,options);
    elseif method == 23
        [t,X] = ode23tb(bubble,tspan,init,options);
    elseif method == 45
        [t,X] = ode45(bubble,tspan,init,options);
    else
        error('f_odesolve not set for the run, check casefile');
    end
    
end

function status = stopOnWallTime(~, ~, flag, startTime, maxwalltime)
    status = 0;
    if isempty(flag) && toc(startTime) >= maxwalltime
        status = 1;
    end
end
