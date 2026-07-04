% file f_odesolve.m
% brief contains function f_odesolve

% brief This function computes time marching for the ODE and PDE system of
% equations computed as part of the IMR solver. The function features known
% integration functions: ODE15, ODE23tb (most stable), and ODE45.

% event_fcn (optional): a MATLAB ODE event function handle, e.g.
%   @(t,X) deal(X(2), 1, 1)
% passed straight through to odeset('Events', event_fcn). Omit or pass []
% to get the previous (no early stopping) behavior unchanged. This is used
% to stop the integration exactly at the collapse minimum radius (first
% Rdot zero-crossing from negative to positive) instead of relying on
% post-processing a fixed output time grid -- see f_imr_fd.m's
% f_Rmin_event for the concrete event function used for that purpose.
function [t,X] = f_odesolve(bubble, init, method, divisions, tspan, event_fcn)
    if nargin < 6
	event_fcn = [];
    end

    if divisions == 0
        options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    else
        options = odeset('MaxStep',tspan(end)/divisions,'RelTol',1e-12);
    end

    if ~isempty(event_fcn)
	options = odeset(options, 'Events', event_fcn);
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
