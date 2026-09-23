% file f_odesolve.m
% brief contains function f_odesolve

% brief This function computes time marching for the ODE and PDE system of
% equations computed as part of the IMR solver. The function features known
% integration functions: ODE15, ODE23tb (most stable), and ODE45. Optional
% inputs: Sp, either a Jacobian sparsity pattern (JPattern) or an analytic
% Jacobian function handle, and the relative and absolute tolerances.
function [t,X] = f_odesolve(bubble, init, method, divisions, tspan, Sp, rtol, atol)
    
    if nargin < 6
        Sp   = [];
    end
    if nargin < 7
        rtol = 1e-8;
    end
    if nargin < 8
        atol = 1e-8;
    end
    
    if divisions == 0
        options = odeset('RelTol',rtol,'AbsTol',atol);
    else
        options = odeset('MaxStep',tspan(end)/divisions,'RelTol',rtol);
    end
    if ~isempty(Sp)
        if isa(Sp,'function_handle')
            options = odeset(options,'Jacobian',Sp);
        else
            options = odeset(options,'JPattern',Sp);
        end
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
