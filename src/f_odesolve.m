% file f_odesolve.m
% brief contains function f_odesolve

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function computes time marching for the ODE and PDE system of
% equations computed as part of the IMR solver with the stiff integrators
% ode15s, ode23tb (most stable) and CVODES; the thermal, mass transfer and
% stress fields make the system stiff, so non-stiff integrators are not
% offered. method is 15 (ode15s), 23 (ode23tb) or 3 (the compiled SUNDIALS
% integrator CVODES through the MATLAB ode object, solver cvodesstiff,
% MATLAB R2024a or later). Optional inputs: Sp, either a Jacobian sparsity
% pattern (JPattern) or a Jacobian function handle, and the relative and
% absolute tolerances.
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
    
    if method == 3
        [t,X] = f_cvodes(bubble, init, divisions, tspan, Sp, rtol, atol);
        return
    end
    
    if divisions == 0
        options = odeset('RelTol',rtol,'AbsTol',atol);
    else
        options = odeset('MaxStep',tspan(end)/divisions,'RelTol',rtol, ...
            'AbsTol',atol);
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
    else
        error(['f_odesolve: method must be 3 (CVODES), 15 (ode15s) or 23 ' ...
            '(ode23tb)']);
    end
    
end

% CVODES (SUNDIALS, compiled) through the ode object. Properties are set by
% assignment so that the file also parses where the ode object is absent.
% With two times in tspan the solver's own steps are returned, otherwise
% the solution at the times of tspan, as with ode15s and ode23tb.
function [t,X] = f_cvodes(bubble, init, divisions, tspan, Sp, rtol, atol)
    
    if exist('ode', 'class') ~= 8
        error(['f_odesolve: method 3 (CVODES) needs the MATLAB ode object ' ...
            'with the cvodesstiff solver (MATLAB R2024a or later), use ' ...
            'method 15 (ode15s) or 23 (ode23tb)']);
    end
    F = ode;
    F.ODEFcn = bubble;
    F.InitialTime = tspan(1);
    F.InitialValue = init;
    F.Solver = 'cvodesstiff';
    F.RelativeTolerance = rtol;
    F.AbsoluteTolerance = atol;
    if ~isempty(Sp)
        if isa(Sp,'function_handle')
            F.Jacobian = odeJacobian('Jacobian', Sp);
        else
            F.Jacobian = odeJacobian('SparsityPattern', Sp);
        end
    end
    if divisions ~= 0
        F.SolverOptions.MaxStep = tspan(end)/divisions;
    end
    if numel(tspan) == 2
        sol = solve(F, tspan(1), tspan(2));
    else
        sol = solve(F, tspan);
    end
    t = sol.Time(:);
    X = sol.Solution.';
    
end
