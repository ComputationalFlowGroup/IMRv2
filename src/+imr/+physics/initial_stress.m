% file initial_stress.m
% brief contains function imr.physics.initial_stress

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function computes the stress value at the bubble maximum
% radius by a shooting method on the initial wall velocity. The model does
% not contain heat or mass transfer. alphax is the strain-stiffening
% parameter of the quadratic neo-Hookean model (0 for neo-Hookean).
function [S0] = initial_stress(Req, Re, Ca, De, We, CL, Pv_star, alphax)
    if nargin < 8
        alphax = 0;
    end
    % the maximum radius grows with the initial wall velocity and equals
    % Req at rest, so bracket the root of Rmax - 1 by doubling. The
    % Keller-Miksis growth is singular at Rdot = CL, so the bracket stays
    % below the sound speed
    f_res = @(Rdoti) f_maxwell_growth_iter(Rdoti, Req, Re, Ca, De, We, CL) - 1;
    Rdot_max = (1 - 1e-6)*CL;
    lo = 0;
    hi = min(1, Rdot_max);
    while f_res(hi) < 0
        if hi == Rdot_max
            error(['imr.physics.initial_stress: no initial wall velocity below the ' ...
                'sound speed grows the bubble from Req = %g to R = 1; ' ...
                'give the initial stress with ''stress0'''], Req);
        end
        lo = hi;
        hi = min(2*hi, Rdot_max);
    end
    Rdotic = fzero(f_res, [lo, hi], optimset('TolX', 1e-10));
    tf_nd = 2;
    tspan = [0, tf_nd];
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
    X0 = [Req, Rdotic, 0];
    % root-find the maximum radius with an integrator event on Rdot = 0
    opts = odeset(opts,'Events',@(t,X) f_rmax_event(X));
    [~,X,~,Xe,~] = ode23tb(@(t,X) f_RP_bub_collapse_comp(X, Req, Re, Ca, ...
        De, We, CL), tspan, X0, opts);
    % the first zero may not be the global maximum, and the trajectory
    % itself is the fallback if no event fires
    Xall = vertcat(Xe, X);
    [~,Imax] = max(abs(Xall(:,1)));
    S0 = Xall(Imax,3);

    function Rmax = f_maxwell_growth_iter(Rdoti, Req, Re, Ca, De, We, CL)
        X0 = [Req, Rdoti, 0];
        % final time
        tf_nd = 2;
        tspan = [0, tf_nd];
        opts = odeset('RelTol',1e-6,'AbsTol',1e-6, ...
            'Events',@(t,X) f_rmax_event(X));
        [~,X,~,Xe,~] = ode23tb(@(t,X) f_RP_bub_collapse_comp(X, Req, Re, ...
            Ca, De, We, CL), tspan, X0, opts);
        % no event fires when the bubble never turns (Xe is then empty)
        Rmax = max(abs(X(:,1)));
        if ~isempty(Xe)
            Rmax = max(Rmax, max(abs(Xe(:,1))));
        end
    end

    % event locating the maximum radius, growth -> collapse crossing
    function [value,isterminal,direction] = f_rmax_event(X)
        % wall velocity
        value = X(2);
        % record and continue (first zero may not be global max)
        isterminal = 0;
        direction = -1;
    end

    % non-dimensional rayleigh plesset solver
    function [Xdot] = f_RP_bub_collapse_comp(X, Req, Re, Ca, De, We, CL)

        % bubble radius R(t)
        R = X(1);
        % bubble wall velocity Rdot(t)
        Rdot = X(2);
        % stress integral
        S = X(3);

        Pgo = 1-Pv_star+1/(We*Req);
        Pb = Pv_star+Pgo*(Req/R)^(3);
        Pbdot = -3*Pgo*(Req/R)^(3)*Rdot/R;
        if De == 0 || De == Inf
            S = (-0.5/Ca*(5-4*Req/R - (Req/R)^4) - 4*Rdot/(Re*R))*R^3;
            Sdot = (2*Rdot*((Req/R)^2 - (Req/R)^5)/Ca - 4/Re*Rdot^2/R^2)*R^3;
            Rddot = ((1+Rdot/CL)*(Pb - 1 -1/(We*R) + S/R^3) + R/CL*(Pbdot + 1/We*Rdot/R^2+ ...
                (Sdot-3*Rdot*R^2*S/R^3)/R^3) ...
                - 1.5*(1-Rdot/(3*CL))*Rdot^2)/((1-Rdot/CL)*R + 4/Re/CL);
        else % SLS or Maxwell
            % elastic stress of the (quadratic) neo-Hookean model and its
            % rate; alphax = 0 gives the neo-Hookean model
            Rst = Req/R;
            Se = (3*alphax - 1)*(5 - Rst^4 - 4*Rst)/(2*Ca) + ...
                (2*alphax/Ca)*(27/40 + Rst^8/8 + Rst^5/5 + Rst^2 - 2/Rst);
            dSe = (Rdot/R)*((3*alphax - 1)/(2*Ca))*(4*Rst^4 + 4*Rst) - ...
                (2*alphax/Ca)*(Rdot/R)*(Rst^8 + Rst^5 + 2*Rst^2 + 2/Rst);
            Sdot = -4*Rdot/(Re*R*De) + dSe - (S - Se)/De;
            Rddot = ((1+Rdot/CL)...
                *(Pb - 1 + S -1/(We*R)) + R/CL*(Pbdot + 1/We*Rdot/R^2 + Sdot) ...
                - 1.5*(1-Rdot/(3*CL))*Rdot^2)/((1-Rdot/CL)*R + 4/Re/CL);
        end
        Xdot = [Rdot;
        Rddot;
        Sdot];
    end

end
