% file f_wall_mass_bc.m
% brief contains function f_wall_mass_bc

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function evaluates the bubble-wall heat flux balance with
% vapor mass transfer and its gradient. The balance couples the wall
% temperature Tw, the saturated vapor concentration kvw(Tw, P), and the
% latent heat flux; through the saturated vapor pressure
% Pv(Tw) = C1*exp(C2/Tw) it is transcendental in Tw and has no closed-form
% root. f_imr_nodal solves the balance for the wall vapor mass fraction
% (f_wall_mass_closure); this theta_w form and its gradient give the
% implicit derivative of theta_w in f_jacobian_nodal.
% Inputs: v = [theta_w, Tm(iwall_m), theta(iwall_b), kv(end), kv(iwall_b),
% P], where iwall_m and iwall_b are the other nodes of the wall rows of the
% medium and bubble first-derivative matrices (two each for finite
% differences, every node for Chebyshev), and the constants c (wall-row
% coefficients grad_Tm, grad_Tr, grad_C with the wall node first, alpha_g,
% alpha_v, C1, C2, Rvg, Rva, Rg). Outputs: the residual r and its gradient
% dr with respect to v. The sums run in node order so that the finite
% difference result does not depend on the row length.
function [r, dr] = f_wall_mass_bc(v, c)
    
    nm = numel(c.grad_Tm) - 1;
    nb = numel(c.grad_Tr) - 1;
    im = 1 + (1:nm);
    ib = 1 + nm + (1:nb);
    ie = 2 + nm + nb;
    ik = ie + (1:nb);
    ip = ie + nb + 1;
    
    theta_w = v(1);
    kvE = v(ie);
    P = v(ip);
    
    % wall temperature from theta with the local mixture coefficient
    alpha_m = kvE*c.alpha_v + (1 - kvE)*c.alpha_g;
    u = sqrt(1 + 2*theta_w*alpha_m);
    Tw = (alpha_m - 1 + u)/alpha_m;
    % saturated vapor pressure and wall vapor concentration
    Pv = c.C1*exp(c.C2/Tw);
    dem = 1 + c.Rvg*(P/Pv - 1);
    kvw = 1/dem;
    % latent heat flux term G*H
    M = kvw*c.Rva + c.Rg;
    N = Tw*(1 - kvw);
    G = P/(M*N);
    H = c.grad_C(1)*kvw;
    for k = 1:nb
        H = H + c.grad_C(k+1)*v(ik(k));
    end
    
    r = c.grad_Tm(1)*Tw;
    for k = 1:nm
        r = r + c.grad_Tm(k+1)*v(im(k));
    end
    r = r + c.grad_Tr(1)*theta_w;
    for k = 1:nb
        r = r + c.grad_Tr(k+1)*v(ib(k));
    end
    r = r + G*H;
    
    if nargout > 1
        nv = numel(v);
        % d(Tw)/d(theta_w) and d(Tw)/d(kv(end)) through alpha_m
        dTw = zeros(1,nv);
        dTw(1) = 1/u;
        dTw(ie) = (1 + theta_w/u - Tw)/alpha_m*(c.alpha_v - c.alpha_g);
        % d(kvw) through Tw and P
        dPv_dTw = -Pv*c.C2/Tw^2;
        dkvw = kvw^2*c.Rvg*P/Pv^2*dPv_dTw*dTw;
        dkvw(ip) = dkvw(ip) - kvw^2*c.Rvg/Pv;
        % d(G) = G*(dP/P - dM/M - dN/N)
        dM = c.Rva*dkvw;
        dN = (1 - kvw)*dTw - Tw*dkvw;
        dG = -G*(dM/M + dN/N);
        dG(ip) = dG(ip) + G/P;
        dH = c.grad_C(1)*dkvw;
        dH(ik) = dH(ik) + c.grad_C(2:end);
        dr = c.grad_Tm(1)*dTw + dG*H + G*dH;
        dr(1) = dr(1) + c.grad_Tr(1);
        dr(im) = dr(im) + c.grad_Tm(2:end);
        dr(ib) = dr(ib) + c.grad_Tr(2:end);
    end
    
end
