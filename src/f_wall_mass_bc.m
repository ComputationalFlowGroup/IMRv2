% file f_wall_mass_bc.m
% brief contains function f_wall_mass_bc

% brief This function evaluates the bubble-wall heat flux balance with
% vapor mass transfer and its gradient. The balance couples the wall
% temperature Tw, the saturated vapor concentration kvw(Tw, P), and the
% latent heat flux; through the saturated vapor pressure
% Pv(Tw) = C1*exp(C2/Tw) it is transcendental in Tw and has no closed-form
% root, so the solver finds theta_w by Newton iteration on this residual.
% Inputs: v = [theta_w, Tm(2), Tm(3), theta(end-1), theta(end-2), kv(end),
% kv(end-1), kv(end-2), P] and the constants c (grad_Tm, grad_Tr, grad_C
% stencil coefficients, alpha_g, alpha_v, C1, C2, Rvg, Rva, Rg). Outputs:
% the residual r and its gradient dr (1x9) with respect to v.
function [r, dr] = f_wall_mass_bc(v, c)
    
    theta_w = v(1);
    kvE = v(6);
    P = v(9);
    
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
    H = c.grad_C(1)*kvw + c.grad_C(2)*v(7) + c.grad_C(3)*v(8);
    
    r = c.grad_Tm(1)*Tw + c.grad_Tm(2)*v(2) + c.grad_Tm(3)*v(3) + ...
        c.grad_Tr(1)*theta_w + c.grad_Tr(2)*v(4) + c.grad_Tr(3)*v(5) + G*H;
    
    if nargout > 1
        % d(Tw)/d(theta_w) and d(Tw)/d(kv(end)) through alpha_m
        dTw = zeros(1,9);
        dTw(1) = 1/u;
        dTw(6) = (1 + theta_w/u - Tw)/alpha_m*(c.alpha_v - c.alpha_g);
        % d(kvw) through Tw and P
        dPv_dTw = -Pv*c.C2/Tw^2;
        dkvw = kvw^2*c.Rvg*P/Pv^2*dPv_dTw*dTw;
        dkvw(9) = dkvw(9) - kvw^2*c.Rvg/Pv;
        % d(G) = G*(dP/P - dM/M - dN/N)
        dM = c.Rva*dkvw;
        dN = (1 - kvw)*dTw - Tw*dkvw;
        dG = -G*(dM/M + dN/N);
        dG(9) = dG(9) + G/P;
        dH = c.grad_C(1)*dkvw;
        dH(7) = dH(7) + c.grad_C(2);
        dH(8) = dH(8) + c.grad_C(3);
        dr = c.grad_Tm(1)*dTw + dG*H + G*dH;
        dr(1) = dr(1) + c.grad_Tr(1);
        dr(2) = dr(2) + c.grad_Tm(2);
        dr(3) = dr(3) + c.grad_Tm(3);
        dr(4) = dr(4) + c.grad_Tr(2);
        dr(5) = dr(5) + c.grad_Tr(3);
    end
    
end
