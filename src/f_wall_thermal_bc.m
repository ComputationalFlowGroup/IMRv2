% file f_wall_thermal_bc.m
% brief contains function f_wall_thermal_bc

% brief This function computes the bubble-wall value of the auxiliary
% temperature theta for the thermal-only (no mass transfer) heat flux
% balance c1*Tw + c2*theta_w + K = 0, where Tw = (ag-1+u)/ag and
% theta_w = (u^2-1)/(2*ag) with u = sqrt(1+2*ag*theta_w). The balance is
% quadratic in u, c2*u^2 + 2*c1*u + q = 0; the positive root (u = 1 at
% equilibrium) is taken in cancellation-free form. dtheta_dK is
% d(theta_w)/dK, where K collects the interior stencil contributions.
function [theta_w, dtheta_dK] = f_wall_thermal_bc(K, c1, c2, alpha_g)
    q = 2*c1*(alpha_g - 1) - c2 + 2*alpha_g*K;
    u = -q/(c1 + sqrt(c1^2 - c2*q));
    theta_w = (u^2 - 1)/(2*alpha_g);
    dtheta_dK = -u/(c2*u + c1);
end
