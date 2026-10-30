% file zener_field.m
% brief contains function imr.physics.zener_field

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function returns the time derivative of the local polymer
% stress difference sp = tau_p,rr - tau_p,thetatheta of the linear Maxwell,
% Jeffreys, Zener, 2nd-order Fung Zener and Mooney-Rivlin Zener models
% (stress 3, 4, 9, 11) at the points y = r/R, for the medium heating. The
% field obeys, at fixed y,
%   d(sp)/dt = adv - (sp - sKV)/De - 3*(Rdot/R)*(1 - y^-3)*sp
% where adv = -ydot*d(sp)/dy with ydot = (Rdot/R)*(y^-2 - y) is the
% advection computed by the caller on its own grid, and
%   sKV = se - 6*etap*a,   a = (Rdot/R)*y^-3,   etap = (1 - LAM)/Re8
% is the local Kelvin-Voigt target, with the hoop stretch lam = r/r0 and
% the elastic stress difference se: neo-Hookean (lam^-4 - lam^2)/Ca times
% 1 + alphax*(I1-3) for stress 4 (alphax = 0 for stress 3) or
% 1 + alphax*(I1-3) + (alphax^2/2)*(I1-3)^2 for the 2nd-order Fung spring
% (stress 9), and (1/Ca)*[(lam^-4-lam^2) - C2*(lam^4-lam^-2)] for the
% Mooney-Rivlin spring (stress 11). This law is the local field whose
% integral 2*int_1^Inf sp/y dy is exactly Z1/R^3 of imr.physics.material_stress, so it
% changes the distribution of the heating but not R(t). The far-field
% point (y = Inf), where sp = 0, is left to the caller. R and Rdot may be
% row vectors with sp holding one state per column.
function spdot = zener_field(sp, adv, y, R, Rdot, Req, Ca, De, Re8, ...
    LAM, stress, alphax, C2)

Rr = Rdot./R;
iy3 = y.^-3;
etap = (1 - LAM)/Re8;
% elastic stress difference from the hoop stretch, x2 = (r0/R)^2
x2 = (y.^3 - 1 + (Req./R).^3).^(2/3);
lam2 = y.^2./x2;
lam4i = x2.^2./y.^4;
if stress == 11
    lam4 = lam2.^2;
    lam2i = x2./y.^2;
    se = (lam4i - lam2)/Ca - (C2/Ca)*(lam4 - lam2i);
else
    % alphax only enters stress 4 (quadratic) and 9 (2nd-order Fung); it
    % is otherwise ignored, so stress 3 (plain linear Zener) is unaffected
    % by whatever alphax happens to be set to
    a4 = alphax*any(stress == [4 9]);
    x = lam4i + 2*lam2 - 3;
    if stress == 9
        se = (lam4i - lam2)/Ca.*(1 + a4*x + 0.5*a4.^2.*x.^2);
    else
        se = (lam4i - lam2)/Ca.*(1 + a4*x);
    end
end
% local Kelvin-Voigt target
sKV = se - 6*etap*Rr.*iy3;
spdot = adv - (sp - sKV)/De - 3*Rr.*(1 - iy3).*sp;

end
