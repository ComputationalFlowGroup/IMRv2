% file f_zener_field.m
% brief contains function f_zener_field

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function returns the time derivative of the local polymer
% stress difference sp = tau_p,rr - tau_p,thetatheta of the linear Maxwell,
% Jeffreys and Zener models (stress 3 and 4) at the points y = r/R,
% for the medium heating. The field obeys, at fixed y,
%   d(sp)/dt = adv - (sp - sKV)/De - 3*(Rdot/R)*(1 - y^-3)*sp
% where adv = -ydot*d(sp)/dy with ydot = (Rdot/R)*(y^-2 - y) is the
% advection computed by the caller on its own grid, and
%   sKV = se - 6*etap*a,   a = (Rdot/R)*y^-3,   etap = (1 - LAM)/Re8
% is the local Kelvin-Voigt target, with the neo-Hookean elastic stress
% difference se = (lam^-4 - lam^2)/Ca of the hoop stretch lam = r/r0, times
% 1 + alphax*(I1 - 3) for the quadratic law (alphax = 0 for stress 3).
% This law is the local field whose integral 2*int_1^Inf sp/y dy is exactly
% Z1/R^3 of f_stress, so it changes the distribution of the heating but not
% R(t). The far-field point (y = Inf), where sp = 0, is left to the caller.
% R and Rdot may be row vectors with sp holding one state per column.
function spdot = f_zener_field(sp, adv, y, R, Rdot, Req, Ca, De, Re8, ...
    LAM, alphax)

Rr = Rdot./R;
iy3 = y.^-3;
etap = (1 - LAM)/Re8;
% elastic stress difference from the hoop stretch, x2 = (r0/R)^2
x2 = (y.^3 - 1 + (Req./R).^3).^(2/3);
lam2 = y.^2./x2;
lam4i = x2.^2./y.^4;
se = (lam4i - lam2)/Ca.*(1 + alphax*(lam4i + 2*lam2 - 3));
% local Kelvin-Voigt target
sKV = se - 6*etap*Rr.*iy3;
spdot = adv - (sp - sKV)/De - 3*Rr.*(1 - iy3).*sp;

end
