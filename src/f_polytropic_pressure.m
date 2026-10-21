% file f_polytropic_pressure.m
% brief contains function f_polytropic_pressure

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief Dimensionless polytropic gas pressure plus constant vapor pressure.
% Pg0 is the gas contribution at R = 1. Inputs may be arrays of states.
% The derivative acts on the gas alone because Pv is constant.
function [P, Pdot] = f_polytropic_pressure(R, Rdot, Pg0, Pv, kappa)
    Pg = Pg0.*R.^(-3*kappa);
    P = Pg + Pv;
    Pdot = -3*kappa*Pg.*Rdot./R;
end
