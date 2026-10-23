% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Frozen pre-consolidation physical defaults. Do not compute from the converter.
function expected = default_physics()
    expected = struct();
    expected.Req = 0.2;
    expected.Rdotzero = 0;
    expected.Cstar = 152.071002414908;
    expected.GAMa = 3009.257340241796;
    expected.We = 36.1875;
    expected.Re8 = 623.3655014295973;
    expected.DRe = 0;
    expected.Ca = 101.325;
    expected.De = 0.01951719889306811;
    expected.LAM = 0.1;
    expected.v_lambda_star = 0;
    expected.Foh = 0.0002533864560818251;
    expected.Br = 7.639413243495116e-05;
    expected.alpha_g = 0.5761007932661034;
    expected.beta_g = 0.4263395891806357;
    expected.alpha_v = 0.3600629957913147;
    expected.beta_v = 0.6374966217619463;
    expected.chi = 0.1647900432365406;
    expected.Km_star = 20.12762009007293;
    expected.Fom = 0.04959728111106163;
    expected.L_heat_star = 23781.93575129533;
    expected.Rv_star = 1444.950104049612;
    expected.Ra_star = 898.6805465194676;
    expected.Pv_sat = 0.03076018438379325;
    expected.thetav = 17.44088546033876;
    expected.ee = 0;
    expected.om = 0;
    expected.tw = 0;
    expected.dt = 0;
end
