% file s_test_suite_mie_gruneisen.m
% brief checks imr.physics.mie_gruneisen_eos directly, independent of the
% bubble-dynamics ODE that is its only caller (imr.physics.radial_acceleration)

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Run from tests/ after imr_test_setup. Uses base MATLAB only.
clear;

% Parameter sets to check: water's own schema defaults (nstate = 7.15,
% hugoniot_s = 1.65, so nog = (nstate-1)/2 = 3.075), the s = 0 (purely
% linear/acoustic) limit recommended when no material-specific Hugoniot
% fit is available, and an arbitrary third point away from both.
param_sets = {
[1.0, 1.65, 3.075]   % Cstar, hugoniot_s, nog: water's schema defaults
[0.735, 0, 3.075]    % hugoniot_s = 0: linear/acoustic EoS, any fluid
[2.1, 0.9, 1.4]      % an arbitrary third point
};

for row = 1:numel(param_sets)
    p = param_sets{row};
    Cstar = p(1);
    hugoniot_s = p(2);
    nog = p(3);
    eos = @(pval) imr.physics.mie_gruneisen_eos(pval, Cstar, hugoniot_s, nog);

    % At the ambient reference state (p = 1), the EoS must reduce to the
    % ambient sound speed with zero enthalpy and unit specific volume,
    % for any hugoniot_s/nog: this is what makes coupling a new fluid
    % profile's ambient (rho8, c8) to this EoS always consistent at the
    % reference point (docs/model/fluids.md).
    [C1, hB1, hH1] = eos(1);
    assert(abs(C1/Cstar - 1) < 1e-13);
    assert(abs(hB1) < 1e-13);
    assert(abs(hH1 - 1) < 1e-13);

    % Thermodynamic identity: hB(p) = integral from 1 to p of hH(p') dp',
    % i.e. d(hB)/dp = hH (hH = 1/rho, hB = enthalpy). Check the closed
    % form against independent central differences over a pressure range
    % that stays on the physical (hH > 0) branch.
    pvals = linspace(0.3, 3, 25);
    [~, ~, hH] = eos(pvals);
    assert(all(hH > 0) && all(isreal(hH)));
    step = 1e-6;
    [~, hB_plus] = eos(pvals + step);
    [~, hB_minus] = eos(pvals - step);
    dhB_dp = (hB_plus - hB_minus)/(2*step);
    assert(max(abs(dhB_dp./hH - 1)) < 1e-6);

    % Vectorized and scalar evaluation must agree elementwise.
    [C_vec, hB_vec, hH_vec] = eos(pvals);
    for k = 1:numel(pvals)
        [Ck, hBk, hHk] = eos(pvals(k));
        assert(Ck == C_vec(k) && hBk == hB_vec(k) && hHk == hH_vec(k));
    end

    fprintf('Mie-Gruneisen EoS set %d (Cstar=%.3g, s=%.3g, nog=%.3g): reference state, thermodynamic identity, vectorization PASSED.\n', ...
        row, Cstar, hugoniot_s, nog);
end
fprintf('All Mie-Gruneisen EoS checks PASSED.\n');
