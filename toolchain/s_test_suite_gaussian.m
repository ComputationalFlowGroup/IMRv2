% file s_test_suite_gaussian.m
% brief contains script to test the response to a Gaussian pressure pulse

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This script drives a bubble at equilibrium with a Gaussian pressure
% pulse and checks that it returns to its equilibrium radius, in the three
% solver modes: f_imr_nodal with finite differences and with Chebyshev
% collocation, and f_imr_spectral. It covers every radial, thermal
% (without mass transfer), vapor and viscous stress model that the
% compatibility matrix allows (f_test_cases). It needs no stored results.
% Run it from toolchain/.
clc;
clear;

addpath('../toolchain/');
addpath('../src/');

% largest allowed difference between the final and initial radius
threshold = 1e-2;

cases = f_test_cases('gaussian');
failed = false(1, numel(cases));
msg = cell(1, numel(cases));
for k = 1:numel(cases)
    [R, msg{k}] = f_test_run(cases(k));
    if isempty(msg{k})
        err = abs(R(end) - R(1));
        failed(k) = ~(err <= threshold);
        msg{k} = sprintf('|R(end) - R(1)| = %.3e', err);
    else
        failed(k) = true;
    end
    fprintf('%-10s %-28s %s\n', cases(k).mode, cases(k).key, msg{k});
end

f_test_report('Gaussian pulse', cases, failed, msg);
