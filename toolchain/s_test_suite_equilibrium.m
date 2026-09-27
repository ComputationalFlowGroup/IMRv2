% file s_test_suite_equilibrium.m
% brief contains script to test that the equilibrium state is kept

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This script checks that a bubble started at equilibrium (R0 = Req,
% no driving) stays there, in the three solver modes: f_imr_nodal with
% finite differences and with Chebyshev collocation, and f_imr_spectral.
% It covers every radial, thermal (without mass transfer), vapor and stress
% model that the compatibility matrix allows (f_test_cases). It needs no
% stored results. Run it from toolchain/.
clc;
clear;

addpath('../toolchain/');
addpath('../src/');

% largest allowed deviation of R from 1: rounding only
threshold = 1e-13;

cases = f_test_cases('equilibrium');
failed = false(1, numel(cases));
msg = cell(1, numel(cases));
for k = 1:numel(cases)
    [R, msg{k}] = f_test_run(cases(k));
    if isempty(msg{k})
        err = max(abs(R - 1));
        failed(k) = ~(err <= threshold);
        msg{k} = sprintf('max |R - 1| = %.3e', err);
    else
        failed(k) = true;
    end
    fprintf('%-10s %-28s %s\n', cases(k).mode, cases(k).key, msg{k});
end

f_test_report('Equilibrium', cases, failed, msg);
