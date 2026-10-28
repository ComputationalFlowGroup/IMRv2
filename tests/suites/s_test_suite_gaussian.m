% file s_test_suite_gaussian.m
% brief contains script to test the response to a Gaussian pressure pulse

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This script drives a bubble at equilibrium with a Gaussian pressure
% pulse and checks that it returns to its equilibrium radius, in the three
% solver modes: imr.solve_nodal with finite differences and with Chebyshev
% collocation, and imr.solve_spectral. It covers every radial, thermal
% (without mass transfer), vapor and viscous stress model that the
% compatibility matrix allows (f_test_cases). It needs no stored results.
% Run it from tests/ after imr_test_setup.
clc;
clear;


% largest allowed difference between the final and initial radius
threshold = 1e-2;

cases = f_test_cases('gaussian');
failed = false(1, numel(cases));
msg = cell(1, numel(cases));
for k = 1:numel(cases)
    [R, exec_msg, t] = f_test_run(cases(k));
    if isempty(exec_msg)
        % Check completion first
        t_request = cases(k).args{find(strcmpi(cases(k).args, 'tfin'), 1) + 1};
        [is_complete, completion_msg] = imrtest.completion_check(t, t_request);

        if ~is_complete
            failed(k) = true;
            msg{k} = completion_msg;
        else
            % Check physical condition (return to equilibrium)
            err = abs(R(end) - R(1));
            failed(k) = ~(err <= threshold);
            msg{k} = sprintf('%s | |ΔR| = %.3e', completion_msg, err);
        end
    else
        failed(k) = true;
        msg{k} = exec_msg;
    end
    fprintf('%-10s %-28s %s\n', cases(k).mode, cases(k).key, msg{k});
end

f_test_report('Gaussian pulse', cases, failed, msg);
