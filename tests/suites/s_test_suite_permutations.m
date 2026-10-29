% file s_test_suite_permutations.m
% brief contains script to run every accepted combination of the solver
% options

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This script runs every combination of the solver options that the
% other suites do not cross, in the three solver modes: imr.solve_nodal with
% finite differences and with Chebyshev collocation, and imr.solve_spectral.
% It covers the non-Newtonian models with every viscous stress model, the
% driving waveforms with every radial model, the integrators and Jacobians
% with every thermal configuration, the limits Re, We, Ca = Inf, collapse
% with every stress model, and the remaining input forms (f_test_cases).
% A combination that imr.configuration.validate_case rejects is left out; every other one must
% run to the final time with finite results, so an input that passes the
% checks and then fails is caught here. The runs are short and on coarse
% grids, and the suite needs no stored results. Run it from tests/ after imr_test_setup.
clear;


cases = f_test_cases('permutations');
failed = false(1, numel(cases));
msg = cell(1, numel(cases));
for k = 1:numel(cases)
    [R, exec_msg, t] = f_test_run(cases(k));
    if isempty(exec_msg)
        % Expected output count: tvector's length if given (last occurrence
        % wins, matching how imr.nondimensionalize itself resolves a
        % repeated name), else unconstrained (a bare tfin has no fixed count).
        kt = find(strcmp(cases(k).args, 'tvector'), 1, 'last');
        if isempty(kt)
            nt = numel(R);
        else
            nt = numel(cases(k).args{kt + 1});
        end

        % Check completion
        [is_complete, completion_msg] = imrtest.completion_check(t, cases(k).t_request);

        if ~is_complete
            failed(k) = true;
            msg{k} = completion_msg;
        elseif numel(R) ~= nt || nt < 2
            failed(k) = true;
            msg{k} = sprintf('stopped after %d of %d output times', numel(R), nt);
        elseif ~all(isfinite(R))
            failed(k) = true;
            msg{k} = 'non-finite radius';
        else
            msg{k} = 'ran';
        end
    else
        failed(k) = true;
        msg{k} = exec_msg;
    end
    fprintf('%-10s %-28s %s\n', cases(k).mode, cases(k).key, msg{k});
end

f_test_report('Permutations', cases, failed, msg);
