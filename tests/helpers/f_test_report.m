% file f_test_report.m
% brief contains function f_test_report

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function ends a test suite: it prints the number of cases per
% solver mode and every failure with its mode and case key, and raises an
% error if any case failed. In a batch run (matlab -batch, as in CI) the
% error exits MATLAB with a nonzero code, which fails the workflow;
% interactively it only stops the script. name is the suite name, cases
% the list from f_test_cases, and failed the logical vector of failures
% with their messages in msg.
function f_test_report(name, cases, failed, msg)
    if isempty(cases)
        error('IMR:emptyTestSuite', '%s: no accepted test cases were generated.', name);
    end
    
    modes = unique({cases.mode});
    for m = 1:numel(modes)
        in = strcmp({cases.mode}, modes{m});
        fprintf('%s, %s: %d of %d cases passed\n', name, modes{m}, ...
            sum(in & ~failed), sum(in));
    end
    if any(failed)
        fprintf('\nFailed cases:\n');
        for k = find(failed)
            fprintf('  %s %s: %s\n', cases(k).mode, cases(k).key, msg{k});
        end
        error('%s: %d of %d cases FAILED', name, sum(failed), numel(cases));
    end
    fprintf('All %d cases PASSED.\n', numel(cases));
    
end
