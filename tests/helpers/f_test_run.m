% file f_test_run.m
% brief contains function f_test_run

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function runs one test case from f_test_cases with the solver's
% console output hidden. It returns the radius history R, or [] and the
% error message if the run fails or any physical output is nonfinite, so
% that one failing case does not stop a
% suite. Optional t is the actual returned time vector, including an early
% stop, so completion checks do not infer success from radius length alone.
function [R, msg, t] = f_test_run(c)
    
    R = [];
    t = [];
    msg = '';
    % solver and args are used inside the evalc string, which the code
    % analyzer does not read
    solver = c.solver; %#ok<NASGU>
    args = c.args; %#ok<NASGU>
    try
        % char(59) is the statement's semicolon, kept out of the string
        % literal because s_format splits lines with two semicolons
        out = cell(1,7);
        evalc(['[out{:}] = solver(args{:})' char(59)]);
        t = out{1};
        R = out{2};
        labels = {'time','radius','velocity','pressure','bubble temperature', ...
            'medium temperature','vapor fraction'};
        for k = 1:numel(out)
            v = out{k};
            if ~isreal(v) || any(~isfinite(v(:)))
                error('IMR:testOutput','Nonfinite or complex %s.',labels{k});
            end
        end
    catch err
        msg = err.message;
    end
    
end
