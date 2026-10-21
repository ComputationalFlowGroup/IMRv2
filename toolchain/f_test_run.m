% file f_test_run.m
% brief contains function f_test_run

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function runs one test case from f_test_cases with the solver's
% console output hidden. It returns the radius history R, or [] and the
% error message if the run fails, so that one failing case does not stop a
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
        evalc(['[t, R] = solver(args{:})' char(59)]);
    catch err
        msg = err.message;
    end
    
end
