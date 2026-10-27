% file s_lint.m
% brief contains script to run a lint checker throughout the code

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This script runs the Matlab linter to check for inefficiencies in
% the code, then a cheap self-check of the formatter (f_format_lines,
% called by s_format.m): checkcode lints a file as it already is, so it
% cannot by itself catch a formatter that corrupts code while reindenting
% it. The self-check instead runs the formatter on a small fixture that
% holds a semicolon inside a string and inside a comment, and fails the
% same way a lint issue does if either one is mangled.
clc;
clear;

% this script lives in toolchain/, next to src/
toolchainPath = fileparts(mfilename('fullpath'));
srcPath = fullfile(toolchainPath, '..', 'src');
rootPath = fileparts(toolchainPath);
testPath = fullfile(rootPath, 'tests');

files = [dir(fullfile(toolchainPath, '**', '*.m'));
dir(fullfile(srcPath, '**', '*.m'));
dir(fullfile(testPath, '**', '*.m'));
dir(fullfile(rootPath, 'imr_setup.m'))];

errorFound = false; % Track if any issues are detected

% A line must come out exactly as many lines as it went in, with its
% quoted and commented text untouched: only the leading whitespace may
% change. fixture holds one representative line of each kind that a
% naive, string-unaware formatter could split or corrupt on a ';'.
fixture = {
'function out = f_format_lines_fixture(x)'
'    message = ''left;right'';'
'    if x > 0'
'        disp(''a;b;c''); % a trailing comment; with its own semicolon'
'    end'
'    out = message;'
'end'
};
addpath(toolchainPath);
formatted = f_format_lines(fixture);
if numel(formatted) ~= numel(fixture)
    fprintf(['Formatter self-check FAILED: %d lines in, %d out. ' ...
        'f_format_lines must not split or merge a quoted string.\n'], ...
        numel(fixture), numel(formatted));
    errorFound = true;
else
    for k = 1:numel(fixture)
        if ~strcmp(strtrim(formatted{k}), strtrim(fixture{k}))
            fprintf(['Formatter self-check FAILED on line %d: got "%s", ' ...
                'expected "%s" (only its indentation may change).\n'], ...
                k, formatted{k}, strtrim(fixture{k}));
            errorFound = true;
        end
    end
end
if ~errorFound
    fprintf('Formatter self-check passed: strings and comments are preserved.\n');
end

for i = 1:length(files)
    filePath = fullfile(files(i).folder, files(i).name);
    fprintf('Checking: %s\n', filePath);
    
    issues = checkcode(filePath, '-id');
    
    if ~isempty(issues)
        fprintf('Linting issues found in: %s\n', filePath);
        for j = 1:length(issues)
            fprintf('%s (Line %d): %s\n', issues(j).id, issues(j).line, issues(j).message);
        end
        errorFound = true; % Mark that an issue was found
    end
end

% Exit MATLAB with an error code if any issues were detected
if errorFound
    fprintf('❌ Linting issues detected. Failing the workflow.\n');
    exit(1); % Fail the workflow
else
    fprintf('✅ No linting issues found. Passing the workflow.\n');
    exit(0); % Pass the workflow
end
