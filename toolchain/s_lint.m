% file s_lint.m
% brief contains script to run a lint checker throughout the code

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This script runs the Matlab linter to check for inefficiencies in
% the code
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
