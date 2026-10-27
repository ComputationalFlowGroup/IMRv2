% file s_format.m
% brief contains script to run formatting throughout the toolchain and src
% folders

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This script formats the toolchain and src code
close;
clear;
clc;

rootPath = fileparts(fileparts(mfilename('fullpath')));
folders = fullfile(string(rootPath), ["toolchain","src","docs","benchmarks","tests"]);
fileList = cell(100,1);
count = 1;
thisScriptFile = strcat('s_format', '.m');

for folder = folders
    files = dir(fullfile(folder, '**', '*.m'));
    for f = files'
        filePath = fullfile(f.folder, f.name);

        [~, fileName, ext] = fileparts(filePath);
        if strcmp([fileName, ext], thisScriptFile)
            fprintf('Skipping self: %s\n', filePath);
            continue;
        end

        fileList{count} = filePath;
        count = count + 1;
    end
end

lastNonEmptyIdx = find(~cellfun(@isempty, fileList), 1, 'last');
fileList = fileList(1:lastNonEmptyIdx);

for i = 1:length(fileList)
    filePath = fileList{i};

    try
        fid = fopen(filePath, 'r');
        lines = cell(5000,1);
        count = 1;
        while ~feof(fid)
            line = fgetl(fid);
            lines{count} = line;
            count = count + 1;
        end
        fclose(fid);

        lastNonEmptyIdx = find(~cellfun(@isempty, lines), 1, 'last');
        lines = lines(1:lastNonEmptyIdx);

        % f_format_lines only re-indents: it never rewrites, splits or
        % searches a line's content, so a quoted string or a comment
        % (including one holding a ';') is never mangled.
        formattedLines = f_format_lines(lines);

        fid = fopen(filePath, 'w');
        fprintf(fid, '%s\n', formattedLines{:});
        fclose(fid);

        fprintf('Auto-indented: %s\n', filePath);
    catch e
        fprintf('Skipping: %s (Error: %s)\n', filePath, e.message);
    end
end

fprintf('Auto-indentation complete.\n');
