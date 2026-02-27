% file s_compile_c_files.m
% brief contains script to recursively find and compile all .c files in project
clc;
close all;
clear;

rootDir = fileparts(pwd);
% recursively search through project for any c files
cfiles = dir(fullfile(rootDir, '**', '*.c'));
for k = 1:numel(cfiles)
    mexPath = fullfile(cfiles(k).folder, cfiles(k).name);
    try
        mex(mexPath);
        fprintf('Compiled: %s\n', mexPath);
    catch ME
        warning('Failed to compile %s: %s\n', mexPath, ME.message);
    end
end
