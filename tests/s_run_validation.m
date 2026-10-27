% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Edit these settings, then run this script from MATLAB or with -batch.
% Modes: tests, generate, generate-and-test. Tests never generate baselines.
mode = 'tests';
suites = 'all';
overwrite = false;

% The controller keeps its workspace separate from suites containing clear.
addpath(fileparts(mfilename('fullpath')));
validation_report = run_validation(mode,suites,overwrite);
