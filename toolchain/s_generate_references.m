% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Compatibility script for golden-data generation. Run from toolchain/.
% Set overwrite = true only to replace existing stored results deliberately.
% The implementation is shared with tests/run_validation.m.
overwrite = false;
rootPath = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootPath,'tests'));
imr_test_setup;
generate_references(overwrite);
