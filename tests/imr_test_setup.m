% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Add test suites, helpers, reference functions and source once per session.
% From tests/: imr_test_setup; s_test_suite_settings
% The public imr_setup deliberately does not add test-only paths.
function root = imr_test_setup()
    testRoot = fileparts(mfilename('fullpath'));
    root = fileparts(testRoot);
    addpath(testRoot, fullfile(testRoot,'suites'), ...
        fullfile(testRoot,'helpers'), fullfile(root,'src'));
end
