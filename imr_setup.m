% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Run once per MATLAB session from the checkout root: imr_setup
% Adds the public source and examples, without recursively exposing test
% fixtures or package internals. Does not save or otherwise persist paths.
function root = imr_setup()
    root = fileparts(mfilename('fullpath'));
    addpath(root, fullfile(root,'src'), fullfile(root,'examples'));
end
