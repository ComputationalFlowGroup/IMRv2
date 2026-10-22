% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Compatibility adapter for existing scripts. Prefer imr.default_inputs.
% The implementation lives only in the imr namespace.
function varargout = f_imr_inputs(varargin)
    [varargout{1:nargout}] = imr.default_inputs(varargin{:});
end
