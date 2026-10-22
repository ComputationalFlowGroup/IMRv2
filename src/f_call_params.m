% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Compatibility adapter for existing scripts. Prefer imr.configuration.legacy_options.
% The implementation lives only in the imr namespace.
function varargout = f_call_params(varargin)
    [varargout{1:nargout}] = imr.configuration.legacy_options(varargin{:});
end
