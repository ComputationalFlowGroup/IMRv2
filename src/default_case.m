% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Compatibility adapter for existing scripts. Prefer imr.configuration.defaults.
% The implementation lives only in the imr namespace.
function varargout = default_case(varargin)
    [varargout{1:nargout}] = imr.configuration.defaults(varargin{:});
end
