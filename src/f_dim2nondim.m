% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Compatibility adapter for existing scripts. Prefer imr.nondimensionalize.
% The implementation lives only in the imr namespace.
function varargout = f_dim2nondim(varargin)
    [varargout{1:nargout}] = imr.nondimensionalize(varargin{:});
end
