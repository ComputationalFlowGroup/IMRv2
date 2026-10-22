% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Compatibility adapter for existing scripts. Prefer imr.fluid_properties.
% The implementation lives only in the imr namespace.
function varargout = f_fluid_properties(varargin)
    [varargout{1:nargout}] = imr.fluid_properties(varargin{:});
end
