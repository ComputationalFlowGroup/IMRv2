% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Compatibility adapter for existing scripts. Prefer imr.solve_nodal.
% The implementation lives only in the imr namespace.
function varargout = f_imr_nodal(varargin)
    [varargout{1:nargout}] = imr.solve_nodal(varargin{:});
end
