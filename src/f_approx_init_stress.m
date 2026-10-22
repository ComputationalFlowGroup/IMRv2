% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Compatibility adapter for existing scripts. Prefer imr.approximate_initial_stress.
% The implementation lives only in the imr namespace.
function varargout = f_approx_init_stress(varargin)
    [varargout{1:nargout}] = imr.approximate_initial_stress(varargin{:});
end
