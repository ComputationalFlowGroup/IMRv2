% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Compatibility adapter for existing scripts. Prefer imr.sweep.
% The implementation lives only in the imr namespace.
function varargout = f_imr_sweep(varargin)
    [varargout{1:nargout}] = imr.sweep(varargin{:});
end
