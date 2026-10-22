% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Compatibility adapter for existing scripts. Prefer imr.solve_spectral.
% The implementation lives only in the imr namespace.
function varargout = f_imr_spectral(varargin)
    [varargout{1:nargout}] = imr.solve_spectral(varargin{:});
end
