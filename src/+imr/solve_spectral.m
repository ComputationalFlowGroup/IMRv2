% file solve_spectral.m
% brief contains function imr.solve_spectral

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function features a Chebyshev spectral collocation solver of
% the PDEs involving thermal transport and viscoelasticity coupled to the
% radial (Rayleigh-Plesset family) dynamics
function varargout = solve_spectral(varargin)
    [varargout{1:nargout}] = imr.solvers.run('spectral',varargin{:});
end
