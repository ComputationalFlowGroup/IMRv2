% file solve_nodal.m
% brief contains function imr.solve_nodal

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function features the nodal solver of the PDEs involving
% thermal transport, vapor mass transfer and viscoelasticity coupled to the
% radial (Rayleigh-Plesset family) dynamics. The fields are carried as node
% values, and the spatial discretization is selected by the differentiation
% matrices of imr.numerics.differentiation: scheme = 0 second-order finite differences,
% scheme = 1 Chebyshev collocation. Everything else is shared.
function varargout = solve_nodal(varargin)
    [varargout{1:nargout}] = imr.solvers.run('nodal',varargin{:});
end
