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
% Outputs: [t,R,Rdot,P,T,Tm,kv,inputs,Rddot,S,Sdot]; acceleration, the
% stress integral and its time derivative are evaluated only when
% requested (outputs 9-11). Use imr.solve to select the solver by name.
function varargout = solve_nodal(varargin)
    [varargout{1:nargout}] = imr.solvers.run('nodal',varargin{:});
end
