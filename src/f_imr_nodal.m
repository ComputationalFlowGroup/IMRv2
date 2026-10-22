% file f_imr_nodal.m
% brief contains function f_imr_nodal

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function features the nodal solver of the PDEs involving
% thermal transport, vapor mass transfer and viscoelasticity coupled to the
% radial (Rayleigh-Plesset family) dynamics. The fields are carried as node
% values, and the spatial discretization is selected by the differentiation
% matrices of f_diff_mat: scheme = 0 second-order finite differences,
% scheme = 1 Chebyshev collocation. Everything else is shared.
function varargout = f_imr_nodal(varargin)
    prepared = f_prepare_case('nodal',varargin{:});
    system = f_imr_system(prepared);
    n = prepared.numerics;
    [t,X] = f_odesolve(system.rhs,system.initial,n.method,n.divisions, ...
        prepared.time,system.jacobian,n.rtol,n.atol);
    varargout = system.outputs(t,X);
end
