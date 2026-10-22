% file f_imr_spectral.m
% brief contains function f_imr_spectral

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function features a Chebyshev spectral collocation solver of
% the PDEs involving thermal transport and viscoelasticity coupled to the
% radial (Rayleigh-Plesset family) dynamics
function varargout = f_imr_spectral(varargin)
    prepared = f_prepare_case('spectral',varargin{:});
    system = f_imr_system(prepared);
    n = prepared.numerics;
    [t,X] = f_odesolve(system.rhs,system.initial,n.method,n.divisions, ...
        prepared.time,system.jacobian,n.rtol,n.atol);
    varargout = system.outputs(t,X);
end
