% file solve_spectral.m
% brief contains function imr.solve_spectral

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function features a Chebyshev spectral collocation solver of
% the PDEs involving thermal transport and viscoelasticity coupled to the
% radial (Rayleigh-Plesset family) dynamics
% Outputs: [t,R,Rdot,P,T,Tm,kv,inputs,Rddot,S,Sdot]; acceleration, the
% stress integral and its time derivative are evaluated only when
% requested (outputs 9-11). Use imr.solve to select the solver by name.
function varargout = solve_spectral(varargin)
    [varargout{1:nargout}] = imr.solvers.run('spectral',varargin{:});
end
