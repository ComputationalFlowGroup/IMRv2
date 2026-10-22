% file legacy_options.m
% brief contains function imr.configuration.legacy_options

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Legacy grouped-vector interface for dimensionless case preparation.
% Delegates parsing, defaults and validation to imr.configuration.prepare_case, then returns
% the original twelve outputs in their original order. Forward solvers use
% the named case directly. Input structs and name/value overrides retain
% their existing precedence; use imr.nondimensionalize before passing SI quantities.
function [eqns_opts, solve_opts, init_opts, init_stress, tspan_opts, out_opts, ...
    acos_opts, wave_opts, sigma_opts, thermal_opts, mass_opts, inputs] = ...
    legacy_options(solver, varargin)

c = imr.configuration.prepare_case(solver,varargin{:});
% equation settings
eqns_opts = [c.model.radial c.model.bubtherm c.model.medtherm c.model.stress ...
    c.model.eps3 c.model.masstrans];
% solver options
solve_opts = [c.numerics.method c.numerics.spectral c.numerics.divisions c.numerics.Nv ...
    c.numerics.Nt c.numerics.Mt c.numerics.Lv c.numerics.Lt ...
    c.numerics.rtol c.numerics.atol c.numerics.jacobian ...
    c.numerics.scheme];
% dimensionless initial conditions
init_opts = [c.initial.R c.initial.Rdot c.initial.P c.initial.Pv c.initial.Req c.stress.alphax];
% dimensionless initial stress
init_stress = c.initial.stress;
% time span options
tspan_opts = c.time;
% output options
out_opts = c.output.progress;

% physical parameters

% acoustic parameters
acos_opts = [c.acoustic.Cstar c.acoustic.GAMa c.acoustic.kappa c.acoustic.nstate c.acoustic.hugoniot_s];
% dimensionless waveform parameters
wave_opts = {[c.forcing.om c.forcing.ee c.forcing.tw c.forcing.dt c.forcing.mn c.forcing.wave_type], c.forcing.wave_poly, c.forcing.wave_dpoly};
% dimensionless viscoelastic
sigma_opts = [c.stress.We c.stress.Re8 c.stress.DRe c.stress.v_a ...
    c.stress.v_nc c.stress.Ca c.stress.LAM c.stress.De ...
    c.stress.JdotA c.stress.nu_model c.stress.v_lambda_star c.stress.iDRe];
% dimensionless thermal
thermal_opts = [c.thermal.Foh c.thermal.Br c.thermal.alpha_g c.thermal.beta_g ...
    c.thermal.alpha_v c.thermal.beta_v c.thermal.chi c.thermal.iota];
% dimensionless mass transfer
mass_opts = [c.mass.Fom c.mass.kv0 c.mass.Rv_star c.mass.Ra_star ...
    c.mass.L_heat_star c.mass.mv0 c.mass.ma0 c.mass.thetav];

inputs = c.inputs;
end
