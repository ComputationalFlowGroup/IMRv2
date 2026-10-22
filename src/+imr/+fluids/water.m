% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Legacy reference properties, in SI. This is NOT a temperature-dependent
% water EoS calibration. Copy this data-only function to define a new fluid.
function p = water()
    p.name = 'water (legacy reference)';
    p.source = 'IMRv2 reference defaults, saturation: Preston (2004)';
    p.commit = '9ddbd70';
    p.temperature_range = [0 Inf];
    p.pressure_range = [0 Inf];
    p.sound_temperature_range = [0 Inf];
    p.liquid.model = 'constant';
    p.liquid.rho = 1064;
    p.liquid.c = 1484;
    p.liquid.Cp = 4181;
    p.liquid.q = 0;
    p.liquid.tait_pressure = 3049.13e5; % Pa, separate Tait branch parameter.
    p.vapor.model = 'ideal';
    p.vapor.Rs = 8.3144598/18.01528e-3;
    p.vapor.gamma = 1.33;
    % Approximate ideal-vapor diagnostic; not a fit.
    % Constant latent heat reproduces the historical converter exactly.
    p.caloric.model = 'constant_latent_heat';
    p.caloric.latent_heat = 2264.76e3;
    p.saturation.model = 'exp_inverse_temperature';
    p.saturation.c1 = 1.17e11;
    p.saturation.c2 = -5200;
    p.saturation.temperature_range = [0 Inf];
    p.capabilities.radial = 1:6;
    p.capabilities.mass_transfer = true;
    p.capabilities.collapse = true;
    p.derived_inputs = {};
    p.required = {};
    p.required_viscous = {};
    p.required_thermal = {};
    p.required_medium = {};
    % Reference transport data, not enthalpy coefficients.
    p.transport.surft = 0.07; % N/m, surface tension.
    p.transport.mu = 8.3283e-4; % Pa s, viscosity.
    p.transport.km = 0.55; % W/(m K), liquid conductivity.
    p.transport.dmass = 24.2e-6; % m^2/s, vapor mass diffusivity.
    p.transport.atv = 3.30e-5; % W/(m K^2), vapor conductivity slope.
    p.transport.btv = 1.742e-2; % W/(m K), vapor conductivity intercept.
end
