% file s_run_mercury.m
% brief demonstrates mercury properties in pressure-form Keller-Miksis

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Run in base MATLAB. No Python, toolbox or mercury_eos checkout is needed.
% This example deliberately uses zero viscosity, elasticity and surface
% tension, not assumed transport coefficients for real mercury. Supply
% measured SI mu and surft for a physical material case.
clear;
% Run imr_setup once from the checkout root before this example.

[args, ref] = imr.nondimensionalize('fluid',@imr.fluids.mercury,'t8',298.15,'p8',101325, ...
    'radial',2,'stress',1,'g',0,'lambda1',0,'lambda2',0, ...
    'r0',50e-6,'req',40e-6,'mu',0,'surft',0,'kappa',1.4, ...
    'vapor',1,'tvector',linspace(0,12e-6,200));
[t, R] = imr.solve_nodal(args{:});
fprintf('Hg ambient density: %.6f kg/m^3\n', ref.rho8);
fprintf('Hg ambient sound speed: %.6f m/s\n', ref.fluid.properties.liquid.c);
fprintf('Hg ambient saturation pressure: %.6f Pa\n', ref.fluid.Pv);
plot(t*ref.t0*1e6, R*ref.R0*1e6);
xlabel('Time (microseconds)');
ylabel('Bubble radius (micrometers)');
