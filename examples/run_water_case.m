% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% From the checkout root, run: imr_setup; run_water_case
% Pressure-form Keller-Miksis, polytropic gas and Newtonian water reference.
[args,ref] = imr.nondimensionalize('fluid',@imr.fluids.water, ...
    'r0',50e-6,'req',50e-6,'u0',0.1,'g',0, ...
    'radial',2,'stress',1,'bubtherm',0,'medtherm',0, ...
    'vapor',0,'masstrans',0,'tfin',20e-6);
[t,R,U,P] = imr.solve_nodal(args{:});
plot(t*ref.t0*1e6,R*ref.R0*1e6);
xlabel('Time (microseconds)');
ylabel('Radius (micrometers)');
