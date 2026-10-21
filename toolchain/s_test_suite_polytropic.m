% file s_test_suite_polytropic.m
% brief checks the constant-vapor polytropic pressure law

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Run from toolchain/. Requires base MATLAB, no stored reference files.
clear;
addpath('../src');

R = [0.7 1 1.3];
U = [-0.2 0 0.3];
Pg0 = 0.8;
Pv = 0.2;
kappa = 1.4;
[P, Pdot] = f_polytropic_pressure(R, U, Pg0, Pv, kappa);
step = 1e-6;
Pplus = f_polytropic_pressure(R + step*U, U, Pg0, Pv, kappa);
Pminus = f_polytropic_pressure(R - step*U, U, Pg0, Pv, kappa);
assert(max(abs(Pdot - (Pplus - Pminus)/(2*step))) < 1e-8);
[Pgas, Dgas] = f_polytropic_pressure(R, U, Pg0, 0, kappa);
assert(max(abs(P - Pgas - Pv)) < 1e-14);
assert(isequal(Pdot, Dgas));
[~, Dvapor] = f_polytropic_pressure(R, U, 0, Pv, kappa);
assert(all(Dvapor == 0));
fprintf('Polytropic gas derivative and constant-vapor invariance PASSED.\n');

% Returned pressure is an integrated state. Check it against the algebraic
% pressure law during actual motion, including the spectral analytic path.
solvers = {@f_imr_nodal, @f_imr_nodal, @f_imr_spectral};
extra = {{'scheme',0}, {'scheme',1}, {}};
for mode = 1:numel(solvers)
    args = {'radial',2,'stress',0,'vapor',1,'pv',Pv,'req',0.8, ...
        'we',Inf,'tvector',linspace(0,0.2,41),'reltol',1e-10, ...
        'abstol',1e-10,'nt',8,'mt',8,'jacobian',2};
    [~, radius, ~, pressure] = solvers{mode}(args{:}, extra{mode}{:});
    initial_gas = (1 - Pv)*0.8^(3*1.47);
    expected = initial_gas*radius.^(-3*1.47) + Pv;
    assert(max(abs(radius - radius(1))) > 1e-4);
    assert(max(abs(pressure - expected)) < 1e-7);
    fprintf('Mode %d polytropic pressure trajectory PASSED.\n', mode);
end
fprintf('All polytropic checks PASSED.\n');
