% file s_test_suite_polytropic.m
% brief checks the constant-vapor polytropic pressure law

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Run from tests/ after imr_test_setup. Requires base MATLAB, no stored reference files.
clear;

R = [0.7 1 1.3];
U = [-0.2 0 0.3];
Pg0 = 0.8;
Pv = 0.2;
kappa = 1.4;
[P, Pdot] = imr.physics.polytropic_pressure(R, U, Pg0, Pv, kappa);
step = 1e-6;
Pplus = imr.physics.polytropic_pressure(R + step*U, U, Pg0, Pv, kappa);
Pminus = imr.physics.polytropic_pressure(R - step*U, U, Pg0, Pv, kappa);
assert(max(abs(Pdot - (Pplus - Pminus)/(2*step))) < 1e-8);
[Pgas, Dgas] = imr.physics.polytropic_pressure(R, U, Pg0, 0, kappa);
assert(max(abs(P - Pgas - Pv)) < 1e-14);
assert(isequal(Pdot, Dgas));
[~, Dvapor] = imr.physics.polytropic_pressure(R, U, 0, Pv, kappa);
assert(all(Dvapor == 0));
fprintf('Polytropic gas derivative and constant-vapor invariance PASSED.\n');

% Returned pressure is an integrated state. Check it against the algebraic
% pressure law during actual motion, including the spectral analytic path.
solvers = {@imr.solve_nodal, @imr.solve_nodal, @imr.solve_spectral};
extra = {{'scheme',0}, {'scheme',1}, {}};
tvec = linspace(0,0.2,41);
for mode = 1:numel(solvers)
    args = {'radial',2,'stress',0,'vapor',1,'pv',Pv,'req',0.8, ...
        'we',Inf,'tvector',tvec,'reltol',1e-10, ...
        'abstol',1e-10,'nt',8,'mt',8,'jacobian',2};
    [t, radius, ~, pressure] = solvers{mode}(args{:}, extra{mode}{:});
    % An early stop (NaN, event, integration failure) would otherwise show
    % up only as a confusing size mismatch below, not a clear diagnosis.
    [is_complete, completion_msg] = imrtest.completion_check(t, tvec(end));
    assert(is_complete, 'Mode %d did not complete: %s', mode, completion_msg);
    initial_gas = (1 - Pv)*0.8^(3*1.47);
    expected = initial_gas*radius.^(-3*1.47) + Pv;
    assert(max(abs(radius - radius(1))) > 1e-4);
    assert(max(abs(pressure - expected)) < 1e-7);
    fprintf('Mode %d polytropic pressure trajectory PASSED.\n', mode);
end
fprintf('All polytropic checks PASSED.\n');
