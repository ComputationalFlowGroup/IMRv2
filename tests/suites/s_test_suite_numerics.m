% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Focused regressions for pressure forcing, time stepping and input limits.
% Run from tests/ after imr_test_setup in base MATLAB. No stored reference files are needed.
clear;
n = 7.15;
gama = 3009.26;
sam = 1 + gama;
no = (n-1)/n;
Cstar = sqrt(n*sam);
nog = (n-1)/2;
s = 1.65;

% Equal wall and far-field pressures, including equal time derivatives,
% must give zero acceleration at rest for every enthalpy model. The old
% pressure subtraction violates this identity whenever forcing is nonzero.
forcing = [-0.5 0 5 100];
for radial = 3:6
    accel = imr.physics.radial_acceleration(radial,1+forcing,0.2,forcing,0.2,0, ...
        1,0,0,0,Cstar,sam,no,gama,n,nog,s,0,0,0);
    assert(max(abs(accel)) < 1e-12);
    % Independently integrate dP/rho between the two pressures. Check both
    % H and its far-field derivative in a moving, driven state.
    Pwall = 500;
    Pfar = 101;
    U = 0.3;
    Pdot = -0.7;
    Fdot = 0.4;
    if radial <= 4
        rho = @(p) ((p+gama)/sam).^(1/n);
        local_C = sqrt(n*(Pwall+gama)/rho(Pwall));
    else
        rho = @(p) f_density(p,Cstar,nog,s);
        dp = 1e-2;
        local_C = sqrt(2*dp/(rho(Pwall+dp)-rho(Pwall-dp)));
    end
    C = Cstar;
    if mod(radial,2) == 0
        C = local_C;
    end
    H = integral(@(p) 1./rho(p),Pfar,Pwall,'AbsTol',1e-10,'RelTol',1e-12);
    Hdot = Pdot/rho(Pwall) - Fdot/rho(Pfar);
    expected = ((1+U/C)*H + Hdot/C - 1.5*(1-U/(3*C))*U^2)/(1-U/C);
    actual = imr.physics.radial_acceleration(radial,Pwall,Pdot,Pfar-1,Fdot,0, ...
        1,U,0,0,Cstar,sam,no,gama,n,nog,s,0,0,0);
    assert(abs(actual-expected)/max(1,abs(expected)) < 1e-8);
end
fprintf('Driven enthalpy balance and independent quadrature PASSED.\n');

% The Tait integral has a finite logarithmic limit at nstate = 1. Nearby
% exponents also need a stable evaluation rather than subtracting powers.
for exponent = [1 1-1e-8 1+1e-8 7.15]
    power = (exponent-1)/exponent;
    rho = @(p) ((p+gama)/sam).^(1/exponent);
    Pwall = 1 + 1e-6;
    Pfar = 1 - 1e-6;
    expected = integral(@(p) 1./rho(p),Pfar,Pwall,'AbsTol',1e-15);
    for radial = 3:4
        actual = imr.physics.radial_acceleration(radial,Pwall,0,Pfar-1,0,0, ...
            1,0,0,0,Cstar,sam,power,gama,exponent,nog,s,0,0,0);
        assert(isfinite(actual) && abs(actual-expected) < 1e-14);
        balanced = imr.physics.radial_acceleration(radial,[1 101],0,[0 100],0,0, ...
            1,0,0,0,Cstar,sam,power,gama,exponent,nog,s,0,0,0);
        assert(all(balanced == 0));
    end
end
fprintf('Tait logarithmic limit and near-ambient pressures PASSED.\n');

% Histotripsy supports zero and fractional pulse exponents. Edge values
% follow the explicit pointwise convention, with no 0*Inf derivatives.
om = 2;
shift = 3;
amplitude = -2;
edges = shift + [-pi/om pi/om];
times = [edges(1)-1 edges(1) shift edges(2) edges(2)+1];
for exponent = [0 0.25 0.5 1 4]
    args = {[om amplitude 0 shift exponent 2],[],[]};
    [pulse, slope] = imr.physics.far_field_pressure(times,args);
    assert(isequal(size(pulse),size(times)) && all(isfinite(slope)));
    assert(pulse(3) == amplitude && slope(3) == 0);
    assert(all(pulse([1 5]) == 0) && all(slope([1 2 4 5]) == 0));
    if exponent == 0
        assert(all(pulse([2 4]) == amplitude));
    else
        assert(all(pulse([2 4]) == 0));
    end
    interior = shift + [-1 -0.2 0.4 1];
    step = 1e-6;
    [~, slope] = imr.physics.far_field_pressure(interior,args);
    numerical = (imr.physics.far_field_pressure(interior+step,args)-imr.physics.far_field_pressure(interior-step,args))/(2*step);
    assert(max(abs(slope-numerical)) < 1e-8);
    [column, ~] = imr.physics.far_field_pressure(times',args);
    assert(isequal(column,pulse'));
    [matrix, ~] = imr.physics.far_field_pressure(reshape(interior,2,2),args);
    assert(isequal(size(matrix),[2 2]));
    args{1}(2) = 0;
    [pulse, slope] = imr.physics.far_field_pressure(times,args);
    assert(all(pulse == 0) && all(slope == 0));
end
% The old 1+cos expression rounds to zero at these interior points.
args = {[om amplitude 0 shift 0.25 2],[],[]};
[pulse, slope] = imr.physics.far_field_pressure(edges + [1e-10 -1e-10],args);
assert(all(isfinite(slope)) && all(pulse ~= 0));
fprintf('Histotripsy edges and independent derivative checks PASSED.\n');

% A nonzero start time must not enlarge MaxStep. Constant solutions remove
% accuracy-driven step restrictions so the configured maximum is exercised.
methods = [15 23];
if exist('ode','class') == 8 && ~isMATLABReleaseOlderThan('R2024a')
    methods(end+1) = 3;
end
for method = methods
    [t, X] = imr.numerics.integrate(@(~,x) zeros(size(x)),1,method,40,[1000 1001]);
    assert(t(1) == 1000 && t(end) == 1001);
    assert(max(diff(t)) <= 1/40 + 1e-10);
    assert(max(abs(X-1)) < 1e-12);
end
fprintf('Nonzero-start maximum time step PASSED.\n');

% The implemented thermal transform needs positive gas slope. Preserve
% existing kappa=1 limits and allow unused zero vapor conductivity.
imr.configuration.prepare_case('nodal','bubtherm',0,'kappa',1,'method',23);
imr.configuration.prepare_case('nodal','bubtherm',1,'kappa',1,'method',23);
f_expect_failure(@() imr.configuration.prepare_case('nodal','bubtherm',1,'alpha_g',0,'method',23));
imr.configuration.prepare_case('nodal','bubtherm',1,'alpha_v',0,'beta_v',0,'method',23);
f_expect_failure(@() f_test_report('Empty',[],false(1,0),{}));
c = struct('solver',@f_early_return,'args',{{}});
[R, msg, t] = f_test_run(c);
assert(isempty(msg) && isequal(t,[0 0.5]') && isequal(R,[1 1]'));
fprintf('Thermal guards and test-runner time reporting PASSED.\n');

% Invalid grids/time spans must fail during validation, before allocation
% or ODE setup. Real, increasing row and column time vectors remain valid.
for name = {'Nt','Mt','Nv','divisions'}
    for invalid = {Inf,NaN,1+1i,[1 2],-1,0.5}
        f_expect_failure(@() imr.configuration.prepare_case('nodal','method',23,name{1},invalid{1}));
    end
end
for invalid = {[0 1 2 3; 1 2 3 4], [0 1+1i], [0 Inf], [0 0], [1 0]}
    f_expect_failure(@() imr.configuration.prepare_case('nodal','method',23,'tvector',invalid{1}));
end
imr.configuration.prepare_case('nodal','tvector',[0 1 2],'method',23);
imr.configuration.prepare_case('nodal','tvector',[0 1 2]','method',23);
fprintf('Finite integer and real time-vector validation PASSED.\n');

% Disabled thermal fields allow small/zero grid counts. They must not
% trigger invalid stencils or zero-order collocation during initialization.
solvers = {@imr.solve_nodal,@imr.solve_spectral};
for k = 1:numel(solvers)
    solver = solvers{k};
    schemes = 0;
    if k == 1
        schemes = [0 1];
    end
    for scheme = schemes
        for unused = 0:2
            [t,R,U,P,T,Tm] = solver('bubtherm',0,'medtherm',0, ...
                'masstrans',0,'Req',1,'ee',0,'Nt',unused,'Mt',unused, ...
                'scheme',scheme,'tfin',0.01,'progdisplay',0,'method',23);
            assert(t(end) == 0.01 && all(isfinite(P)) && isempty(T));
            assert(max(abs(R-1)) < 1e-8 && max(abs(U)) < 1e-8);
            assert(all(Tm == 1));
        end
        [t,R] = solver('bubtherm',1,'medtherm',0,'masstrans',0, ...
            'Req',1,'ee',0,'Nt',5,'Mt',0,'scheme',scheme, ...
            'tfin',0.01,'progdisplay',0,'method',23);
        assert(t(end) == 0.01 && max(abs(R-1)) < 1e-8);
    end
end
fprintf('Inactive thermal grids in both solvers PASSED.\n');

function rho = f_density(p,C,nog,s)
    A = p/C^2;
    a = A*s^2 - nog;
    b = -2*A*s - 1;
    rho = 1 + 2*A./(-b+sqrt(b.^2-4*a.*A));
end

function f_expect_failure(fun)
    try
        fun();
    catch
        return
    end
    error('Expected invalid input to fail.');
end

function [t,R,U,P,T,Tm,kv] = f_early_return()
    t = [0;
    0.5];
    R = [1;
    1];
    U = [0 0]';
    P = R;
    T = [];
    Tm = R;
    kv = [];
end
