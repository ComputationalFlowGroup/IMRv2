% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Focused regressions for pressure forcing, time stepping and input limits.
% Run from toolchain/ in base MATLAB. No stored reference files are needed.
clear;
addpath('../src');
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
    accel = f_radial_eq(radial,1+forcing,0.2,forcing,0.2,0, ...
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
    actual = f_radial_eq(radial,Pwall,Pdot,Pfar-1,Fdot,0, ...
        1,U,0,0,Cstar,sam,no,gama,n,nog,s,0,0,0);
    assert(abs(actual-expected)/max(1,abs(expected)) < 1e-8);
end
fprintf('Driven enthalpy balance and independent quadrature PASSED.\n');

% A nonzero start time must not enlarge MaxStep. Constant solutions remove
% accuracy-driven step restrictions so the configured maximum is exercised.
methods = [15 23];
if exist('ode','class') == 8 && ~isMATLABReleaseOlderThan('R2024a')
    methods(end+1) = 3;
end
for method = methods
    [t, X] = f_odesolve(@(~,x) zeros(size(x)),1,method,40,[1000 1001]);
    assert(t(1) == 1000 && t(end) == 1001);
    assert(max(diff(t)) <= 1/40 + 1e-10);
    assert(max(abs(X-1)) < 1e-12);
end
fprintf('Nonzero-start maximum time step PASSED.\n');

% The implemented thermal transform needs positive gas slope. Preserve
% existing kappa=1 limits and allow unused zero vapor conductivity.
f_call_params('nodal','bubtherm',0,'kappa',1);
f_call_params('nodal','bubtherm',1,'kappa',1);
f_expect_failure(@() f_call_params('nodal','bubtherm',1,'alpha_g',0));
f_call_params('nodal','bubtherm',1,'alpha_v',0,'beta_v',0);
f_expect_failure(@() f_test_report('Empty',[],false(1,0),{}));
c = struct('solver',@f_early_return,'args',{{}});
[R, msg, t] = f_test_run(c);
assert(isempty(msg) && isequal(t,[0 0.5]') && isequal(R,[1 1]'));
fprintf('Thermal guards and test-runner time reporting PASSED.\n');

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

function [t,R] = f_early_return()
    t = [0;
    0.5];
    R = [1;
    1];
end
