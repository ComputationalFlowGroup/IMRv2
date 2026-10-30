% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Explicit acceptance/rejection cases: never filter through the checker.
% Run from tests/ after imr_test_setup in base MATLAB, with no stored trajectory files.
clear;

% Malformed values must produce input errors before arithmetic/allocation.
bad = {
{'req',[1 2]}
{'req',0}
{'we',0}
{'kappa',1+1i}
{'eps3',0.1+1i}
{'eps3',0.6}
{'v_a',Inf,'nu_model',2,'stress',1}
{'v_a',[1 2],'nu_model',7,'stress',1}
{'p0star',NaN}
{'p0star',1+1i}
{'p0star',[1 2]}
{'p0star','x'}
{'stress',3,'stress0',NaN}
{'stress',3,'stress0',1+1i}
{'tfin','1'}
{'tfin',[1 2]}
{'method',["ode15s" "ode23tb"]}
{'bubtherm',1,'p0star',0,'vapor',0}
{'bubtherm',1,'vapor',1,'masstrans',1,'p0star',0}
{'bubtherm',1,'vapor',1,'rv',0}
{'bubtherm',1,'vapor',1,'ra',0}
{'bubtherm',1,'vapor',1,'p0star',0,'alpha_v',0,'beta_v',0}
};
for k = 1:numel(bad)
    f_input_error(@() imr.configuration.prepare_case('nodal',bad{k}{:}));
end
f_input_error(@() imr.configuration.prepare_case('spectral', ...
    'stress',5,'spectral',1,'Nv',2,'stress0',ones(2)));
p = imr.default_inputs();
f_input_error(@() imr.configuration.prepare_case('nodal',[p p]));
fprintf('Malformed and coupled physical input guards PASSED.\n');

% Integers are valid numeric inputs, but their arithmetic must use double.
c = imr.configuration.prepare_case('nodal','req',int32(1), ...
    'we',int32(10),'re',int32(10),'tvector',int32([0 1]));
assert(isa(c.initial.Pg,'double') && abs(c.initial.Pg-1.1) < 1e-14);
assert(isa(c.time,'double') && isa(c.stress.We,'double'));

% Independent forward map at known temperatures, including tiny and zero
% slopes: the previous inverse could give T=0 at theta=0 for positive a.
known_T = [0.5 1 2];
for a = [0 1e-20 1e-12 0.5]
    b = 1;
    theta = a*(known_T.^2-1)/2+b*(known_T-1);
    [T,K,Ta,Tb] = imr.physics.thermal_inverse(theta,a,b);
    assert(max(abs(T-known_T)) < 1e-14);
    assert(max(abs(K-(a*known_T+b))) < 1e-14);
    assert(all(isfinite(Ta)) && all(isfinite(Tb)));
    assert(abs(Ta(2))+abs(Tb(2)) == 0);
end
fprintf('Floating-point input normalization and thermal inverse limits PASSED.\n');

% Accepted boundaries are enumerated explicitly, rather than silently
% dropping any case that the checker rejects. All three discretizations,
% every compatible radial/thermal family and all supported Jacobians.
checked = 0;
thermal_cases = [0 0;
1 0;
1 1];
for mode = 1:3
    solver = 'nodal';
    scheme = mode-1;
    jacobians = 0:3;
    if mode == 3
        solver = 'spectral';
        scheme = 0;
        jacobians = 0:2;
    end
    for radial = 1:6
        thermal_count = 1+(radial <= 4)+(radial <= 2);
        for flags = thermal_cases(1:thermal_count,:)'
            for jacobian = jacobians
                args = {'Req',1,'radial',radial,'bubtherm',flags(1), ...
                    'medtherm',flags(2),'masstrans',0,'vapor',1,'pv',0, ...
                    'rv',0,'ra',0,'Nt',3*flags(1),'Mt',3*flags(2), ...
                    'scheme',scheme,'jacobian',jacobian,'kappa',1, ...
                    'alpha_g',1e-20,'beta_g',1,'alpha_v',0,'beta_v',0};
                c = imr.configuration.prepare_case(solver,args{:});
                sys = imr.solvers.build(c);
                f_finite_system(c,sys);
                assert(c.mass.kv0 == 0);
                assert(norm(sys.rhs(c.time(1),sys.initial),Inf) < 1e-8);
                checked = checked+1;
            end
        end
    end
    % Pure vapor without diffusion remains a valid frozen thermal mixture
    % when it has nonzero conductivity, including a constant conductivity.
    c = imr.configuration.prepare_case(solver,'Req',1,'bubtherm',1, ...
        'medtherm',0,'masstrans',0,'vapor',1,'p0star',0,'pv',1, ...
        'alpha_v',0,'beta_v',1,'Nt',3,'Mt',0,'scheme',scheme);
    sys = imr.solvers.build(c);
    f_finite_system(c,sys);
    assert(c.mass.kv0 == 1);
end
fprintf('%d explicitly accepted boundary configurations PASSED.\n',checked);

% Every non-Newtonian law at zero wall shear and with a disabled excess
% viscosity. Infinite exponents must not turn a zero contribution into NaN.
for nu = 1:7
    for stress = 1:5
        for excess = [0 10]
            c = imr.configuration.prepare_case('nodal','Req',1,'stress',stress, ...
                'nu_model',nu,'dre',excess,'v_lambda',1,'v_a',2,'v_nc',0.5, ...
                'bubtherm',0,'medtherm',0,'Nt',0,'Mt',0);
            f_finite_system(c,imr.solvers.build(c));
        end
    end
end
fprintf('Non-Newtonian zero-shear and zero-contribution combinations PASSED.\n');

% DRe=0 or a zero solvent share removes the excess viscosity exactly.
% Evaluating its unused law first can overflow and create 0*Inf=NaN.
for solver = {'nodal','spectral'}
    for stress = [1 3]
        args = {'Req',1,'stress',stress,'rdotzero',0.1,'bubtherm',0, ...
            'medtherm',0,'Nt',0,'Mt',0,'nu_model',1, ...
            'v_nc',1000,'v_lambda',1000,'dre',0};
        c = imr.configuration.prepare_case(solver{1},args{:});
        disabled = imr.solvers.build(c);
        c0 = imr.configuration.prepare_case(solver{1},args{:},'nu_model',0);
        reference = imr.solvers.build(c0);
        assert(isequal(disabled.rhs(0,disabled.initial),reference.rhs(0,reference.initial)));
        if stress == 3
            c = imr.configuration.prepare_case(solver{1},args{:},'dre',10,'lam',0);
            f_finite_system(c,imr.solvers.build(c));
        end
    end
end
fprintf('Inactive viscosity avoids overflow PASSED.\n');

% stress = 6 (PTT) and 7 (Giesekus) have no closed-form reduction: rejected
% in imr.solve_nodal, and forced to the spectral stress field in
% imr.solve_spectral (a warning, since spectral is left at its default 0).
f_input_error(@() imr.configuration.prepare_case('nodal','stress',6));
f_input_error(@() imr.configuration.prepare_case('nodal','stress',7));
warning('off','IMR:parameters');
c6 = imr.configuration.prepare_case('spectral','Req',1,'stress',6, ...
    'Nv',6,'Nt',0,'Mt',0);
warning('on','IMR:parameters');
assert(c6.numerics.spectral == 1);
f_finite_system(c6,imr.solvers.build(c6));
fprintf('PTT/Giesekus nodal rejection and forced spectral field PASSED.\n');

% Giesekus (stress 7) at eps3 = 0 must reduce exactly to UCM/Oldroyd-B
% (stress 5): the same right-hand side and Jacobian, not just a similar
% one, since the exponential and quadratic terms both vanish identically.
c5 = imr.configuration.prepare_case('spectral','Req',1,'stress',5, ...
    'spectral',1,'Nv',6,'Nt',0,'Mt',0,'jacobian',2);
warning('off','IMR:parameters');
c7 = imr.configuration.prepare_case('spectral','Req',1,'stress',7, ...
    'eps3',0,'Nv',6,'Nt',0,'Mt',0,'jacobian',2);
warning('on','IMR:parameters');
sys5 = imr.solvers.build(c5);
sys7 = imr.solvers.build(c7);
assert(isequal(sys5.initial,sys7.initial));
state = sys5.initial;
state(c5.layout.stress) = state(c5.layout.stress)+1e-5;
assert(isequal(sys5.rhs(0.1,state),sys7.rhs(0.1,state)));
assert(isequal(full(sys5.jacobian(0.1,state)),full(sys7.jacobian(0.1,state))));
fprintf('Giesekus eps3 = 0 reduces exactly to UCM/Oldroyd-B PASSED.\n');

% Scalar compatibility alone cannot establish that the assembled radial
% equation has a real, finite initial derivative. Check both builders.
for kind = {'nodal','spectral'}
    solver = kind{1};
    common = {'Req',1,'stress',0,'we',Inf,'Nt',0,'Mt',0};
    c = imr.configuration.prepare_case(solver,common{:}, ...
        'radial',2,'cstar',100,'rdotzero',100);
    f_initial_error(@() imr.solvers.build(c));
    c = imr.configuration.prepare_case(solver,common{:}, ...
        'radial',3,'gama',1,'ee',-3);
    f_initial_error(@() imr.solvers.build(c));
    c = imr.configuration.prepare_case(solver,common{:}, ...
        'radial',3,'gama',1,'we',0.1,'p0star',0.1);
    f_initial_error(@() imr.solvers.build(c));
end
fprintf('Initial sonic and EoS-domain guards PASSED.\n');

% Collapse initialization solves the positive equilibrium radius. Check
% the pressure balance independently, including boiling/zero-tension limits.
for pv = [0.1 1 1.1]
    for we = [10 Inf]
        args = {'collapse',1,'bubtherm',1,'medtherm',1,'masstrans',1, ...
            'vapor',1,'stress',0,'pv',pv,'we',we};
        target = 0.25;
        pg = (1-pv)*target^3+target^2/we;
        if pv >= 1 && isinf(we)
            f_input_error(@() imr.configuration.prepare_case('nodal',args{:},'p0star',0.01));
        else
            c = imr.configuration.prepare_case('nodal',args{:},'p0star',pg);
            assert(abs(c.initial.Req-target) < 1e-10);
            residual = pg/c.initial.Req^3+pv-1-1/(we*c.initial.Req);
            assert(abs(residual) < 1e-9);
        end
    end
end
fprintf('Bracketed positive collapse equilibria PASSED.\n');

% The shared test runner must notice failures outside the radius history.
c = struct('solver',@f_bad_pressure,'args',{{}});
[~,msg] = f_test_run(c);
assert(contains(msg,'pressure'));
fprintf('Full physical-output test-runner checks PASSED.\n');

function f_finite_system(c,sys)
    rate = sys.rhs(c.time(1),sys.initial);
    assert(numel(rate) == c.layout.count && isreal(rate) && all(isfinite(rate)));
    if isa(sys.jacobian,'function_handle')
        J = full(sys.jacobian(c.time(1),sys.initial));
        assert(isequal(size(J),[c.layout.count c.layout.count]));
        assert(isreal(J) && all(isfinite(J(:))));
    end
end

function f_input_error(fun)
    try
        fun();
    catch err
        assert(startsWith(err.message,'INPUT ERROR'), ...
            'Expected a deliberate input error, got: %s',err.message);
        return
    end
    error('Expected an input error.');
end

function f_initial_error(fun)
    try
        fun();
    catch err
        assert(strcmp(err.identifier,'IMR:initialState'), ...
            'Expected the assembled-state guard, got: %s',err.message);
        return
    end
    error('Expected an initial-state error.');
end

function [t,R,U,P,T,Tm,kv] = f_bad_pressure()
    t = [0 1]';
    R = [1 1]';
    U = [0 0]';
    P = [1 NaN]';
    T = [];
    Tm = R;
    kv = [];
end
