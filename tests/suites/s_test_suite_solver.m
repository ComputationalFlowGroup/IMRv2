% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Unified dispatch, replay and optional acceleration, without golden data.
clear;

% Selection follows ordinary struct/pair precedence; defaults stay neutral.
assert(~isfield(imr.default_inputs(),'solver'));
c = imr.configuration.prepare_case('');
assert(strcmp(c.solver,'nodal') && strcmp(c.inputs.solver,'nodal'));
p = imr.default_inputs();
p.solver = "SPECTRAL";
c = imr.configuration.prepare_case('',p);
assert(strcmp(c.solver,'spectral') && strcmp(c.inputs.solver,'spectral'));
c = imr.configuration.prepare_case('',p,'SOLVER','nodal');
assert(strcmp(c.solver,'nodal'));
c = imr.configuration.prepare_case('','solver','nodal','solver','spectral');
assert(strcmp(c.solver,'spectral'));
for invalid = {0,[],struct(),{'nodal'},'',"",'unknown',["nodal" "spectral"]}
    f_input_error(@() imr.solve('solver',invalid{1}));
end
f_input_error(@() imr.solve_nodal('solver','spectral'));
f_input_error(@() imr.solve_spectral('solver','nodal'));
% Incompatible features reach the same checker through the unified API.
f_input_error(@() imr.solve('solver','spectral','bubtherm',1,'vapor',1,'masstrans',1));
f_input_error(@() imr.solve('solver','spectral','jacobian',3));
f_input_error(@() imr.solve('solver','spectral','scheme',1));
f_input_error(@() imr.solve('solver','nodal','stress',5,'medtherm',1,'bubtherm',1));
try
    f_too_many_outputs();
    error('Expected an output-count error.');
catch err
    assert(strcmp(err.identifier,'IMR:outputs'));
end
fprintf('Unified solver selection and compatibility guards PASSED.\n');

% Explicit material conversion must preserve the new routing option.
[pairs,ref] = imr.nondimensionalize('fluid',@imr.fluids.water, ...
    'solver','spectral','req',50e-6,'r0',50e-6,'u0',0.1);
c = imr.configuration.prepare_case('',pairs{:});
assert(strcmp(c.solver,'spectral'));
assert(abs(c.initial.Rdot*ref.Uc-0.1) < 1e-14);
assert(abs(ref.Uc^2/ref.R0-ref.R0/ref.t0^2) < 1e-12*ref.Uc^2/ref.R0);

wrappers = {@imr.solve_nodal,@imr.solve_nodal,@imr.solve_spectral};
solvers = {'nodal','nodal','spectral'};
for mode = 1:3
    scheme = double(mode == 2);
    args = {'Req',1,'radial',2,'rdotzero',0.01,'stress',1,'ca',Inf, ...
        'tvector',linspace(0,0.03,7),'scheme',scheme};
    unified = cell(1,9);
    wrapped = cell(1,9);
    previous = cell(1,8);
    if mode == 1
        [unified{:}] = imr.solve(args{:});
    else
        [unified{:}] = imr.solve(args{:},'solver',solvers{mode});
    end
    [wrapped{:}] = wrappers{mode}(args{:});
    [previous{:}] = wrappers{mode}(args{:});
    assert(isequaln(unified,wrapped) && isequaln(unified(1:8),previous));
    assert(isequal(size(unified{9}),size(unified{1})));
    replayed = cell(1,9);
    [replayed{:}] = imr.solve(unified{8});
    assert(isequaln(unified,replayed));
end
fprintf('Wrapper equivalence, unchanged output positions and replay PASSED.\n');

% Independent RP and pressure-form KM equations at nonzero absolute times.
% Include surface tension and a Gaussian pulse to check both pressure terms.
times = linspace(0.2,0.3,121);
methods = [15 23];
if exist('ode','class') == 8 && ~isMATLABReleaseOlderThan('R2024a')
    methods = [methods 3];
end
for solver = {'nodal','spectral'}
    for radial = 1:2
        for method = methods
            args = {'solver',solver{1},'radial',radial,'method',method, ...
                'Req',1,'rdotzero',0.03,'stress',0,'vapor',0, ...
                'bubtherm',0,'medtherm',0,'masstrans',0,'we',10, ...
                'cstar',100,'kappa',1.4,'wave_type',1,'ee',0.03, ...
                'tw',0.02,'dt',0.25,'tvector',times, ...
                'reltol',1e-10,'abstol',1e-12};
            [t,R,U,P,~,~,~,used,A] = imr.solve(args{:});
            pulse = exp(-((t-0.25)/0.02).^2);
            pf = -0.03*pulse;
            pf_dot = 2*0.03*(t-0.25)/0.02^2.*pulse;
            pressure_rate = -3*1.4*P.*U./R;
            if radial == 1
                expected = (P-1-pf-0.1./R-1.5*U.^2)./R;
            else
                expected = ((1+U/100).*(P-1-pf-0.1./R) ...
                    + R/100.*(pressure_rate+0.1*U./R.^2-pf_dot) ...
                    - 1.5*(1-U/300).*U.^2)./((1-U/100).*R);
            end
            assert(t(1) == times(1) && t(end) == times(end));
            assert(max(abs(A-expected)) < 1e-10*max(1,norm(expected,Inf)));
            h = times(2)-times(1);
            sampled = (U(3:end)-U(1:end-2))/(2*h);
            assert(max(abs(sampled-A(2:end-1))) < 2e-4);
            assert(strcmp(used.solver,solver{1}));
        end
    end
    % Two-time input returns adaptive steps, also aligned with acceleration.
    [t,~,~,~,~,~,~,~,A] = imr.solve('solver',solver{1}, ...
        'Req',1,'tfin',0.03);
    assert(isequal(size(t),size(A)) && t(end) == 0.03);
end
fprintf('Independent radial equations, sampled velocity and adaptive times PASSED.\n');

% Every radial closure returns zero acceleration at equilibrium, including
% the enthalpy and Mie-Gruneisen models in both solver implementations.
for solver = {'nodal','spectral'}
    for radial = 1:6
        [~,~,~,~,~,~,~,~,A] = imr.solve('solver',solver{1}, ...
            'radial',radial,'Req',1,'tvector',[0 0.01 0.02]);
        assert(isreal(A) && all(isfinite(A)) && norm(A,Inf) < 1e-9);
    end
end

% Coupled wall reconstruction, vapor diffusion and local polymer stresses.
for mode = 1:4
    solver = 'nodal';
    scheme = double(mode == 2);
    mass = double(mode <= 2);
    stress = 3;
    if mode >= 3
        solver = 'spectral';
        stress = 3+2*(mode == 4);
    end
    args = {'solver',solver,'Req',1,'radial',2,'rdotzero',0.01, ...
        'bubtherm',1,'medtherm',1,'vapor',1,'masstrans',mass, ...
        'stress',stress,'spectral',double(mode == 4),'scheme',scheme, ...
        'Nt',8,'Mt',8,'Nv',6,'tvector',linspace(0,0.02,9)};
    current = cell(1,9);
    previous = cell(1,8);
    [current{:}] = imr.solve(args{:});
    [previous{:}] = imr.solve(args{:});
    assert(isequaln(current(1:8),previous));
    t = current{1};
    T = current{5};
    Tm = current{6};
    used = current{8};
    A = current{9};
    assert(isreal(A) && all(isfinite(A)) && isequal(size(A),size(t)));
    assert(all(isfinite(T(:))) && all(isfinite(Tm(:))));
    c = imr.configuration.prepare_case('',used,'progdisplay',1);
    sys = imr.solvers.build(c);
    % Quiet output evaluation uses the same physics without progress lines.
    captured = evalc(['sys.output_rhs(t(1),sys.initial)' char(59)]);
    assert(isempty(captured));
    rate = sys.output_rhs(t(1),sys.initial);
    assert(abs(A(1)-rate(c.layout.Rdot)) < 1e-9*max(1,abs(A(1))));
end
fprintf('Equilibrium, coupled fields and quiet acceleration evaluation PASSED.\n');

function f_input_error(fun)
    try
        fun();
    catch err
        assert(startsWith(err.message,'INPUT ERROR'), ...
            'Expected an input error, got: %s',err.message);
        return
    end
    error('Expected an input error.');
end

function f_too_many_outputs()
    out = cell(1,10);
    [out{:}] = imr.solve();
    assert(numel(out) == 10);
end
