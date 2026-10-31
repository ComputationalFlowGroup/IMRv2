% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Independent contracts for preparation, field activation and Jacobians.
% Run from tests/ after imr_test_setup using base MATLAB, without trajectory baseline files.
clear;
schema = imr.configuration.schema();
public = imr.default_inputs();
d = imr.configuration.defaults();
assert(numel(unique({schema.name})) == numel(schema));
assert(numel(schema) == numel(fieldnames(public)));
for k = 1:numel(schema)
    item = schema(k);
    assert(isfield(d,item.field));
    assert(isequaln(public.(item.name),d.(item.field)));
end
% Parameter names in the guide must cover the executable schema. Defaults
% stay numeric in MATLAB, while the guide may round them for readability.
testRoot = fileparts(fileparts(mfilename('fullpath')));
guide = fileread(fullfile(testRoot,'..','docs','model','parameters.md'));
rows = regexp(guide,'(?m)^\|\s*([^|]+)\|','tokens');
tokens = regexp(strjoin([rows{:}],' '),'`([^`]+)`','tokens');
documented = cellfun(@(token) lower(token{1}),tokens,'UniformOutput',false);
assert(all(ismember({schema.name},documented)));
% Independently differentiate the shared inverse with respect to all three
% arguments; do not compare the helper with a second copy of its algebra.
for theta = [-0.01 0 0.02]
    a = 0.4;
    b = 0.7;
    [~,K,da,db] = imr.physics.thermal_inverse(theta,a,b);
    h = 1e-6;
    dt = (imr.physics.thermal_inverse(theta+h,a,b)-imr.physics.thermal_inverse(theta-h,a,b))/(2*h);
    na = (imr.physics.thermal_inverse(theta,a+h,b)-imr.physics.thermal_inverse(theta,a-h,b))/(2*h);
    nb = (imr.physics.thermal_inverse(theta,a,b+h)-imr.physics.thermal_inverse(theta,a,b-h))/(2*h);
    assert(max(abs([dt-1/K na-da nb-db])) < 1e-8);
end
p = public;
p.req = 1;
p.tvector = linspace(0,0.1,8);
p.tfin = p.tvector(end);
p.method = 'ode15s';
from_struct = imr.configuration.prepare_case('nodal',p,'ee',0.02);
from_pairs = imr.configuration.prepare_case('nodal','req',1,'tvector',p.tvector, ...
    'method','ode15s','ee',0.02);
assert(isequaln(from_struct,from_pairs));
assert(from_struct.numerics.method == 15);
% Compatibility serialization keeps the original output order and values.
[eq,~,initial,S,t,~,acoustic,wave] = imr.configuration.legacy_options('nodal',p,'ee',0.02);
assert(isequal(eq,[1 0 0 0 0 0]) && isempty(S));
assert(isequal(t,p.tvector) && initial(1) == 1 && initial(5) == 1);
assert(acoustic(3) == p.kappa && wave{1}(2) == 0.02);
fprintf('Public input schema, struct overrides and compatibility adapter PASSED.\n');
% An active stress expansion needs a mode; an inactive field may use zero.
try
    imr.configuration.prepare_case('spectral','stress',9,'spectral',1,'Nv',0);
    error('Expected active zero-mode stress to fail.');
catch err
    assert(strcmp(err.identifier,'IMR:stressGrid'));
end

checked = 0;
thermal_cases = [0 0 0;
1 0 0;
1 1 0;
1 0 1;
1 1 1];
for kind = {'nodal','spectral'}
    solver = kind{1};
    schemes = 0;
    thermal_count = 3;
    if strcmp(solver,'nodal')
        schemes = [0 1];
        thermal_count = 5;
    end
    for scheme = schemes
        for flags = thermal_cases(1:thermal_count,:)'
            for stress = [0 1 3 4 5 6 7 8 9 10 11]
                medium = flags(2);
                if stress == 9 && medium && strcmp(solver,'nodal')
                    continue
                end
                % PTT (10) and Giesekus (11) have no closed-form reduction
                % and run in imr.solve_spectral only
                if (stress == 10 || stress == 11) && strcmp(solver,'nodal')
                    continue
                end
                gies_args = {};
                if stress == 11
                    % a nonzero mobility factor, so the Jacobian check below
                    % exercises the Giesekus quadratic term, not its eps3 = 0
                    % limit (stress 9)
                    gies_args = {'eps3',0.3};
                end
                c = imr.configuration.prepare_case(solver,'Req',1,'radial',2, ...
                    'bubtherm',flags(1),'medtherm',medium, ...
                    'masstrans',flags(3),'vapor',flags(3), ...
                    'stress',stress,'spectral',double(strcmp(solver,'spectral')), ...
                    'Nv',6,'Nt',7,'Mt',7,gies_args{:}, ...
                    'scheme',scheme,'jacobian',2,'progdisplay',0);
                sys = imr.solvers.build(c);
                idx = c.layout;
                assert(numel(sys.initial) == idx.count);
                blocks = [idx.R idx.Rdot idx.P idx.bubble idx.medium ...
                    idx.vapor idx.stress idx.polymer];
                assert(isequal(blocks,1:idx.count));
                assert((sys.operatorSizes(1) > 0) == logical(flags(1)));
                assert((sys.operatorSizes(2) > 0) == logical(medium));
                state = sys.initial;
                state(idx.R) = 0.99;
                state(idx.Rdot) = -0.01;
                state(idx.P) = 1.01*state(idx.P);
                if flags(1)
                    if strcmp(solver,'nodal')
                        y = linspace(0,1,numel(idx.bubble))';
                        state(idx.bubble) = 1e-4*(1-y.^2);
                    else
                        state(idx.bubble(1:2)) = [1e-4 -1e-4]';
                    end
                end
                state(idx.stress) = state(idx.stress)+1e-5;
                time = 0.1;
                rhs = sys.rhs(time,state);
                assert(all(isfinite(rhs)) && numel(rhs) == idx.count);
                quiet_rhs = sys.output_rhs(time,state);
                assert(max(abs(quiet_rhs-rhs)) < 1e-9*max(1,norm(rhs,Inf)));
                if isa(sys.jacobian,'function_handle')
                    J = full(sys.jacobian(time,state));
                    independent = f_central(sys.rhs,time,state);
                    scale = max(1,max(abs(independent),[],2));
                    relative_error = max(abs(J-independent)./scale,[],'all');
                    assert(relative_error < 5e-4,'Jacobian mismatch %s, stress %d: %g',solver,stress,relative_error);
                    checked = checked+1;
                end
                if strcmp(solver,'nodal')
                    batch = sys.rhs(time,[state,state]);
                    assert(max(abs(batch-rhs),[],'all') < 1e-8*max(1,norm(rhs,Inf)));
                end
                out = sys.outputs([time;
                time+0.01],[state';
                state']);
                assert(numel(out) == 8 && isequal(size(out{2}),[2 1]));
                if medium
                    if strcmp(solver,'nodal')
                        assert(max(abs(out{5}(:,end)-out{6}(:,1))) < 1e-12);
                    else
                        assert(all(isfinite(out{6}), 'all'));
                    end
                end
                if flags(3)
                    Tw = out{5}(:,end);
                    pv = c.initial.Pv*exp(c.mass.thetav*(1-1./Tw));
                    expected = 1./(1+c.mass.Rv_star/c.mass.Ra_star*(out{4}./pv-1));
                    assert(max(abs(out{7}(:,end)-expected)) < 1e-10);
                end
            end
        end
    end
end
fprintf('State layouts, enabled operators, outputs and %d independent Jacobians PASSED.\n',checked);

% Both discretizations must reduce to the same pressure/stress ODE system
% when fields are disabled, including every radial closure.
for radial = 1:6
    options = {'radial',radial,'Req',1,'bubtherm',0,'medtherm',0, ...
        'masstrans',0,'stress',1,'Nt',0,'Mt',0,'jacobian',0};
    nodal = imr.solvers.build(imr.configuration.prepare_case('nodal',options{:}));
    spectral = imr.solvers.build(imr.configuration.prepare_case('spectral',options{:}));
    state = nodal.initial;
    state(1:2) = [0.98 -0.02]';
    a = nodal.rhs(0.2,state);
    b = spectral.rhs(0.2,state);
    assert(max(abs(a-b)) < 1e-12);
    assert(all(nodal.operatorSizes == 0) && all(spectral.operatorSizes == 0));
end
% Small thermal oscillations in the shared gas-only, cold-medium model.
% Different spatial bases should give converged radius and pressure traces.
time = linspace(0,0.1,20);
options = {'Req',1,'rdotzero',0.02,'bubtherm',1,'medtherm',0, ...
    'vapor',0,'masstrans',0,'stress',0,'tvector',time,'jacobian',2};
[tn,Rn,~,Pn] = imr.solve_nodal(options{:},'Nt',40,'Mt',0);
[ts,Rs,~,Ps] = imr.solve_spectral(options{:},'Nt',12,'Mt',0);
assert(tn(end) == time(end) && ts(end) == time(end));
assert(max(abs(Rn-Rs)) < 5e-4 && max(abs(Pn-Ps)) < 5e-4);
fprintf('Shared ODE physics and thermal nodal/spectral comparison PASSED.\n');

function J = f_central(rhs,t,x)
    J = zeros(numel(x),numel(x));
    for k = 1:numel(x)
        step = eps^(1/3)*max(1,abs(x(k)));
        plus = x;
        minus = x;
        plus(k) = plus(k)+step;
        minus(k) = minus(k)-step;
        J(:,k) = (rhs(t,plus)-rhs(t,minus))/(2*step);
    end
end
