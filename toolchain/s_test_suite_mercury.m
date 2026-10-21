% file s_test_suite_mercury.m
% brief checks fluid profiles, saturation, SI preparation and pressure KM

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Run from toolchain/. Uses base MATLAB and bundled .m references only.
clear;
addpath('../src','../tests/fixtures');
[grid, observations] = imrtest.mercury_reference();
p = fluids.mercury();
properties = @(P,T) f_fluid_properties(p,P,T);
saturation = @(T) f_fluid_properties(p,[],T);
prepare = @(varargin) f_dim2nondim('fluid',@fluids.mercury,varargin{:});
hg = properties(grid(:,1),grid(:,2));
actual = [hg.liquid.rho hg.liquid.c hg.liquid.h hg.vapor.v ...
    hg.vapor.c hg.vapor.h hg.liquid.K0 hg.liquid.K0_prime];
assert(max(abs(actual./grid(:,3:end) - 1),[],'all') < 5e-12);
observed = properties(observations(:,2),observations(:,1));
assert(max(abs(observed.liquid.rho./observations(:,4) - 1)) < 5e-12);
assert(max(abs(observed.liquid.c./observations(:,5) - 1)) < 5e-12);
rms_percent = 100*sqrt(mean((observed.liquid.c./observations(:,3) - 1).^2));
assert(abs(rms_percent - p.liquid.cs_rms_canonical_percent) < 1e-9);
fprintf('Mercury reference states and 115 observations PASSED, RMS %.9f%%.\n', rms_percent);

mesh = properties([0 1e8 1e9],[298.15 412.5]');
assert(isequal(size(mesh.liquid.rho),[2 3]));
assert(isequal(size(mesh.liquid.h),[2 3]));
for row = 1:2
    for col = 1:3
        scalar = properties(mesh.P(row,col),mesh.T(row,col));
        assert(abs(scalar.liquid.c/mesh.liquid.c(row,col) - 1) < 1e-14);
    end
end
assert(all(isinf(mesh.vapor.v(:,1))));
assert(all(mesh.vapor.rho(:,1) == 0));
assert(max(abs(mesh.liquid.rho(:,1) - mesh.liquid.rho0(:,1))) < 1e-10);
assert(max(abs(mesh.liquid.c(:,1) - mesh.liquid.c0(:,1))) < 1e-10);
at_ref = properties(101325,p.liquid.Tref);
assert(abs(at_ref.liquid.K0/p.liquid.K0_at_Tref - 1) < 1e-14);
assert(abs(at_ref.liquid.K0_prime/p.liquid.K0_prime_at_Tref - 1) < 1e-14);
f_expect_error(@() properties(-1,298.15), 'IMR:fluidPressure');
f_expect_error(@() properties(7e9,298.15), 'IMR:fluidPressure');
f_expect_error(@() properties(101325,274), 'IMR:fluidTemperature');
f_expect_failure(@() properties(NaN,298.15));
f_expect_failure(@() properties(101325,Inf));
f_expect_failure(@() properties(1i,298.15));
f_expect_failure(@() properties(101325,[]));
old_warning = warning('on','IMR:fluidExtrapolation');
restore_warning = onCleanup(@() warning(old_warning));
lastwarn('');
properties(101325,550);
[~, warning_id] = lastwarn;
assert(strcmp(warning_id,'IMR:fluidExtrapolation'));
clear restore_warning;
fprintf('Mercury vectorization and domain checks PASSED.\n');

sat_nodes = saturation(p.saturation.T);
Pv_nodes = sat_nodes.pressure;
assert(max(abs(Pv_nodes./p.saturation.P - 1)) < 1e-13);
middle_T = 2/(1/350 + 1/400);
assert(abs(f_saturation_pressure(p,middle_T)/sqrt(10.148*138.75) - 1) < 1e-13);
assert(all(diff(f_saturation_pressure(p,linspace(298.15,550,100))) > 0));
f_expect_error(@() saturation(298), 'IMR:fluidSaturation');
f_expect_error(@() saturation(551), 'IMR:fluidSaturation');
fprintf('Mercury saturation interpolation PASSED.\n');

base = {'stress',1,'g',0,'lambda1',0,'lambda2',0, ...
    'mu',0,'surft',0,'kappa',1.4,'r0',50e-6, ...
    'req',40e-6,'tvector',linspace(0,1e-6,41)};
[args, ref] = prepare(base{:});
ambient = properties(101325,298.15);
assert(abs(ref.rho8/13532.443297336 - 1) < 1e-12);
assert(abs(ref.Uc/sqrt(101325/ref.rho8) - 1) < 1e-14);
assert(abs(ref.t0/(ref.R0/ref.Uc) - 1) < 1e-14);
assert(f_pair(args,'radial') == 2 && f_pair(args,'kappa') == 1.4);
assert(f_pair(args,'pv') == 0);
assert(isinf(f_pair(args,'re')) && isinf(f_pair(args,'we')));
assert(isinf(f_pair(args,'ca')));
assert(abs(f_pair(args,'cstar')*ref.Uc/ambient.liquid.c - 1) < 1e-14);
[again, replay] = f_dim2nondim(ref.fluid.input_pairs{:});
assert(isequaln(args,again) && isequaln(ref,replay));
[wet, wet_ref] = prepare(base{:},'vapor',1);
assert(abs(f_pair(wet,'pv')*wet_ref.P8/0.26117 - 1) < 1e-13);
assert(abs(f_pair(wet,'rv')*wet_ref.P8/(wet_ref.rho8*wet_ref.T8)/p.vapor.Rs - 1) < 1e-13);
% Synthetic transport values test plumbing, not a physical Hg calibration.
transport = {'atg',1e-4,'btg',0.01,'atv',2e-4,'btv',0.02};
[warm, warm_ref] = prepare(base{:},transport{:}, ...
    'bubtherm',1,'medtherm',1,'km',5);
expected_dm = 5/(warm_ref.rho8*p.liquid.Cp);
assert(abs(f_pair(warm,'foh')*warm_ref.Uc*warm_ref.R0/expected_dm - 1) < 1e-13);
f_expect_error(@() prepare(base{:},transport{:},'bubtherm',1, ...
    'atg',0), 'IMR:fluidThermal');
f_expect_error(@() prepare(base{:},transport{:}, ...
    'bubtherm',1,'medtherm',1), 'IMR:fluidMissingInput');
f_expect_error(@() prepare('stress',1,'surft',0,'kappa',1.4), 'IMR:fluidMissingInput');
f_expect_error(@() prepare('mu',0,'kappa',1.4), 'IMR:fluidMissingInput');
f_expect_error(@() prepare('mu',0,'surft',0), 'IMR:fluidMissingInput');
f_expect_error(@() prepare(base{:},'c8',1484), 'IMR:fluidOverride');
f_expect_error(@() prepare(base{:},'rho8',1064), 'IMR:fluidOverride');
f_expect_error(@() prepare(base{:},'cstar',150), 'IMR:fluidOverride');
f_expect_error(@() prepare(base{:},'radial',5), 'IMR:fluidRadial');
f_expect_error(@() prepare(base{:},'masstrans',1), 'IMR:fluidMassTransfer');
f_expect_error(@() prepare(base{:},'collapse',1), 'IMR:fluidMassTransfer');
f_expect_error(@() prepare(base{:},'bubtherm',1), 'IMR:fluidMissingInput');
f_expect_error(@() prepare(base{:},'vapor',1,'t8',290), 'IMR:fluidSaturation');
f_expect_error(@() prepare(base{:},'vapor',1,'p8',0.1), 'IMR:fluidBoiling');
% Check table derivative away from knots against independent differences.
Tprobe = reshape([320 375 425 475],[2 2]);
st = saturation(Tprobe);
step = 1e-3;
fd = (f_saturation_pressure(p,Tprobe+step) - f_saturation_pressure(p,Tprobe-step))/(2*step);
assert(isequal(size(st.pressure),size(Tprobe)));
assert(max(abs(st.dpressure_dT./fd - 1),[],'all') < 1e-8);
% Compare with independently supplied ambient properties, without the helper.
manual = f_dim2nondim('radial',2,'stress',1,'g',0,'lambda1',0,'lambda2',0, ...
    base{:},'rho8',ambient.liquid.rho,'c8',ambient.liquid.c, ...
    'cp',p.liquid.Cp,'rv',p.vapor.Rs,'pv',0,'thetav',0, ...
    'lheat',ambient.latent_heat);
solvers = {@f_imr_nodal, @f_imr_nodal, @f_imr_spectral};
extra = {{'scheme',0}, {'scheme',1}, {}};
for mode = 1:numel(solvers)
    [time, radius] = solvers{mode}(args{:}, extra{mode}{:});
    [manual_time, manual_radius] = solvers{mode}(manual{:}, extra{mode}{:});
    assert(isequal(time, manual_time));
    assert(max(abs(radius - manual_radius)) < 1e-10);
    assert(all(isfinite(radius)) && numel(radius) == 41);
    equilibrium = prepare(base{:},'req',50e-6,'vapor',1);
    [~, equilibrium_radius] = solvers{mode}(equilibrium{:}, extra{mode}{:});
    assert(max(abs(equilibrium_radius - 1)) < 1e-10);
    fprintf('Mercury KM mode %d: explicit-input agreement and equilibrium PASSED.\n', mode);
end
fprintf('All mercury checks PASSED.\n');

function value = f_pair(args,name)
    k = find(strcmp(args(1:2:end),name),1,'last');
    value = args{2*k};
end

function f_expect_error(fun,id)
    try
        fun();
    catch err
        assert(strcmp(err.identifier,id), 'Unexpected error: %s', err.message);
        return
    end
    error('Expected error %s was not raised.', id);
end

function f_expect_failure(fun)
    try
        fun();
    catch
        return
    end
    error('Expected invalid property input to fail.');
end

function pressure = f_saturation_pressure(profile,T)
    sat = f_fluid_properties(profile,[],T);
    pressure = sat.pressure;
end
