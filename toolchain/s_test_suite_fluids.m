% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Generic fluid/default checks, independent of any optional material profile.
clear;
addpath('../src','../tests/fixtures');
[water, water_ref] = f_dim2nondim();
assert(water_ref.rho8 == 1064);
assert(abs(f_pair(water,'cstar')*water_ref.Uc - 1484) < 1e-12);


% Explicit water selection must reproduce the no-profile reference set.
[explicit, explicit_ref] = f_dim2nondim('fluid',@fluids.water);
assert(isequaln(explicit,water));
assert(explicit_ref.rho8 == water_ref.rho8 && explicit_ref.t0 == water_ref.t0);
si = {'rho8',998,'c8',1490,'mu',1e-3,'surft',0.072,'cp',4200, ...
    't8',300,'pv',1000,'tfin',1e-5,'g',0,'wave_type',1,'pa',100};
[legacy_override, legacy_ref] = f_dim2nondim(si{:});
[profile_override, profile_ref] = f_dim2nondim('fluid',@fluids.water,si{:});
assert(isequaln(legacy_override,profile_override));
assert(legacy_ref.rho8 == profile_ref.rho8);
assert(profile_ref.fluid.SI.rho8 == 998);
ws = f_fluid_properties(@fluids.water,[],[290 298.15 310]);
expected = 1.17e11*exp(-5200./[290 298.15 310]);
assert(max(abs(ws.pressure./expected - 1)) < 1e-14);
assert(max(abs(ws.dpressure_dT./(expected.*5200./[290 298.15 310].^2) - 1)) < 1e-14);
assert(f_pair(explicit,'thetav') == 5200/298.15);
% A new material using existing model families needs only a data profile.
custom = fluids.water();
custom.name = 'synthetic fluid';
custom.liquid.rho = 2000;
custom.liquid.c = 900;
custom.saturation.c1 = 2e10;
[custom_args, custom_ref] = f_dim2nondim('fluid',custom,'mu',0.002);
assert(custom_ref.rho8 == 2000);
assert(abs(f_pair(custom_args,'cstar')*custom_ref.Uc - 900) < 1e-12);
assert(custom_ref.fluid.SI.mu == 0.002);
f_expect_error(@() f_dim2nondim('fluid',custom,'re',10), 'IMR:fluidOverride');
custom.saturation.model = 'table_logp_inverseT';
f_expect_error(@() f_dim2nondim('fluid',custom), 'IMR:fluidProfile');
fprintf('Water compatibility, saturation derivatives and custom profile PASSED.\n');

% Independently frozen pre-refactor numbers detect accidental default drift.
d = default_case();
expected = imrtest.default_physics();
fields = fieldnames(expected);
for k = 1:numel(fields)
    name = fields{k};
    assert(abs(d.(name)-expected.(name)) <= 5e-13*max(1,abs(expected.(name))));
end
% Every converted physical group must agree with the defaults exposed publicly.
public = f_imr_inputs();
for k = 1:2:numel(water)
    name = water{k};
    assert(abs(public.(name)-water{k+1}) <= 5e-13*max(1,abs(water{k+1})));
end
assert(~isfield(fluids.water(),'defaults') && ~isfield(fluids.water(),'units'));
fprintf('Default physical groups and single-source water data PASSED.\n');
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
