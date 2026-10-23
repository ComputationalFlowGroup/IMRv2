% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Generic fluid/default checks, independent of any optional material profile.
clear;
[water, water_ref] = imr.nondimensionalize();
assert(water_ref.rho8 == 1064);
assert(abs(f_pair(water,'cstar')*water_ref.Uc - 1484) < 1e-12);


% Explicit water selection must reproduce the no-profile reference set.
[explicit, explicit_ref] = imr.nondimensionalize('fluid',@imr.fluids.water);
assert(isequaln(explicit,water));
assert(explicit_ref.rho8 == water_ref.rho8 && explicit_ref.t0 == water_ref.t0);
si = {'rho8',998,'c8',1490,'mu',1e-3,'surft',0.072,'cp',4200, ...
    't8',300,'pv',1000,'tfin',1e-5,'g',0,'wave_type',1,'pa',100};
[legacy_override, legacy_ref] = imr.nondimensionalize(si{:});
[profile_override, profile_ref] = imr.nondimensionalize('fluid',@imr.fluids.water,si{:});
assert(isequaln(legacy_override,profile_override));
assert(legacy_ref.rho8 == profile_ref.rho8);
assert(profile_ref.fluid.SI.rho8 == 998);
ws = imr.fluid_properties(@imr.fluids.water,[],[290 298.15 310]);
expected = 1.17e11*exp(-5200./[290 298.15 310]);
assert(max(abs(ws.pressure./expected - 1)) < 1e-14);
assert(max(abs(ws.dpressure_dT./(expected.*5200./[290 298.15 310].^2) - 1)) < 1e-14);
assert(f_pair(explicit,'thetav') == 5200/298.15);
% A new material using existing model families needs only a data profile.
custom = imr.fluids.water();
custom.name = 'synthetic fluid';
custom.liquid.rho = 2000;
custom.liquid.c = 900;
custom.saturation.c1 = 2e10;
[custom_args, custom_ref] = imr.nondimensionalize('fluid',custom,'mu',0.002);
assert(custom_ref.rho8 == 2000);
assert(abs(f_pair(custom_args,'cstar')*custom_ref.Uc - 900) < 1e-12);
assert(custom_ref.fluid.SI.mu == 0.002);
f_expect_error(@() imr.nondimensionalize('fluid',custom,'re',10), 'IMR:fluidOverride');
custom.saturation.model = 'table_logp_inverseT';
f_expect_error(@() imr.nondimensionalize('fluid',custom), 'IMR:fluidProfile');
fprintf('Water compatibility, saturation derivatives and custom profile PASSED.\n');

% User-supplied saturation tables may be rows or columns. Check against
% the analytic derivative of log(P) linear in inverse temperature.
table = imr.fluids.water();
table.saturation.model = 'table_logp_inverseT';
table.saturation.temperature_range = [280 340];
table.saturation.T = [280 300 340];
table.saturation.P = [100 1000 10000];
query = [280 290 300 320 340];
row = imr.fluid_properties(table,[],query);
table.saturation.T = table.saturation.T';
table.saturation.P = table.saturation.P';
column = imr.fluid_properties(table,[],query');
assert(isequal(row.pressure,column.pressure'));
assert(isequal(row.dpressure_dT,column.dpressure_dT'));
for interval = 1:2
    T1 = table.saturation.T(interval);
    T2 = table.saturation.T(interval+1);
    logP1 = log(table.saturation.P(interval));
    logP2 = log(table.saturation.P(interval+1));
    slope = (logP2-logP1)/(1/T2-1/T1);
    T = (T1+T2)/2;
    pressure = exp(logP1+slope*(1/T-1/T1));
    actual = imr.fluid_properties(table,[],T);
    assert(abs(actual.pressure/pressure-1) < 1e-14);
    assert(abs(actual.dpressure_dT/(-slope*pressure/T^2)-1) < 1e-14);
    at_knot = imr.fluid_properties(table,[],T1);
    assert(abs(at_knot.dpressure_dT/(-slope*at_knot.pressure/T1^2)-1) < 1e-14);
end
fprintf('Row/column saturation tables and analytic derivatives PASSED.\n');

% Independently frozen pre-refactor numbers detect accidental default drift.
d = imr.configuration.defaults();
expected = imrtest.default_physics();
fields = fieldnames(expected);
for k = 1:numel(fields)
    name = fields{k};
    assert(abs(d.(name)-expected.(name)) <= 5e-13*max(1,abs(expected.(name))));
end
% Every converted physical group must agree with the defaults exposed publicly.
public = imr.default_inputs();
for k = 1:2:numel(water)
    name = water{k};
    assert(abs(public.(name)-water{k+1}) <= 5e-13*max(1,abs(water{k+1})));
end
assert(~isfield(imr.fluids.water(),'defaults') && ~isfield(imr.fluids.water(),'units'));
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
