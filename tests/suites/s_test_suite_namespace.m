% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Namespace resolution, measured data and documented compatibility adapters.
clear;
testRoot = fileparts(fileparts(mfilename('fullpath')));
assert(strcmp(which('f_test_cases'),fullfile(testRoot,'helpers','f_test_cases.m')));
assert(strcmp(which('imrtest.default_physics'),fullfile(testRoot,'+imrtest','default_physics.m')));
assert(strcmp(which('s_test_suite_structure'),fullfile(testRoot,'suites','s_test_suite_structure.m')));
assert(isequaln(imr.default_inputs(),f_imr_inputs()));
assert(isequaln(imr.configuration.defaults(),default_case()));
assert(isequaln(imr.fluids.water(),fluids.water()));
[modern,modern_ref] = imr.nondimensionalize('fluid',@imr.fluids.water);
[legacy,legacy_ref] = f_dim2nondim('fluid',@fluids.water);
assert(isequaln(modern,legacy) && isequaln(modern_ref,legacy_ref));
assert(isequaln(imr.fluid_properties(@imr.fluids.water,[],300), ...
    f_fluid_properties(@fluids.water,[],300)));
for wave_type = -3:-1
    c = imr.configuration.prepare_case('nodal','wave_type',wave_type,'method',23);
    forcing = c.forcing;
    assert(~isempty(forcing.wave_poly) && ~isempty(forcing.wave_dpoly));
    args = {[forcing.om forcing.ee forcing.tw forcing.dt forcing.mn wave_type], ...
        forcing.wave_poly,forcing.wave_dpoly};
    [pressure,derivative] = imr.physics.far_field_pressure(0.1,args);
    assert(isfinite(pressure) && isfinite(derivative));
end
% New public entry points and old adapters must preserve every output.
modern_solvers = {@imr.solve_nodal,@imr.solve_spectral};
legacy_solvers = {@f_imr_nodal,@f_imr_spectral};
args = {'Req',1,'rdotzero',0.01,'stress',1,'ca',Inf, ...
    'bubtherm',0,'medtherm',0,'vapor',0,'masstrans',0, ...
    'tvector',linspace(0,0.1,6),'method',23};
for k = 1:numel(modern_solvers)
    modern = cell(1,8);
    legacy = cell(1,8);
    [modern{:}] = modern_solvers{k}(args{:});
    [legacy{:}] = legacy_solvers{k}(args{:});
    assert(isequaln(modern,legacy));
end
new_sweep = imr.sweep(@imr.solve_nodal,{'method',23},'req',[1 0.99]);
old_sweep = f_imr_sweep(@f_imr_nodal,{'method',23},'req',[1 0.99]);
assert(isequaln(new_sweep,old_sweep));
% A namespace resolver must not require package folders on the search path.
entries = strsplit(path,pathsep);
for k = 1:numel(entries)
    assert(isempty(regexp(entries{k},'[\\/]\+imr([\\/]|$)','once')));
end
fprintf('Namespace calls, measured data and compatibility adapters PASSED.\n');
