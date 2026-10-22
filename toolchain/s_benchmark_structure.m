% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Local profiling aid. Run from toolchain/ in base MATLAB; compare the same
% machine, release and settings. Construction is a single cold measurement;
% timeit repeats warmed RHS/Jacobian evaluations. No speedup is asserted.
clear;
addpath('../src');
fprintf('MATLAB %s, %s\n',version,computer);
for solver = {'nodal','spectral'}
    started = tic;
    c = imr.configuration.prepare_case(solver{1},'Req',1,'radial',2, ...
        'bubtherm',1,'medtherm',1,'Nt',12,'Mt',12,'jacobian',2);
    sys = imr.solvers.build(c);
    construction = toc(started);
    state = sys.initial;
    state(2) = 0.01;
    rhs_time = timeit(@() sys.rhs(0,state));
    jacobian_time = timeit(@() sys.jacobian(0,state));
    fprintf('%s: %d states, construct %.6g s, RHS %.6g s, Jacobian %.6g s\n', ...
        solver{1},numel(state),construction,rhs_time,jacobian_time);
end
