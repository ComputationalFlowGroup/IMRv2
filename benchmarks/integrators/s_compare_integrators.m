% file s_compare_integrators.m
% brief contains script to compare the time integrators and Jacobians

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This script compares the stiff time integrators of the nodal solver,
% ode23tb and ode15s (interpreted MATLAB) and CVODES (compiled SUNDIALS
% code through the ode object, MATLAB R2024a or later), each with the
% analytic Jacobian (jacobian = 2) and the batched finite-difference
% Jacobian (jacobian = 3). Every run is checked against a reference at
% tight tolerances on the same output times, and timed as the fastest of
% nrep runs after one warm-up run. The table gives the error in R, the
% wall time and the speedup over ode23tb with the analytic Jacobian (the
% default). Results are saved to integrators.mat and integrators.csv.
clc;
clear;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src'));

% settings
% repetitions per run, the fastest counts
nrep = 3;
% tolerances of the compared runs and of the reference
tol = 1e-8;
tol_ref = 1e-12;
% output times, the same for every run
tvector = linspace(0, 3.903439778613622, 2001);

% test cases: name and inputs (one per row)
cases = {
'polytropic, Rayleigh-Plesset', {'radial', 1}
'thermal, Keller-Miksis, Kelvin-Voigt', {'radial', 2, 'bubtherm', 1, ...
    'medtherm', 1, 'stress', 1}
'thermal, Keller-Miksis, Zener', {'radial', 2, 'bubtherm', 1, ...
    'medtherm', 1, 'stress', 3}
'thermal, Carreau viscosity', {'radial', 2, 'bubtherm', 1, ...
    'medtherm', 1, 'stress', 1, 'nu_model', 1, 'dre', 100, 'v_nc', 0.5, ...
    'v_lambda', 0.1}
'mass transfer, finite difference', {'radial', 2, 'bubtherm', 1, ...
    'medtherm', 1, 'stress', 1, 'vapor', 1, 'masstrans', 1}
'mass transfer, Chebyshev', {'radial', 2, 'bubtherm', 1, ...
    'medtherm', 1, 'stress', 1, 'vapor', 1, 'masstrans', 1, 'scheme', 1}
};

% integrators and Jacobians (one per row: label, method, jacobian)
combos = {
'ode23tb, analytic J', 'ode23tb', 2
'ode23tb, batched FD J', 'ode23tb', 3
'ode15s, analytic J', 'ode15s', 2
'ode15s, batched FD J', 'ode15s', 3
'cvodes, analytic J', 'cvodes', 2
'cvodes, batched FD J', 'cvodes', 3
};

% CVODES needs the ode object with the SUNDIALS solvers
has_cvodes = exist('ode', 'class') == 8;
if has_cvodes
    try
        F = ode;
        F.Solver = 'cvodesstiff';
    catch
        has_cvodes = false;
    end
end
if ~has_cvodes
    fprintf(['CVODES (the ode object with the cvodesstiff solver) is not ' ...
        'available in this MATLAB, its runs are skipped\n\n']);
end

ncase = size(cases, 1);
ncombo = size(combos, 1);
err = NaN(ncase, ncombo);
wall = NaN(ncase, ncombo);
msg = cell(ncase, ncombo);
common = {'tvector', tvector, 'reltol', tol, 'abstol', tol};
% one solver run on run_args, its console output hidden by evalc (char(59)
% is the statement's semicolon, kept out of the string because s_format
% splits lines with two semicolons)
run_cmd = ['[~, R] = imr.solve_nodal(run_args{:})' char(59)];

for ic = 1:ncase
    args = cases{ic, 2};
    fprintf('case %d: %s\n', ic, cases{ic, 1});
    
    % reference at tight tolerances, ode15s or else ode23tb
    R_ref = [];
    for ref_method = {'ode15s', 'ode23tb'}
        ref = {'tvector', tvector, 'reltol', tol_ref, 'abstol', tol_ref, ...
            'method', ref_method{1}, 'jacobian', 2};
        run_args = [args, ref];
        try
            evalc(run_cmd);
            R_ref = R;
            break
        catch
        end
    end
    if isempty(R_ref)
        fprintf('  no reference at tolerance %g, case skipped\n\n', tol_ref);
        msg(ic, :) = {'no reference'};
        continue
    end
    scale = max(abs(R_ref));
    
    for jc = 1:ncombo
        method = combos{jc, 2};
        if strcmp(method, 'cvodes') && ~has_cvodes
            msg{ic, jc} = 'skipped';
            continue
        end
        opts = [common, {'method', method, 'jacobian', combos{jc, 3}}];
        try
            % warm-up, then the fastest of nrep runs
            run_args = [args, opts];
            evalc(run_cmd);
            best = Inf;
            for k = 1:nrep
                t0 = tic;
                evalc(run_cmd);
                best = min(best, toc(t0));
            end
            wall(ic, jc) = best;
            err(ic, jc) = max(abs(R - R_ref))/scale;
            msg{ic, jc} = '';
        catch e
            msg{ic, jc} = e.message;
        end
        fprintf('  %-24s', combos{jc, 1});
        if isempty(msg{ic, jc})
            fprintf('error %8.1e   time %8.3f s   speedup %6.2f\n', ...
                err(ic, jc), wall(ic, jc), wall(ic, 1)/wall(ic, jc));
        else
            fprintf('failed: %s\n', msg{ic, jc});
        end
    end
    fprintf('\n');
end

% summary: geometric mean speedup over the default, and worst error
fprintf('%-24s %14s %12s\n', 'integrator', 'mean speedup', 'worst error');
for jc = 1:ncombo
    ok = isfinite(wall(:, jc)) & isfinite(wall(:, 1));
    if any(ok)
        sp = exp(mean(log(wall(ok, 1)./wall(ok, jc))));
        fprintf('%-24s %14.2f %12.1e\n', combos{jc, 1}, sp, max(err(ok, jc)));
    else
        fprintf('%-24s %14s %12s\n', combos{jc, 1}, '-', '-');
    end
end

% save
names = cases(:, 1);
labels = combos(:, 1);
save('integrators.mat', 'names', 'labels', 'err', 'wall', 'msg', 'tol', ...
    'tol_ref', 'nrep');
fid = fopen('integrators.csv', 'w');
fprintf(fid, 'case,integrator,error,time_s\n');
for ic = 1:ncase
    for jc = 1:ncombo
        fprintf(fid, '"%s","%s",%.3e,%.4f\n', names{ic}, labels{jc}, ...
            err(ic, jc), wall(ic, jc));
    end
end
fclose(fid);
