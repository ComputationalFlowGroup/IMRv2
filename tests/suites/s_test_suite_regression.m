% file s_test_suite_regression.m
% brief contains script to test the solvers against stored results

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This script runs the regression cases of f_test_cases in the three
% solver modes (imr.solve_nodal with finite differences and with Chebyshev
% collocation, and imr.solve_spectral) and compares R(t) with the stored
% results in tests/reference/<mode>.mat, written by s_generate_references.
% The cases cover every radial, thermal (with and without mass transfer)
% and stress model that the compatibility matrix allows, for two
% materials, and the non-Newtonian models, driving waveforms, the stress
% field representation and collapse. A case fails if its run fails, if it
% has no stored result, or if the L2 norm of its relative error over the
% output times exceeds the threshold. Run it from tests/ after imr_test_setup.
clc;
clear;


% largest allowed L2 norm of R./R_ref - 1 over the output times
threshold = 1e-4;
testRoot = fileparts(fileparts(mfilename('fullpath')));
refdir = fullfile(testRoot,'reference');

cases = f_test_cases('regression');
modes = unique({cases.mode});

% stored results, one file per mode
ref = struct();
for m = 1:numel(modes)
    file = fullfile(refdir, [modes{m} '.mat']);
    if exist(file, 'file')
        ref.(modes{m}) = load(file);
    else
        ref.(modes{m}) = struct('keys', {{}}, 'R', {{}});
        fprintf('No stored results for %s: %s is missing\n', modes{m}, file);
    end
end

failed = false(1, numel(cases));
msg = cell(1, numel(cases));
for k = 1:numel(cases)
    c = cases(k);
    j = find(strcmp(ref.(c.mode).keys, c.key), 1);
    if isempty(j)
        failed(k) = true;
        msg{k} = 'no stored result, run s_generate_references';
    else
        [R, msg{k}] = f_test_run(c);
        R_ref = ref.(c.mode).R{j};
        if ~isempty(msg{k})
            failed(k) = true;
        elseif ~isequal(size(R), size(R_ref))
            failed(k) = true;
            msg{k} = 'output size differs from the stored result';
        else
            err = norm(abs(R./R_ref - 1), 2);
            failed(k) = ~(err <= threshold);
            msg{k} = sprintf('L2 relative error = %.3e', err);
        end
    end
    fprintf('%-10s %-30s %s\n', c.mode, c.key, msg{k});
end

% stored results that no case uses any more
for m = 1:numel(modes)
    stale = setdiff(ref.(modes{m}).keys, {cases(strcmp({cases.mode}, ...
        modes{m})).key});
    if ~isempty(stale)
        fprintf('Note: %d stored %s results match no case, e.g. %s\n', ...
            numel(stale), modes{m}, stale{1});
    end
end

f_test_report('Regression', cases, failed, msg);
