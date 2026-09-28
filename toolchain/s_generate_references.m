% file s_generate_references.m
% brief contains script to generate the stored results of the regression suite

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This script runs every regression case of f_test_cases and stores
% R(t) in tests/reference/<mode>.mat, one file per solver mode, with the
% case keys, for s_test_suite_regression. Regenerate the results only when
% a change of the results is intended and agreed, since the suite then
% compares against the new ones. Set overwrite = true to replace existing
% files. The runs use a parallel pool when the Parallel Computing Toolbox
% is available (parfor on the default pool) and run in series otherwise.
% Nothing is written unless every case succeeds: the files are written to
% a temporary folder and moved into place at the end. The script also
% reports how far the three solver modes are apart on the cases they
% share, as a sanity check on the stored results. Run it from toolchain/.
clc;
clear;

addpath('../src/');

% replace existing stored results
overwrite = false;
% cross-mode differences above this are listed for a look
report_level = 1e-2;
refdir = fullfile('..', 'tests', 'reference');

cases = f_test_cases('regression');
modes = unique({cases.mode});
for m = 1:numel(modes)
    if exist(fullfile(refdir, [modes{m} '.mat']), 'file') && ~overwrite
        error(['s_generate_references: stored results exist in %s. Set ' ...
            'overwrite = true to replace them.'], refdir);
    end
end

% the runs
n = numel(cases);
Rs = cell(1, n);
msgs = cell(1, n);
fprintf('Running %d cases...\n', n);
parfor k = 1:n
    [Rs{k}, msgs{k}] = f_test_run(cases(k));
end

bad = find(~cellfun(@isempty, msgs));
if ~isempty(bad)
    for k = bad
        fprintf('  %s %s: %s\n', cases(k).mode, cases(k).key, msgs{k});
    end
    error('s_generate_references: %d cases failed, nothing was written', ...
        numel(bad));
end

% how far the modes are apart on the cases they share
fprintf('\nLargest relative difference between modes, per shared case:\n');
nlisted = 0;
for i = 1:numel(modes)
    for j = i+1:numel(modes)
        ii = find(strcmp({cases.mode}, modes{i}));
        worst = 0;
        for k = ii
            jj = find(strcmp({cases.mode}, modes{j}) & ...
                strcmp({cases.key}, cases(k).key), 1);
            if isempty(jj)
                continue
            end
            d = max(abs(Rs{k}./Rs{jj} - 1));
            worst = max(worst, d);
            if d > report_level
                fprintf('  %s vs %s, %s: %.2e\n', modes{i}, modes{j}, ...
                    cases(k).key, d);
                nlisted = nlisted + 1;
            end
        end
        fprintf('%s vs %s: largest %.2e\n', modes{i}, modes{j}, worst);
    end
end
if nlisted > 0
    fprintf(['%d case(s) differ between modes by more than %g. Large ' ...
        'differences are expected for violent collapses at these grid ' ...
        'sizes, but check any that are not.\n'], nlisted, report_level);
end

% write all files to a temporary folder, then move them into place
tmp = tempname;
mkdir(tmp);
for m = 1:numel(modes)
    in = strcmp({cases.mode}, modes{m});
    keys = {cases(in).key};
    R = Rs(in);
    if exist('OCTAVE_VERSION', 'builtin')
        platform = 'GNU Octave';
    else
        platform = 'MATLAB';
    end
    info = sprintf('generated %s with %s %s', datestr(now, 31), ...
        platform, version);
    save(fullfile(tmp, [modes{m} '.mat']), 'keys', 'R', 'info', '-v7');
end
if ~exist(refdir, 'dir')
    mkdir(refdir);
end
for m = 1:numel(modes)
    movefile(fullfile(tmp, [modes{m} '.mat']), ...
        fullfile(refdir, [modes{m} '.mat']), 'f');
end
rmdir(tmp);
fprintf('\nWrote %d stored results for %d modes to %s\n', n, ...
    numel(modes), refdir);
