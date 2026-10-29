% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Generate keyed regression baselines after every case completes successfully.
% files = generate_references(overwrite), default false, writes tests/reference.
% Optional developer arguments refdir, cases and parallel allow isolated tests
% with synthetic cases and temporary files. Production callers use defaults.
% Case arrays follow f_test_cases and specify tvector or tfin in args.
function files = generate_references(overwrite,refdir,cases,parallel)
    if nargin < 1
        overwrite = false;
    end
    validateattributes(overwrite,{'logical'},{'scalar'},mfilename,'overwrite');
    if nargin < 2
        test_root = fileparts(fileparts(mfilename('fullpath')));
        refdir = fullfile(test_root,'reference');
    end
    valid_text = (ischar(refdir) && isrow(refdir)) || (isstring(refdir) && isscalar(refdir));
    if ~valid_text || strlength(string(refdir)) == 0
        error('IMR:referenceDirectory','refdir must be a nonempty directory path.');
    end
    refdir = char(refdir);
    if nargin < 4
        parallel = true;
    end
    validateattributes(parallel,{'logical'},{'scalar'},mfilename,'parallel');
    if nargin < 3
        cases = f_test_cases('regression');
    end
    f_check_cases(cases);
    % Cross-mode differences above this level are listed for review.
    report_level = 1e-2;

    modes = unique({cases.mode});
    files = cell(1,numel(modes));
    for m = 1:numel(modes)
        files{m} = fullfile(refdir,[modes{m} '.mat']);
    end

    f_check_overwrite(files,overwrite);

    % the runs
    n = numel(cases);
    Rs = cell(1, n);
    msgs = cell(1, n);
    fprintf('Running %d cases...\n', n);
    workers = 0;
    if parallel && license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
        workers = Inf;
    end
    parfor (k = 1:n, workers)
        [Rs{k}, msgs{k}] = f_complete_case(cases(k));
    end

    bad = find(~cellfun(@isempty, msgs));
    if ~isempty(bad)
        for k = bad
            fprintf('  %s %s: %s\n', cases(k).mode, cases(k).key, msgs{k});
        end
        error('IMR:referenceCases','%d cases failed, nothing was written', ...
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
    cleanup = onCleanup(@() f_remove_temporary(tmp));
    generated_at = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    info = sprintf('generated %s with MATLAB %s', generated_at, version);
    for m = 1:numel(modes)
        in = strcmp({cases.mode}, modes{m});
        keys = {cases(in).key};
        R = Rs(in);
        save(fullfile(tmp, [modes{m} '.mat']), 'keys', 'R', 'info', '-v7');
    end
    if ~exist(refdir, 'dir')
        mkdir(refdir);
    end
    % Recheck after the simulations, before publishing any staged results.
    f_check_overwrite(files,overwrite);
    for m = 1:numel(modes)
        movefile(fullfile(tmp, [modes{m} '.mat']), ...
            fullfile(refdir, [modes{m} '.mat']), 'f');
    end
    fprintf('\nWrote %d stored results for %d modes to %s\n', n, ...
        numel(modes), refdir);

end

function f_check_cases(cases)
    required = {'mode','key','solver','args'};
    if ~isstruct(cases) || isempty(cases) || ~all(isfield(cases,required))
        error('IMR:referenceCases','Expected a nonempty array of regression cases.');
    end
    for k = 1:numel(cases)
        c = cases(k);
        valid = ischar(c.mode) && isrow(c.mode) && ~isempty(regexp(c.mode,'^[a-zA-Z]\w*$','once')) && ...
            ischar(c.key) && isrow(c.key) && ~isempty(c.key) && ...
            isa(c.solver,'function_handle') && iscell(c.args) && mod(numel(c.args),2) == 0;
        if ~valid
            error('IMR:referenceCases','Invalid regression case %d.',k);
        end
    end
    modes = unique({cases.mode});
    for m = 1:numel(modes)
        keys = {cases(strcmp({cases.mode},modes{m})).key};
        if numel(unique(keys)) ~= numel(keys)
            error('IMR:referenceCases','Duplicate keys in mode %s.',modes{m});
        end
    end
end

function f_check_overwrite(files,overwrite)
    for k = 1:numel(files)
        if exist(files{k},'file') && ~overwrite
            error('IMR:referenceOverwrite', ...
                'Stored results exist: %s. Set overwrite = true to replace them.',files{k});
        end
    end
end

function [R,msg] = f_complete_case(c)
    [R,msg,t] = f_test_run(c);
    if ~isempty(msg)
        return
    end
    names = c.args(1:2:end);
    index = find(strcmpi(names,'tvector'),1,'last');
    if isempty(index)
        index = find(strcmpi(names,'tfin'),1,'last');
        if isempty(index)
            msg = 'Regression cases must explicitly specify tvector or tfin.';
            return
        end
        expected = [0 c.args{2*index}];
    else
        expected = c.args{2*index};
    end
    tolerance = 32*eps(max(1,max(abs(expected(:)))));
    complete = ~isempty(t) && isequal(size(t),size(R)) && numel(t) >= 2 && ...
        all(diff(t) > 0) && abs(t(1)-expected(1)) <= tolerance && ...
        abs(t(end)-expected(end)) <= tolerance;
    if numel(expected) > 2
        complete = complete && numel(t) == numel(expected) && ...
            all(abs(t(:)-expected(:)) <= tolerance);
    end
    if ~complete
        msg = 'Incomplete or mismatched output times. No baseline accepted.';
    end
end

function f_remove_temporary(folder)
    if isfolder(folder)
        rmdir(folder,'s');
    end
end
