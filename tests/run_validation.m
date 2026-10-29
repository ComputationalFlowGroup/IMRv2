% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Run tests and/or explicitly generate golden data, using base MATLAB.
% report = run_validation(mode,suites,overwrite)
% Defaults: 'tests', 'all', false. Suite names may be full script names or
% suffixes, e.g. {'diffmat','solver'}. No generation occurs in tests mode.
% generate runs only the generator; it requires suites = 'all'.
% generate-and-test runs selected independent suites, then generates and
% runs regression only if those stages pass. At least one precheck is needed.
% Prints all results and raises an error after any failure or blocked stage.
function report = run_validation(mode,suites,overwrite)
    if nargin < 1
        mode = 'tests';
    end
    if nargin < 2
        suites = 'all';
    end
    if nargin < 3
        overwrite = false;
    end
    valid_text = (ischar(mode) && isrow(mode)) || (isstring(mode) && isscalar(mode));
    if ~valid_text || ~any(strcmpi(mode,{'tests','generate','generate-and-test'}))
        error('IMR:validationOptions','mode must be tests, generate or generate-and-test.');
    end
    mode = lower(char(mode));
    validateattributes(overwrite,{'logical'},{'scalar'},mfilename,'overwrite');
    if strcmp(mode,'tests') && overwrite
        error('IMR:validationOptions','overwrite is available only in a generation mode.');
    end

    test_root = fileparts(mfilename('fullpath'));
    suite_root = fullfile(test_root,'suites');
    available = dir(fullfile(suite_root,'s_test_suite_*.m'));
    [~,order] = sort({available.name});
    available = available(order);
    names = erase({available.name},'.m');
    selected = f_select(suites,names);
    if isempty(selected)
        error('IMR:validationOptions','No test suites were found.');
    end
    all_selected = numel(selected) == numel(names) && all(ismember(names,selected));
    if strcmp(mode,'generate') && ~all_selected
        error('IMR:validationOptions','generate creates all regression baselines. Use suites = ''all''.');
    end
    regression_file = fullfile(suite_root,'s_test_suite_regression.m');
    if strcmp(mode,'generate-and-test')
        selected(strcmp(selected,'s_test_suite_regression')) = [];
        if isempty(selected) || ~isfile(regression_file)
            error('IMR:validationOptions','generate-and-test needs independent prechecks and the regression suite.');
        end
    end
    files = cellfun(@(name) fullfile(suite_root,[name '.m']),selected,'UniformOutput',false);

    previous_path = path;
    previous_folder = pwd;
    cleanup = onCleanup(@() f_restore(previous_path,previous_folder));
    imr_test_setup;
    fprintf('IMR validation: %s, MATLAB %s\n',mode,version);
    report = execute_validation(mode,files,@() generate_references(overwrite),regression_file);
    fprintf('\nIMR validation summary (%s)\n',mode);
    for k = 1:numel(report)
        fprintf('%-32s %-7s %9.2f s\n',report(k).name, ...
            upper(report(k).status),report(k).seconds);
        if ~isempty(report(k).message)
            fprintf('  %s\n',report(k).message);
        end
    end
    failed = strcmp({report.status},'failed');
    skipped = strcmp({report.status},'skipped');
    fprintf('%d passed, %d failed, %d skipped; %.2f s in stages.\n', ...
        sum(~failed & ~skipped),sum(failed),sum(skipped),sum([report.seconds]));
    if any(failed | skipped)
        error('IMR:validationFailed','Validation failed. See the complete stage summary above.');
    end
end

function selected = f_select(suites,available)
    if ischar(suites) && isrow(suites)
        suites = {suites};
    elseif isstring(suites) && isvector(suites)
        suites = cellstr(suites);
    end
    if ~iscell(suites) || ~isvector(suites) || isempty(suites)
        error('IMR:validationOptions','suites must be ''all'' or a nonempty list of suite names.');
    end
    selected = cell(1,numel(suites));
    for k = 1:numel(suites)
        name = suites{k};
        valid_text = (ischar(name) && isrow(name)) || (isstring(name) && isscalar(name));
        if ~valid_text
            error('IMR:validationOptions','Each suite name must be a text scalar.');
        end
        name = lower(char(name));
        if strcmp(name,'all') && numel(suites) == 1
            selected = available;
            return
        end
        if ~startsWith(name,'s_test_suite_')
            name = ['s_test_suite_' name]; %#ok<AGROW>
        end
        if ~any(strcmp(name,available))
            error('IMR:validationOptions','Unknown suite: %s. Use ''all'' or a script in tests/suites.',name);
        end
        selected{k} = name;
    end
    selected = unique(selected,'stable');
end

function f_restore(previous_path,previous_folder)
    path(previous_path);
    cd(previous_folder);
end
