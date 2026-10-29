% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Execute a checked plan; return failures for the controller's final report.
% Paths and the generation callback are explicit for isolated fixture tests.
function report = execute_validation(mode,suite_files,generator,regression_file)
    generation_index = 0;
    if strcmp(mode,'generate')
        suite_files = {''};
        generation_index = 1;
    elseif strcmp(mode,'generate-and-test')
        generation_index = numel(suite_files)+1;
        suite_files = [suite_files,{''},{regression_file}];
    end
    empty = struct('name','','status','skipped','seconds',0,'message','','identifier','');
    report = repmat(empty,1,numel(suite_files));
    for k = 1:numel(suite_files)
        if k == generation_index
            name = 'generate_references';
            action = generator;
        else
            [~,name] = fileparts(suite_files{k});
            action = @() f_run_script(suite_files{k});
        end
        dependent = strcmp(mode,'generate-and-test') && k >= generation_index;
        blocked = dependent && any(~strcmp({report(1:k-1).status},'passed'));
        if blocked
            report(k).name = name;
            report(k).message = 'Blocked by an earlier failed stage. No dependent action was run.';
        else
            report(k) = f_run_stage(name,action);
        end
    end
end

function result = f_run_stage(name,action)
    % Keep restoration objects outside the workspace the script can clear.
    previous_path = path;
    previous_folder = pwd;
    previous_warnings = warning;
    cleanup = onCleanup(@() f_restore(previous_path,previous_folder,previous_warnings));
    result = struct('name',name,'status','passed','seconds',0,'message','','identifier','');
    fprintf('\nRunning %s...\n',name);
    stage_timer = tic;
    try
        action();
    catch err
        result.status = 'failed';
        result.message = err.message;
        result.identifier = err.identifier;
        fprintf('%s\n',getReport(err,'extended','hyperlinks','off'));
    end
    result.seconds = toc(stage_timer);
end

function f_run_script(script_file)
    % No runner variables or cleanup objects share this script workspace.
    run(script_file);
end

function f_restore(previous_path,previous_folder,previous_warnings)
    path(previous_path);
    cd(previous_folder);
    warning(previous_warnings);
end
