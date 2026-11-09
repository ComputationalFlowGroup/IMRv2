% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Execute a checked plan (see build_verification_plan); return failures for
% the controller's final report. Each stage names a suite script path or
% carries the generation callback directly, and marks whether it is
% dependent on every earlier stage having passed.
function report = execute_verification(plan)
    empty = struct('name','','status','skipped','seconds',0,'message','','identifier','');
    report = repmat(empty,1,numel(plan));
    for k = 1:numel(plan)
        stage = plan(k);
        blocked = stage.dependent && any(~strcmp({report(1:k-1).status},'passed'));
        if blocked
            report(k).name = stage.name;
            report(k).message = 'Blocked by an earlier failed stage. No dependent action was run.';
        else
            report(k) = f_run_stage(stage.name,f_stage_action(stage.run));
        end
    end
end

function action = f_stage_action(run)
    if ischar(run) || isstring(run)
        action = @() f_run_script(char(run));
    else
        action = run;
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
