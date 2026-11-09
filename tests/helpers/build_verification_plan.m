% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Build the ordered list of stages a run_verification mode maps to, each a
% name plus either a suite script path or the generator's function handle,
% and whether it is dependent (grouped by execute_verification if an
% earlier stage fails). Shared by run_verification and its isolated tests.
function stages = build_verification_plan(mode,files,generator,regression_file)
    stages = struct('name',{},'run',{},'dependent',{});
    if strcmp(mode,'generate')
        stages(1) = struct('name','generate_references','run',generator,'dependent',false);
        return
    end
    for k = 1:numel(files)
        stages(k) = struct('name',f_stage_name(files{k}),'run',files{k},'dependent',false);
    end
    if strcmp(mode,'generate-and-test')
        stages(end+1) = struct('name','generate_references','run',generator,'dependent',true);
        stages(end+1) = struct('name',f_stage_name(regression_file),'run',regression_file,'dependent',true);
    end
end

function name = f_stage_name(file)
    [~,name] = fileparts(file);
end
