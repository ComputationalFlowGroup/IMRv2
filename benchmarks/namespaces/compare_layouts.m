% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Compare a pre-namespace checkout with this checkout on the same MATLAB.
% report = compare_layouts(baseline_root,candidate_root,repetitions,limit)
% Default: five alternating paired timeit measurements, 3% slowdown limit.
% The tolerance accommodates measurement noise, not an assertion that 3%
% regression is acceptable for every workload. Lower it to enforce a tighter
% budget. The measured table is printed before any timing assertion fails.
function report = compare_layouts(baseline_root,candidate_root,repetitions,limit)
    if nargin < 3
        repetitions = 5;
    end
    if nargin < 4
        limit = 0.03;
    end
    validateattributes(repetitions,{'numeric'},{'real','finite','scalar','integer','>=',3});
    validateattributes(limit,{'numeric'},{'real','finite','scalar','nonnegative'});
    baseline_root = char(baseline_root);
    candidate_root = char(candidate_root);
    assert(isfolder(baseline_root) && isfolder(candidate_root));
    assert(isfile(fullfile(baseline_root,'src','f_imr_nodal.m')));
    assert(~isfile(fullfile(baseline_root,'src','+imr','solve_nodal.m')), ...
        'The baseline must contain flat implementations, not namespace adapters.');
    assert(isfile(fullfile(candidate_root,'src','+imr','solve_nodal.m')));
    previous_path = path;
    previous_directory = pwd;
    cleanup = onCleanup(@() f_restore(previous_path,previous_directory)); %#ok<NASGU>
    % Absolute paths are required before moving away from either checkout.
    cd(baseline_root);
    baseline_root = pwd;
    cd(previous_directory);
    cd(candidate_root);
    candidate_root = pwd;
    cd(tempdir);
    cases = {'nodal_fd','nodal_cheb','spectral','nodal_mass'};
    results = zeros(numel(cases),3);
    fprintf('MATLAB %s, %s, %d paired measurements, limit %.1f%%\n', ...
        version,computer,repetitions,100*limit);
    for k = 1:numel(cases)
        mode = cases{k};
        args = {'Req',1,'rdotzero',0.02,'radial',2,'stress',1,'ca',Inf, ...
            'bubtherm',1,'medtherm',1,'vapor',0,'masstrans',0, ...
            'tvector',linspace(0,1,64),'jacobian',2};
        if strcmp(mode,'nodal_fd')
            args = [args,{'Nt',32,'Mt',32}];
        elseif strcmp(mode,'nodal_cheb')
            args = [args,{'Nt',12,'Mt',12,'scheme',1}];
        elseif strcmp(mode,'nodal_mass')
            args = [args,{'Nt',20,'Mt',20,'vapor',1,'masstrans',1}];
        else
            args = [args,{'Nt',12,'Mt',12}];
        end
        base_solver = f_select(baseline_root,mode,false);
        base_output = f_once(base_solver,args);
        candidate_solver = f_select(candidate_root,mode,true);
        candidate_output = f_once(candidate_solver,args);
        f_equal(base_output,candidate_output);
        timings = zeros(repetitions,2);
        for repeat = 1:repetitions
            order = [1 2];
            if mod(repeat,2) == 0
                order = [2 1];
            end
            for side = order
                if side == 1
                    solver = f_select(baseline_root,mode,false);
                else
                    solver = f_select(candidate_root,mode,true);
                end
                timings(repeat,side) = timeit(@() f_once(solver,args));
            end
        end
        center = median(timings,1);
        results(k,:) = [center center(2)/center(1)];
        fprintf('%s: baseline %.6g s, namespace %.6g s, ratio %.4f\n', ...
            mode,results(k,1),results(k,2),results(k,3));
    end
    report = array2table(results,'VariableNames',{'baseline_s','namespace_s','ratio'}, ...
        'RowNames',cases);
    disp(report);
    assert(all(results(:,3) <= 1+limit), ...
        'IMR:namespaceTiming','A workload exceeded the configured slowdown limit. Investigate before merging.');
end

function solver = f_select(root,mode,namespaced)
    % Isolate source roots so a compatibility adapter cannot accidentally
    % measure the candidate during the baseline run. Keep MATLAB's own path.
    restoredefaultpath;
    addpath(fullfile(root,'src'));
    if strcmp(mode,'spectral')
        name = 'f_imr_spectral';
        if namespaced
            name = 'imr.solve_spectral';
        end
    else
        name = 'f_imr_nodal';
        if namespaced
            name = 'imr.solve_nodal';
        end
    end
    resolved = which(name);
    assert(startsWith(resolved,fullfile(root,'src')), ...
        'The selected solver resolved outside the requested checkout.');
    solver = str2func(name);
end

function out = f_once(solver,args)
    % Both versions capture the same console output and request eight outputs.
    out = cell(1,8);
    evalc(['[out{:}] = solver(args{:})' char(59)]);
end

function f_equal(a,b)
    for k = 1:7
        assert(isequal(size(a{k}),size(b{k})));
        if ~isempty(a{k})
            scale = max(1,max(abs(a{k}),[],'all'));
            assert(all(isfinite(a{k}),'all') && all(isfinite(b{k}),'all'));
            assert(max(abs(a{k}-b{k}),[],'all') <= 5e-9*scale, ...
                'The layouts produced different numerical outputs.');
        end
    end
    assert(isequaln(a{8},b{8}));
end

function f_restore(previous_path,previous_directory)
    path(previous_path);
    cd(previous_directory);
end
