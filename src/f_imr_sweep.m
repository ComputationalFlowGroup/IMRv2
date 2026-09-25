% file f_imr_sweep.m
% brief contains function f_imr_sweep

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function runs a solver over every combination of the given
% input values (a full grid) and collects the results. Call it as
%   [res, grid] = f_imr_sweep(solver, base, name1, values1, name2, values2)
% with any number of swept inputs, followed by options as name/value pairs.
% solver is a function handle such as @f_imr_nodal or @f_imr_spectral.
% base holds the inputs shared by every run: a cell array of name/value
% pairs ({} for the defaults) or a struct from f_imr_inputs. Each name_k is
% a solver input and values_k the values it takes, a numeric vector or a
% cell array. res is a struct array
% of size [numel(values1), numel(values2), ...] with the fields inputs (the
% swept name/value pairs of the run), out (the solver outputs, {t, R} by
% default) and error (the error message of a failed run, '' otherwise). A
% failed run does not stop the sweep. grid holds the swept names and
% values. Options:
%   'transform' function handle applied to the name/value pairs of each
%               run before the solver, e.g. @f_dim2nondim for SI inputs
%   'outputs'   number of solver outputs to keep (default 2: t, R)
%   'parallel'  true runs the grid with parfor on a parallel pool
%               (Parallel Computing Toolbox), false (default) in series
%   'quiet'     true (default) hides the solver's console output
function [res, grid] = f_imr_sweep(solver, base, varargin)
    
    if ~isa(solver, 'function_handle')
        error('f_imr_sweep: solver must be a function handle, e.g. @f_imr_nodal');
    end
    if isstruct(base)
        % a struct of inputs (f_imr_inputs) ahead of the swept pairs
        base = {base};
    elseif ~iscell(base)
        error(['f_imr_sweep: base must be a struct from f_imr_inputs or a ' ...
            'cell array of name/value pairs']);
    end
    
    % options and swept inputs
    transform = [];
    nout = 2;
    parallel = false;
    quiet = true;
    names = {};
    values = {};
    if mod(numel(varargin), 2) == 1
        error('f_imr_sweep: every swept input and option needs a value');
    end
    for n = 1:2:numel(varargin)
        name = varargin{n};
        v = varargin{n+1};
        switch lower(name)
            case 'transform'
            transform = v;
            case 'outputs'
            nout = v;
            case 'parallel'
            parallel = v;
            case 'quiet'
            quiet = v;
            otherwise
            if isnumeric(v) || islogical(v)
                v = num2cell(v(:)');
            end
            if ~iscell(v) || isempty(v)
                error(['f_imr_sweep: the values of %s must be a nonempty ' ...
                    'vector or cell array'], name);
            end
            names{end+1} = name; %#ok<AGROW>
            values{end+1} = v(:)'; %#ok<AGROW>
        end
    end
    if isempty(names)
        error('f_imr_sweep: give at least one input to sweep');
    end
    if ~isempty(transform) && ~isempty(base) && isstruct(base{1})
        error(['f_imr_sweep: a transform takes name/value pairs, so give ' ...
            'base as a cell array when using one']);
    end
    
    % every combination, in grid order (the first input varies fastest)
    dims = cellfun(@numel, values);
    nrun = prod(dims);
    runs = cell(1, nrun);
    idx = cell(1, numel(dims));
    for k = 1:nrun
        [idx{:}] = ind2sub([dims 1], k);
        pairs = cell(1, 2*numel(names));
        for q = 1:numel(names)
            pairs{2*q-1} = names{q};
            pairs{2*q} = values{q}{idx{q}};
        end
        runs{k} = pairs;
    end
    
    % the runs, on a parallel pool if requested (0 workers runs in series)
    out = cell(1, nrun);
    msg = cell(1, nrun);
    base = base(:)';
    if parallel
        nworkers = Inf;
    else
        nworkers = 0;
    end
    parfor (k = 1:nrun, nworkers)
        [out{k}, msg{k}] = f_sweep_run(solver, [base runs{k}], transform, ...
            nout, quiet);
    end
    
    res = reshape(struct('inputs', runs, 'out', out, 'error', msg), [dims 1]);
    grid = struct('names', {names}, 'values', {values});
    
end

% one run of the sweep, with its error message instead of an error
function [out, msg] = f_sweep_run(solver, args, transform, nout, quiet)
    
    out = cell(1, nout);
    msg = '';
    try
        if ~isempty(transform)
            args = transform(args{:});
        end
        if quiet
            % char(59) is the statement's semicolon, kept out of the string
            % literal because s_format splits lines with two semicolons
            evalc(['[out{1:nout}] = solver(args{:})' char(59)]);
        else
            [out{1:nout}] = solver(args{:});
        end
    catch err
        msg = err.message;
    end
    
end
