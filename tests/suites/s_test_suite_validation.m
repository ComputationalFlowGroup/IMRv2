% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Controller and generator checks use temporary fixtures, never real baselines.
clear;
fixture_root = tempname;
mkdir(fixture_root);
cleanup = onCleanup(@() rmdir(fixture_root,'s'));
passing = fullfile(fixture_root,'passing.m');
failing = fullfile(fixture_root,'failing.m');
after = fullfile(fixture_root,'after.m');
regression = fullfile(fixture_root,'regression.m');
marker = fullfile(fixture_root,'generated.marker');
extra = fullfile(fixture_root,'extra');
mkdir(extra);
f_write(passing,sprintf('clear\nassert(~exist(''mode'',''var''))\nassert(~exist(''report'',''var''))\n'));
f_write(failing,sprintf(['clear\nwarning(''off'',''IMR:fixtureWarning'')\n' ...
    'addpath(fullfile(fileparts(mfilename(''fullpath'')),''extra''))\n' ...
    'cd(fullfile(fileparts(mfilename(''fullpath'')),''extra''))\n' ...
    'error(''IMR:fixtureFailure'',''deliberate fixture failure'')\n']));
f_write(after,sprintf(['clear\n' ...
    'fid = fopen(fullfile(fileparts(mfilename(''fullpath'')),''continued.marker''),''w'')\n' ...
    'fclose(fid)\n']));
f_write(regression,sprintf(['clear\n' ...
    'assert(isfile(fullfile(fileparts(mfilename(''fullpath'')),''generated.marker'')))\n']));
generation = @() f_write(marker,'generated');

% Expected fixture failures must not stop later independent scripts, erase
% the controller, change its path/folder/warnings or invoke generation.
previous_path = path;
previous_folder = pwd;
previous_warning = warning('query','IMR:fixtureWarning');
report = execute_validation('tests',{passing,failing,after},generation,regression);
assert(isequal({report.status},{'passed','failed','passed'}));
assert(strcmp(report(2).identifier,'IMR:fixtureFailure'));
assert(isfile(fullfile(fixture_root,'continued.marker')) && ~isfile(marker));
assert(strcmp(path,previous_path) && strcmp(pwd,previous_folder));
assert(isequal(warning('query','IMR:fixtureWarning'),previous_warning));
assert(all([report.seconds] >= 0));
report = execute_validation('tests',{regression,passing},generation,regression);
assert(isequal({report.status},{'failed','passed'}) && ~isfile(marker));

% All prechecks continue, but generation and regression are blocked if one
% fails. A generation failure also blocks regression, never marking a pass.
report = execute_validation('generate-and-test',{failing,passing},generation,regression);
assert(isequal({report.status},{'failed','passed','skipped','skipped'}));
assert(~isfile(marker));
report = execute_validation('generate-and-test',{passing},@f_fail_generation,regression);
assert(isequal({report.status},{'passed','failed','skipped'}));
assert(~isfile(marker));
report = execute_validation('generate-and-test',{passing},generation,regression);
assert(isequal({report.status},{'passed','passed','passed'}));
assert(isfile(marker));
delete(marker);
report = execute_validation('generate',{failing},generation,regression);
assert(numel(report) == 1 && strcmp(report.status,'passed') && isfile(marker));
fprintf('Validation mode sequencing, workspace isolation and failure reports PASSED.\n');

% Public selection guards run before simulations or any generation.
f_error(@() run_validation('unknown'), 'IMR:validationOptions');
f_error(@() run_validation('tests',{}), 'IMR:validationOptions');
f_error(@() run_validation('tests',{'not_a_suite'}), 'IMR:validationOptions');
f_error(@() run_validation('tests','all',true), 'IMR:validationOptions');
f_error(@() run_validation('generate',{'solver'}), 'IMR:validationOptions');
f_error(@() run_validation('generate-and-test',{'regression'}), 'IMR:validationOptions');
% Exercise the actual controller on a small independent suite, avoiding
% recursion into this suite and dependence on local golden data availability.
report = run_validation('tests',{'diffmat','diffmat'});
assert(numel(report) == 1 && strcmp(report.name,'s_test_suite_diffmat'));
assert(strcmp(report.status,'passed') && strcmp(path,previous_path));
assert(strcmp(pwd,previous_folder));
fprintf('Public controller options, subset selection and deduplication PASSED.\n');

% Synthetic trajectories test publishing and overwrite guards in a temporary
% destination. Serial execution does not need a toolbox or start a pool.
refdir = fullfile(fixture_root,'reference');
a = struct('mode','nodal_fd','key','fixture','solver',@f_fixture_solver, ...
    'args',{{'tvector',[0 0.5 1]}});
b = a;
b.mode = 'spectral';
files = generate_references(false,refdir,[a b],false);
assert(numel(files) == 2 && all(cellfun(@isfile,files)));
data = load(files{1});
assert(isequal(data.keys,{'fixture'}) && isequal(data.R,{ones(3,1)}));
assert(ischar(data.info) && contains(data.info,version));
before = f_bytes(files{1});
f_error(@() generate_references(false,refdir,[a b],false),'IMR:referenceOverwrite');
assert(isequal(f_bytes(files{1}),before));

% No valid history may replace a baseline if another case fails, stops early
% or returns different requested output times, even with overwrite enabled.
bad = b;
bad.mode = 'failed';
bad.solver = @f_bad_pressure;
f_error(@() generate_references(true,refdir,[a bad],false),'IMR:referenceCases');
assert(~isfile(fullfile(refdir,'failed.mat')) && isequal(f_bytes(files{1}),before));
bad.solver = @f_incomplete_solver;
f_error(@() generate_references(true,refdir,[a bad],false),'IMR:referenceCases');
assert(~isfile(fullfile(refdir,'failed.mat')) && isequal(f_bytes(files{1}),before));
bad.solver = @f_wrong_times;
f_error(@() generate_references(true,refdir,[a bad],false),'IMR:referenceCases');
assert(~isfile(fullfile(refdir,'failed.mat')) && isequal(f_bytes(files{1}),before));
f_error(@() generate_references(false,refdir,a([]),false),'IMR:referenceCases');
f_error(@() generate_references(false,refdir,[a a],false),'IMR:referenceCases');

a.args = {'tvector',[0 0.5 1],'radius',2};
generate_references(true,refdir,a,false);
data = load(files{1});
assert(isequal(data.R,{2*ones(3,1)}));
% Two-time spans permit adaptive intermediate steps, with matching endpoints.
a.args = {'tvector',[0 1]};
a.solver = @f_adaptive_solver;
generate_references(true,refdir,a,false);
fprintf('Golden-file format, overwrite protection and completion guards PASSED.\n');
clear cleanup;

function f_write(file,contents)
    fid = fopen(file,'w');
    assert(fid >= 0);
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid,'%s',contents);
end

function bytes = f_bytes(file)
    fid = fopen(file,'r');
    assert(fid >= 0);
    cleanup = onCleanup(@() fclose(fid));
    bytes = fread(fid,Inf,'*uint8');
end

function f_fail_generation()
    error('IMR:fixtureGeneration','deliberate fixture generation failure');
end

function f_error(action,identifier)
    try
        action();
    catch err
        assert(strcmp(err.identifier,identifier), ...
            'Expected %s, got %s: %s',identifier,err.identifier,err.message);
        return
    end
    error('Expected %s.',identifier);
end

function [t,R,U,P,T,Tm,kv] = f_fixture_solver(varargin)
    t = varargin{2}(:);
    radius = 1;
    if nargin > 2
        radius = varargin{4};
    end
    R = radius*ones(size(t));
    U = zeros(size(t));
    P = ones(size(t));
    T = [];
    Tm = ones(size(t));
    kv = [];
end

function [t,R,U,P,T,Tm,kv] = f_bad_pressure(varargin)
    [t,R,U,P,T,Tm,kv] = f_fixture_solver(varargin{:});
    P(2) = NaN;
end

function [t,R,U,P,T,Tm,kv] = f_incomplete_solver(varargin)
    [t,R,U,P,T,Tm,kv] = f_fixture_solver(varargin{:});
    t = t(1:end-1);
    R = R(1:end-1);
    U = U(1:end-1);
    P = P(1:end-1);
    Tm = Tm(1:end-1);
end

function [t,R,U,P,T,Tm,kv] = f_wrong_times(varargin)
    [t,R,U,P,T,Tm,kv] = f_fixture_solver(varargin{:});
    t(2) = t(2)+0.1;
end

function [t,R,U,P,T,Tm,kv] = f_adaptive_solver(varargin)
    args = varargin;
    args{2} = [0 0.4 0.8 1];
    [t,R,U,P,T,Tm,kv] = f_fixture_solver(args{:});
end
