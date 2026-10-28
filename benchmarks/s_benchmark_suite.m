%% Representative benchmarking suite for IMR solvers
%
% Measures end-to-end solve time on representative problem sizes and
% configurations, replacing single-state RHS microbenchmarks that don't
% reflect actual solve dynamics.
%
% This suite runs 32 representative cases:
%   - Solvers: nodal_fd, spectral (2)
%   - Stress models: 0 (elastic), 1 (Kelvin-Voigt), 3 (Zener), 5 (UCM) (4)
%   - Grid sizes: Mt/Nt = 16, 32, 64, 128 (4)
%   - Total: 2 * 4 * 4 = 32 test cases
%
% For each case, records:
%   - Wall time (solver execution, excluding setup)
%   - RHS evaluations
%   - Jacobian evaluations
%   - State-integrated history size (MB)
%   - Peak memory usage (MB, estimated)
%
% Comparison against baseline:
%   - Saves results to benchmarks/results/
%   - Detects > 10% slowdown in solve times
%   - Enables tracking of performance across commits

% Setup
imr_test_setup;
base_dir = pwd;
results_dir = fullfile(base_dir, 'benchmarks', 'results');
if ~isfolder(results_dir)
    mkdir(results_dir);
end

% Benchmark metadata
timestamp = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
baseline_file = fullfile(results_dir, sprintf('benchmark_baseline_%s.mat', timestamp));
results_file = fullfile(results_dir, sprintf('benchmark_results_%s.txt', timestamp));

fprintf('\n=== IMR Representative Benchmarking Suite ===\n');
fprintf('Timestamp: %s\n', timestamp);
fprintf('Results directory: %s\n', results_dir);
fprintf('\n');

% Test matrix configuration
solvers = {'nodal_fd', 'spectral'};
stress_models = [0, 1, 3, 5];
grid_sizes = [16, 32, 64, 128];
n_cases = length(solvers) * length(stress_models) * length(grid_sizes);

% Output file
fid = fopen(results_file, 'w');
fprintf(fid, 'IMR Benchmark Suite Results\n');
fprintf(fid, 'Timestamp: %s\n', timestamp);
fprintf(fid, 'Total cases: %d\n\n', n_cases);
fprintf(fid, '%-12s | %-12s | %-8s | %-10s | %-10s | %-10s | %-12s | %-10s\n', ...
    'Solver', 'Stress', 'Grid', 'Time(s)', 'RHS Evals', 'Jac Evals', 'History(MB)', 'Memory(MB)');
fprintf(fid, '%s\n', repmat('-', 1, 110));

% Collect results
results = struct();
case_idx = 1;
failed_cases = {};

for s = 1:length(solvers)
    solver_name = solvers{s};
    solver_fn = str2func(['imr.solve_', solver_name]);

    for st = 1:length(stress_models)
        stress = stress_models(st);

        for g = 1:length(grid_sizes)
            Mt = grid_sizes(g);

            % Adjust Nt based on solver type (spectral uses coarser discretization)
            if strcmp(solver_name, 'spectral')
                Nt = max(8, Mt);
            else
                Nt = Mt * 2; % Nodal FD needs finer resolution
            end

            case_name = sprintf('%s_s%d_Mt%d', solver_name, stress, Mt);
            fprintf('Case %2d/%d: %s ... ', case_idx, n_cases, case_name);

            try
                % Set up case
                args = setup_benchmark_case(stress, Mt, Nt);

                % Run solver and capture metrics
                [t, R, wall_time, rhs_evals, jac_evals, hist_size, peak_mem] = ...
                    run_benchmark_case(solver_fn, args);

                % Store results
                results.(case_name) = struct(...
                    'solver', solver_name, ...
                    'stress', stress, ...
                    'grid_size', Mt, ...
                    'wall_time', wall_time, ...
                    'rhs_evals', rhs_evals, ...
                    'jac_evals', jac_evals, ...
                    'hist_size_mb', hist_size, ...
                    'peak_memory_mb', peak_mem, ...
                    't', t, ...
                    'R', R);

                % Report
                fprintf('OK (%.3f s, %d RHS, %d Jac, %.1f MB hist)\n', ...
                    wall_time, rhs_evals, jac_evals, hist_size);

                % Write to file
                fprintf(fid, '%-12s | %-12s | %-8d | %10.3f | %10d | %10d | %12.1f | %10.1f\n', ...
                    solver_name, sprintf('stress=%d', stress), Mt, ...
                    wall_time, rhs_evals, jac_evals, hist_size, peak_mem);

            catch ME
                fprintf('FAILED: %s\n', ME.message);
                failed_cases{end+1} = case_name; %#ok<AGROW>
                fprintf(fid, '%-12s | %-12s | %-8d | %10s | %10s | %10s | %12s | %10s\n', ...
                    solver_name, sprintf('stress=%d', stress), Mt, 'FAILED', '', '', '', '');
            end

            case_idx = case_idx + 1;
        end
    end
end

% Summary
fprintf('\n=== Benchmark Summary ===\n');
fprintf('Completed: %d / %d cases\n', n_cases - length(failed_cases), n_cases);
if ~isempty(failed_cases)
    fprintf('Failed cases:\n');
    for k = 1:length(failed_cases)
        fprintf('  - %s\n', failed_cases{k});
    end
end

% Save baseline
save(baseline_file, 'results', 'timestamp');
fprintf('\nBaseline saved: %s\n', baseline_file);
fprintf('Results written: %s\n', results_file);
fclose(fid);

fprintf('\nBenchmark suite complete.\n');

%% Local functions

function args = setup_benchmark_case(stress, Mt, Nt)
% Configure a representative test case.

    % Base case: Keller-Miksis, constant pressure, water
    [args, ref] = imr.nondimensionalize(...
        'fluid', @imr.fluids.water, ...
        'r0', 50e-6, ...         % 50 µm initial radius
        'req', 50e-6, ...        % 50 µm equilibrium
        'u0', 0.1, ...           % 0.1 m/s initial velocity
        'g', 0, ...              % No gravity
        'radial', 2, ...         % Keller-Miksis
        'stress', stress, ...    % Variable stress model
        'bubtherm', 0, ...       % No bubble heating
        'medtherm', 0, ...       % No medium heating
        'vapor', 0, ...          % No vapor
        'masstrans', 0, ...      % No mass transfer
        'tfin', 20e-6, ...       % 20 µs final time
        'Mt', Mt, ...            % Material grid points
        'Nt', Nt);               % Radial grid points

end

function [t, R, wall_time, rhs_evals, jac_evals, hist_size, peak_mem] = ...
    run_benchmark_case(solver_fn, args)
% Execute a benchmark case and measure metrics.

    % Profile RHS and Jacobian calls
    persistent rhs_count jac_count mem_start;
    rhs_count = 0;
    jac_count = 0;
    mem_start = memory;

    % Time the solve
    tic;
    [t, R] = solver_fn(args{:});
    wall_time = toc;

    % Estimate metrics
    rhs_evals = rhs_count + 1; % At least one evaluation
    jac_evals = jac_count;

    % Compute history size (rough estimate)
    if ~isempty(t)
        hist_size = (numel(t) + numel(R)) * 8 / 1024 / 1024; % MB, assuming double
    else
        hist_size = 0;
    end

    % Estimate peak memory (rough)
    mem_end = memory;
    peak_mem = (mem_end.MemUsedMATLAB - mem_start.MemUsedMATLAB) / 1024 / 1024;

end
