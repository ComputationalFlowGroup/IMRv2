% file f_test_cases.m
% brief contains function f_test_cases

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function lists the test cases of a suite, for the three solver
% modes: imr.solve_nodal with finite differences (nodal_fd), imr.solve_nodal with
% Chebyshev collocation (nodal_cheb) and imr.solve_spectral (spectral). The
% cases are built from the model features and kept only when the solver's
% own input checks (imr.configuration.legacy_options and imr.configuration.validate_case) accept them, so every
% suite covers exactly the combinations of the compatibility matrix and
% follows it when it changes. suite is 'equilibrium', 'gaussian',
% 'regression' or 'permutations'. Each case has the fields
%   key     a name that identifies the case within its mode, e.g.
%           'r2_bubble+medium_s3_mat1'; stored results are looked up by it
%   mode    'nodal_fd', 'nodal_cheb' or 'spectral'
%   solver  the solver function handle
%   args    the dimensionless name/value pairs of the run, grids included
% Run imr_test_setup from tests/ before calling this helper.
function cases = f_test_cases(suite)
    
    % solver modes and their grids: the finite difference grid is fine
    % enough for second-order accuracy, and the Chebyshev and spectral
    % expansions converge spectrally on far fewer points. The spectral
    % solver's Nv stress modes (the stress 5-7 field, and the Zener field
    % of stress 3 and 4 with medium heat transfer) are set to 32, well
    % below the default 150
    modes = struct( ...
        'name',   {'nodal_fd', 'nodal_cheb', 'spectral'}, ...
        'solver', {@imr.solve_nodal, @imr.solve_nodal, @imr.solve_spectral}, ...
        'call',   {'nodal', 'nodal', 'spectral'}, ...
        'extra',  {{'scheme', 0}, {'scheme', 1}, {}});
    
    % thermal configurations (one per row: name, bubtherm, medtherm,
    % masstrans), rows separated by line breaks
    thermal = {
    'polytropic', 0, 0, 0
    'bubble', 1, 0, 0
    'bubble+medium', 1, 1, 0
    'bubble+mass', 1, 0, 1
    'bubble+medium+mass', 1, 1, 1
    };
    
    switch suite
        case 'equilibrium'
        % a bubble at equilibrium must stay there: every radial, thermal
        % (without mass transfer), vapor and stress model
        grids = {{'Nt', 150, 'Mt', 150}, {'Nt', 20, 'Mt', 20}, ...
            {'Nt', 12, 'Mt', 12, 'Nv', 32}};
        base = {'Req', 100e-6, 'R0', 100e-6, 't8', 298.15, ...
            'tvector', linspace(0, 10e-6, 100), 'collapse', 0};
        specs = {};
        for radial = 1:6
            for it = 1:3
                for vapor = 0:1
                    for stress = 0:8
                        key = sprintf('r%d_%s_v%d_s%d', radial, ...
                            thermal{it,1}, vapor, stress);
                        specs(end+1,:) = {key, [base, ...
                            f_thermal(thermal(it,:)), {'radial', radial, ...
                            'vapor', vapor, 'stress', stress}, ...
                            f_gies(stress)]}; %#ok<AGROW>
                    end
                end
            end
        end

        case 'gaussian'
        % a bubble returns to equilibrium after a Gaussian pressure pulse:
        % every radial, thermal (without mass transfer), vapor and viscous
        % stress model
        grids = {{'Nt', 50, 'Mt', 50}, {'Nt', 20, 'Mt', 20}, ...
            {'Nt', 12, 'Mt', 12, 'Nv', 32}};
        base = {'Req', 10e-6, 'R0', 10e-6, 't8', 298.15, 'mu', 0.5, ...
            'tvector', linspace(0, 2e-3, 1000), 'collapse', 0, ...
            'wave_type', 1, 'dt', 3e-4, 'tw', 1e-4, 'pA', 1e5};
        specs = {};
        for radial = 1:6
            for it = 1:3
                for vapor = 0:1
                    for stress = 1:8
                        key = sprintf('r%d_%s_v%d_s%d', radial, ...
                            thermal{it,1}, vapor, stress);
                        specs(end+1,:) = {key, [base, ...
                            f_thermal(thermal(it,:)), {'radial', radial, ...
                            'vapor', vapor, 'stress', stress}, ...
                            f_gies(stress)]}; %#ok<AGROW>
                    end
                end
            end
        end

        case 'regression'
        % stored results: every radial, thermal (with mass transfer) and
        % stress model for two materials, then the other features of the
        % compatibility matrix on a Keller-Miksis case
        grids = {{'Nt', 70, 'Mt', 70}, {'Nt', 20, 'Mt', 20}, ...
            {'Nt', 12, 'Mt', 12, 'Nv', 32}};
        % a soft, weakly viscous material and a stiff, viscous one with
        % strain stiffening and relaxation
        mats = {{'mu', 1e-4, 'G', 1e2, 'alphax', 1e-3, 'lambda1', 1e-7}, ...
            {'mu', 1e-1, 'G', 2.5e3, 'alphax', 1, 'lambda1', 1e-3}};
        % the collapse and rebounds of a bubble released from four times
        % its equilibrium radius, and the same with mass transfer for a
        % larger bubble. The polytropic bubble loses no heat, so its
        % collapse is far sharper, and it is released from twice its
        % equilibrium radius. The runs are strong enough to exercise every
        % term but not so violent that the result mostly tests the
        % integrator
        osc = {'R0', 50e-6, 'Req', 12.5e-6, 'tvector', linspace(0, 12e-6, 100)};
        poly = {'R0', 50e-6, 'Req', 25e-6, 'tvector', linspace(0, 12e-6, 100)};
        mass = {'R0', 2e-4, 'Req', 5e-5, 'tvector', linspace(0, 30e-6, 100)};
        % collapse: the first collapse and rebound only
        coll = {'R0', 2e-4, 'Req', 5e-5, 'tvector', linspace(0, 20e-6, 100)};
        common = {'lambda2', 0, 'vapor', 1, 'collapse', 0};
        specs = {};
        for radial = 1:6
            for it = 1:5
                for stress = [0:5 8]
                    for im = 1:2
                        if thermal{it,4}
                            geom = mass;
                        elseif it == 1
                            geom = poly;
                        else
                            geom = osc;
                        end
                        key = sprintf('r%d_%s_s%d_mat%d', radial, ...
                            thermal{it,1}, stress, im);
                        specs(end+1,:) = {key, [common, geom, mats{im}, ...
                            f_thermal(thermal(it,:)), {'radial', radial, ...
                            'stress', stress}]}; %#ok<AGROW>
                    end
                end
            end
        end
        % Keller-Miksis with bubble and medium heat transfer
        km = [common, osc, mats{1}, f_thermal(thermal(3,:)), {'radial', 2}];
        % non-Newtonian viscosity, every model, with Kelvin-Voigt stress
        for nu = 1:7
            specs(end+1,:) = {sprintf('r2_nu%d', nu), [km, {'stress', 1, ...
                'nu_model', nu, 'mu0', 5e-2, 'v_lambda', 1e-6, 'v_nc', 0.5, ...
                'v_a', 2}]}; %#ok<AGROW>
        end
        % driving waveforms: Gaussian, histotripsy and Heaviside
        specs(end+1,:) = {'r2_wave1', [km, {'stress', 1, 'wave_type', 1, ...
            'pA', 1e5, 'tw', 1e-6, 'dt', 3e-6}]};
        specs(end+1,:) = {'r2_wave2', [km, {'stress', 1, 'wave_type', 2, ...
            'pA', 1e6, 'omega', 2*pi*1e6, 'dt', 2e-6, 'mn', 3.7}]};
        specs(end+1,:) = {'r2_wave3', [km, {'stress', 1, 'wave_type', 3, ...
            'pA', 1e5, 'tw', 5e-6}]};
        % UCM / Oldroyd-B with the stress field (spectral solver only)
        specs(end+1,:) = {'r2_s5_field', [km, {'stress', 5, 'spectral', 1}]};
        % collapse, which needs the full model with mass transfer, with the
        % soft material: with the stiff one the quadratic model's collapse
        % is too violent to converge with resolution
        for stress = [0 1 3 4]
            specs(end+1,:) = {sprintf('r2_collapse_s%d', stress), [common, ...
                coll, mats{1}, f_thermal(thermal(5,:)), {'radial', 2, ...
                'stress', stress}]}; %#ok<AGROW>
            specs{end,2}{find(strcmp(specs{end,2}, 'collapse'), 1) + 1} = 1;
        end
        
        case 'permutations'
        % every accepted combination of the solver options that the other
        % suites do not cross, run briefly on coarse grids: the suite only
        % checks that each run finishes with finite results
        grids = {{'Nt', 20, 'Mt', 20}, {'Nt', 12, 'Mt', 12}, ...
            {'Nt', 8, 'Mt', 8, 'Nv', 8}};
        base = {'R0', 50e-6, 'Req', 25e-6, 'mu', 1e-2, 'G', 1e3, ...
            'lambda1', 1e-7, 'lambda2', 3e-8, 'alphax', 0.25, ...
            'tvector', linspace(0, 2e-6, 5), 'collapse', 0, 'vapor', 1};
        nonnewt = {'mu0', 5e-2, 'v_lambda', 1e-7, 'v_nc', 0.5, 'v_a', 2};
        wave = {'pA', 1e5, 'TW', 2e-7, 'DT', 1e-6, 'omega', 2*pi*1e6, ...
            'mn', 2};
        specs = {};
        % non-Newtonian viscosity: every model with every viscous stress
        % model, with and without medium heat transfer
        for nu = 1:7
            for stress = 1:8
                for it = 2:3
                    specs(end+1,:) = {sprintf('nu%d_s%d_%s', nu, stress, ...
                        thermal{it,1}), [base, nonnewt, ...
                        f_thermal(thermal(it,:)), {'radial', 2, ...
                        'stress', stress, 'nu_model', nu}, ...
                        f_gies(stress)]}; %#ok<AGROW>
                end
            end
        end
        % driving waveforms with every radial model
        for w = -3:3
            for radial = 1:6
                specs(end+1,:) = {sprintf('wave%d_r%d', w, radial), ...
                    [base, wave, f_thermal(thermal(1,:)), {'radial', ...
                    radial, 'stress', 1, 'wave_type', w}]}; %#ok<AGROW>
            end
        end
        % integrators and Jacobians with every thermal configuration and
        % the stress model families
        for method = [3 15 23]
            for jacobian = 0:3
                for it = 1:5
                    for stress = [0 1 3 5 6 7 8]
                        specs(end+1,:) = {sprintf('m%d_j%d_%s_s%d', method, ...
                            jacobian, thermal{it,1}, stress), [base, ...
                            f_thermal(thermal(it,:)), {'radial', 2, ...
                            'stress', stress, 'method', method, ...
                            'jacobian', jacobian}, f_gies(stress)]}; %#ok<AGROW>
                    end
                end
            end
        end
        % no viscosity, surface tension or elasticity (Re, We, Ca = Inf)
        % (one pair per row, rows separated by line breaks). PTT (6) and
        % Giesekus (7) multiply their active term by Re8, which is 0*Inf
        % (NaN) from rest at the inviscid limit (mu = 0, Re8 = Inf), so that
        % limit is skipped for them: it is not a case those models are
        % meant to run at. Stress 8's viscous term only divides by Re8, so
        % it is unaffected and stays in the sweep.
        limits = {
        'mu', 0
        'surft', 0
        'G', 0
        };
        for il = 1:3
            for stress = 0:8
                if (stress == 6 || stress == 7) && strcmp(limits{il,1}, 'mu')
                    continue
                end
                specs(end+1,:) = {sprintf('inf_%s_s%d', limits{il,1}, ...
                    stress), [base, f_thermal(thermal(3,:)), {'radial', 2, ...
                    'stress', stress, limits{il,1}, limits{il,2}}, ...
                    f_gies(stress)]}; %#ok<AGROW>
            end
        end
        % collapse with every stress model (stress 5 needs stress0) and the
        % limits above
        for radial = 1:2
            for stress = [0 1 2 3 4 5 8]
                for il = 0:3
                    a = [base, f_thermal(thermal(5,:)), {'radial', radial, ...
                        'stress', stress, 'R0', 2e-4, 'Req', 5e-5, ...
                        'tvector', linspace(0, 1e-6, 5)}];
                    a{find(strcmp(a, 'collapse'), 1) + 1} = 1;
                    if il > 0
                        a = [a, limits(il,:)]; %#ok<AGROW>
                    end
                    if stress == 5
                        a = [a, {'stress0', [0 0]}]; %#ok<AGROW>
                    end
                    specs(end+1,:) = {sprintf('collapse_r%d_s%d_lim%d', ...
                        radial, stress, il), a}; %#ok<AGROW>
                end
            end
        end
        % the stress field of stress 5, the bubble-wall nonlinearity eps3,
        % a moving start, a given initial pressure and stress, the end time
        % given as tfin, and the integrator named instead of numbered
        km = [base, f_thermal(thermal(3,:)), {'radial', 2}];
        specs(end+1,:) = {'s5_field', [km, {'stress', 5, 'spectral', 1}]};
        specs(end+1,:) = {'eps3', [km, {'stress', 1, 'eps3', 0.5}]};
        specs(end+1,:) = {'rdotzero', [km, {'stress', 3, 'U0', -5}]};
        specs(end+1,:) = {'p0', [km, {'stress', 1, 'P0', 5e4}]};
        specs(end+1,:) = {'stress0_s3', [km, {'stress', 3, 'stress0', 1e-2}]};
        specs(end+1,:) = {'stress0_s5', [base, f_thermal(thermal(2,:)), ...
            {'radial', 2, 'stress', 5, 'stress0', [1e-2 2e-2]}]};
        ktv = find(strcmp(km, 'tvector'), 1);
        specs(end+1,:) = {'tfin', [km([1:ktv-1, ktv+2:end]), {'stress', 1, ...
            'tfin', 2e-6}]};
        specs(end+1,:) = {'method_name', [km, {'stress', 1, 'method', ...
            'ode23tb'}]};
        
        otherwise
        error(['f_test_cases: suite must be equilibrium, gaussian, ' ...
            'regression or permutations']);
    end
    
    % every case in every mode, kept when the solver's input checks accept it
    cases = struct('key', {}, 'mode', {}, 'solver', {}, 'args', {}, 't_request', {});
    for m = 1:numel(modes)
        for k = 1:size(specs, 1)
            args = [imr.nondimensionalize(specs{k,2}{:}), grids{m}, modes(m).extra];
            if f_accepted(modes(m).call, args)
                cases(end+1) = struct('key', specs{k,1}, ...
                    'mode', modes(m).name, 'solver', modes(m).solver, ...
                    'args', {args}, 't_request', f_t_request(args)); %#ok<AGROW>
            end
        end
    end

end

% The requested final time of a case: the last entry of tvector if given,
% else tfin, with the last occurrence of either winning (matching how
% imr.nondimensionalize itself resolves a name repeated in args, by
% sequential overwrite). Computed once here, where the full concatenated
% args list for this case/mode/grid is already in hand, instead of each
% test suite re-deriving it downstream with its own ad hoc, easily
% mismatched strcmp/strcmpi lookup.
function t_request = f_t_request(args)
    kt = find(strcmp(args, 'tvector'), 1, 'last');
    if ~isempty(kt)
        t_request = args{kt+1}(end);
        return
    end
    kf = find(strcmpi(args, 'tfin'), 1, 'last');
    if ~isempty(kf)
        t_request = args{kf+1};
        return
    end
    error('f_test_cases: case has neither tvector nor tfin.');
end

% name/value pairs of a row of the thermal table
function pairs = f_thermal(row)
    pairs = {'bubtherm', row{2}, 'medtherm', row{3}, 'masstrans', row{4}};
end

% a nonzero Giesekus mobility factor for stress 7, so its cases exercise
% the quadratic term and are not the eps3 = 0 limit that reduces exactly
% to stress 5; empty (no override) for every other stress
function pairs = f_gies(stress)
    if stress == 7
        pairs = {'eps3', 0.3};
    else
        pairs = {};
    end
end

% true when imr.configuration.legacy_options (and so imr.configuration.validate_case) accepts the inputs, with its
% console output hidden; only input errors reject a case. call and args are
% used inside the evalc string, which the code analyzer does not read
function ok = f_accepted(call, args) %#ok<INUSD>
    ok = true;
    try
        evalc('imr.configuration.legacy_options(call, args{:})');
    catch err
        if strncmp(err.message, 'INPUT ERROR', 11)
            ok = false;
        else
            rethrow(err);
        end
    end
end
