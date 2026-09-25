% file f_test_cases.m
% brief contains function f_test_cases

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function lists the test cases of a suite, for the three solver
% modes: f_imr_nodal with finite differences (nodal_fd), f_imr_nodal with
% Chebyshev collocation (nodal_cheb) and f_imr_spectral (spectral). The
% cases are built from the model features and kept only when the solver's
% own input checks (f_call_params and f_checkers) accept them, so every
% suite covers exactly the combinations of the compatibility matrix and
% follows it when it changes. suite is 'equilibrium', 'gaussian' or
% 'regression'. Each case has the fields
%   key     a name that identifies the case within its mode, e.g.
%           'r2_bubble+medium_s3_mat1'; stored results are looked up by it
%   mode    'nodal_fd', 'nodal_cheb' or 'spectral'
%   solver  the solver function handle
%   args    the dimensionless name/value pairs of the run, grids included
% Run it from toolchain/, with src/ on the path.
function cases = f_test_cases(suite)
    
    % solver modes and their grids: the finite difference grid is fine
    % enough for second-order accuracy, and the Chebyshev and spectral
    % expansions converge spectrally on far fewer points
    modes = struct( ...
        'name',   {'nodal_fd', 'nodal_cheb', 'spectral'}, ...
        'solver', {@f_imr_nodal, @f_imr_nodal, @f_imr_spectral}, ...
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
            {'Nt', 12, 'Mt', 12}};
        base = {'Req', 100e-6, 'R0', 100e-6, 't8', 298.15, ...
            'tvector', linspace(0, 10e-6, 100), 'collapse', 0};
        specs = {};
        for radial = 1:6
            for it = 1:3
                for vapor = 0:1
                    for stress = 0:5
                        key = sprintf('r%d_%s_v%d_s%d', radial, ...
                            thermal{it,1}, vapor, stress);
                        specs(end+1,:) = {key, [base, ...
                            f_thermal(thermal(it,:)), {'radial', radial, ...
                            'vapor', vapor, 'stress', stress}]}; %#ok<AGROW>
                    end
                end
            end
        end
        
        case 'gaussian'
        % a bubble returns to equilibrium after a Gaussian pressure pulse:
        % every radial, thermal (without mass transfer), vapor and viscous
        % stress model
        grids = {{'Nt', 50, 'Mt', 50}, {'Nt', 20, 'Mt', 20}, ...
            {'Nt', 12, 'Mt', 12}};
        base = {'Req', 10e-6, 'R0', 10e-6, 't8', 298.15, 'mu', 0.5, ...
            'tvector', linspace(0, 2e-3, 1000), 'collapse', 0, ...
            'wave_type', 1, 'dt', 3e-4, 'tw', 1e-4, 'pA', 1e5};
        specs = {};
        for radial = 1:6
            for it = 1:3
                for vapor = 0:1
                    for stress = 1:5
                        key = sprintf('r%d_%s_v%d_s%d', radial, ...
                            thermal{it,1}, vapor, stress);
                        specs(end+1,:) = {key, [base, ...
                            f_thermal(thermal(it,:)), {'radial', radial, ...
                            'vapor', vapor, 'stress', stress}]}; %#ok<AGROW>
                    end
                end
            end
        end
        
        case 'regression'
        % stored results: every radial, thermal (with mass transfer) and
        % stress model for two materials, then the other features of the
        % compatibility matrix on a Keller-Miksis case
        grids = {{'Nt', 70, 'Mt', 70}, {'Nt', 20, 'Mt', 20}, ...
            {'Nt', 12, 'Mt', 12}};
        % a soft, weakly viscous material and a stiff, viscous one with
        % strain stiffening and relaxation
        mats = {{'mu', 1e-4, 'G', 1e2, 'alphax', 1e-3, 'lambda1', 1e-7}, ...
            {'mu', 1e-1, 'G', 2.5e3, 'alphax', 1, 'lambda1', 1e-3}};
        % the collapse and rebounds of a bubble released from four times
        % its equilibrium radius, and the same with mass transfer for a
        % larger bubble. They are strong enough to exercise every term but
        % not so violent that the result mostly tests the integrator: GNU
        % Octave's weaker ode15s completes every case
        osc = {'R0', 50e-6, 'Req', 12.5e-6, 'tvector', linspace(0, 12e-6, 100)};
        mass = {'R0', 2e-4, 'Req', 5e-5, 'tvector', linspace(0, 30e-6, 100)};
        common = {'lambda2', 0, 'vapor', 1, 'collapse', 0};
        specs = {};
        for radial = 1:6
            for it = 1:5
                for stress = 0:5
                    for im = 1:2
                        if thermal{it,4}
                            geom = mass;
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
        specs(end+1,:) = {'r2_s5_field', [km, {'stress', 5, 'spectral', 1, ...
            'Nv', 40}]};
        % collapse, which needs the full model with mass transfer
        for stress = [0 1 3 4]
            specs(end+1,:) = {sprintf('r2_collapse_s%d', stress), [common, ...
                mass, mats{2}, f_thermal(thermal(5,:)), {'radial', 2, ...
                'stress', stress}]}; %#ok<AGROW>
            specs{end,2}{find(strcmp(specs{end,2}, 'collapse'), 1) + 1} = 1;
        end
        
        otherwise
        error('f_test_cases: suite must be equilibrium, gaussian or regression');
    end
    
    % every case in every mode, kept when the solver's input checks accept it
    cases = struct('key', {}, 'mode', {}, 'solver', {}, 'args', {});
    for m = 1:numel(modes)
        for k = 1:size(specs, 1)
            args = [f_dim2nondim(specs{k,2}{:}), grids{m}, modes(m).extra];
            if f_accepted(modes(m).call, args)
                cases(end+1) = struct('key', specs{k,1}, ...
                    'mode', modes(m).name, 'solver', modes(m).solver, ...
                    'args', {args}); %#ok<AGROW>
            end
        end
    end
    
end

% name/value pairs of a row of the thermal table
function pairs = f_thermal(row)
    pairs = {'bubtherm', row{2}, 'medtherm', row{3}, 'masstrans', row{4}};
end

% true when f_call_params (and so f_checkers) accepts the inputs, with its
% console output hidden; only input errors reject a case
function ok = f_accepted(call, args)
    ok = true;
    try
        evalc('f_call_params(call, args{:})');
    catch err
        if strncmp(err.message, 'INPUT ERROR', 11)
            ok = false;
        else
            rethrow(err);
        end
    end
end
