% file s_generate_masstransfer.m
% brief contains script to generate golden data for IMR mass transfer cases

% Copyright (C) Mauro Rodriguez
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This script generates the golden data for cases using mass transfer
clc;
clear;
close;

addpath('../toolchain/');
addpath('../src/');
addpath('../tests/');
load('file_ids.mat');

shift = 4*2*2*6*2*2*2*2*2 - 4*1*1*6*2*2*2*2;

% mass transfer test case
tvector = linspace(0,50E-6,100);
masstrans = 1;
vapor = 1;
collapse = 0;
R0 = 2.0e-04;
Req = 3.5e-05;
radial_vec = 1:4;
bubtherm_vec = 1;
medtherm_vec = 1;
stress_vec = 0:5;

% set up combinations
dims = [length(stress_vec), length(medtherm_vec), ...
    length(bubtherm_vec), length(radial_vec)];
total_combinations = prod(dims);

filenames_mass = cell(total_combinations,1);
for idx = 1:total_combinations
    filenames_mass{idx} = sprintf('../tests/%s.mat', ids{idx+shift});
end

% ensure the pool is active
if isempty(gcp('nocreate'))
    parpool('local',12);
end

% parallel dispatch; the enthalpy forms (radial 3, 4) take the liquid at
% the far-field temperature, so their medium heat transfer cases are not
% run (their file indices are kept so the others map to the same files)
[~, med_v, ~, rad_v] = ind2sub(dims, 1:total_combinations);
run_idx = find(~(radial_vec(rad_v) >= 3 & medtherm_vec(med_v) == 1));
futures(numel(run_idx)) = parallel.FevalFuture;

for q = 1:numel(run_idx)
    idx_mass = run_idx(q);
    [stress_idx, medtherm_idx, bubtherm_idx, radial_idx] = ind2sub(dims, idx_mass);
    
    radial   = radial_vec(radial_idx);
    bubtherm = bubtherm_vec(bubtherm_idx);
    medtherm = medtherm_vec(medtherm_idx);
    stress   = stress_vec(stress_idx);
    
    varin = {'radial', radial, ...
        'bubtherm', bubtherm, ...
        'tvector', tvector, ...
        'vapor', vapor, ...
        'medtherm', medtherm, ...
        'stress', stress, ...
        'collapse', collapse, ...
        'r0', R0, ...
        'req', Req, ...
        'masstrans', masstrans};
    % dimensional (SI) values converted to the solver's dimensionless inputs
    varin = f_dim2nondim(varin{:});
    
    futures(q) = parfeval(@f_imr_nodal_wrapper, 1, varin);
end

% save results in the correct order
for q = 1:numel(run_idx)
    [completedQ, result] = fetchNext(futures);
    completedIdx = run_idx(completedQ);
    f_savefile_fd(filenames_mass{completedIdx}, result);
    fprintf('Saved result %d to %s\n', completedIdx, filenames_mass{completedIdx});
end

delete(gcp('nocreate'));

% deterministic setup per worker
function Rm = f_imr_nodal_wrapper(varin)
    % fixed RNG seed
    rng(12345, 'twister');
    [~, Rm] = f_imr_nodal(varin{:}, 'Nt', 70, 'Mt', 70);
end


% savefile function
function f_savefile_fd(filename, data)
    Rm = data;
    save(filename, 'Rm');
end
