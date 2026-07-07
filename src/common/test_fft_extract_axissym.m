clear
close all
clc 

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir)
addpath(fullfile(this_dir, '..', 'Common_functions'))

max_mode = 20;
theta = linspace(0, 2 .* pi, 1201).';
theta(end) = [];

radius_offset = 3.2;
specified_modes = [2; 5; 8; 13; 18];
specified_amplitudes = [0.25; -0.12; 0.04; 0.015; -0.01];

radius = radius_offset .* ones(size(theta));
for k = 1:numel(specified_modes)
    n = specified_modes(k);
    radius = radius + specified_amplitudes(k) .* ...
        spherical_harmonic(cos(theta), 0, n, 0);
end

shuffle_order = randperm(numel(theta));
simdata = [theta(shuffle_order), radius(shuffle_order)];

[mode_extractf, amp_extractf, phase_extractf] = fft_extract_axissym(simdata, max_mode);

expected_modes = (0:max_mode).';
expected_amplitudes = zeros(max_mode + 1, 1);
expected_amplitudes(1) = radius_offset;
expected_amplitudes(specified_modes + 1) = specified_amplitudes;

reconstructed_radius = amp_extractf(1) .* ones(size(theta));
for n = 1:max_mode
    reconstructed_radius = reconstructed_radius + amp_extractf(n + 1) .* ...
        spherical_harmonic(cos(theta), 0, n, 0);
end

epn_axisym = normalize_axisym_amplitudes(amp_extractf, amp_extractf(1));
expected_epn = expected_amplitudes ./ radius_offset;
expected_epn(1) = radius_offset;

amplitude_error = amp_extractf - expected_amplitudes;
max_amplitude_error = max(abs(amplitude_error));
max_normalized_error = max(abs(epn_axisym - expected_epn));
max_reconstruction_error = max(abs(reconstructed_radius - radius));

assert(isequal(mode_extractf, expected_modes), 'Extracted modes are not 0:max_mode.')
assert(max_amplitude_error < 1e-10, ...
    'Axisymmetric amplitude extraction failed. Max error: %.3g', max_amplitude_error)
assert(max_normalized_error < 1e-10, ...
    'Axisymmetric normalization failed. Max error: %.3g', max_normalized_error)
assert(max_reconstruction_error < 1e-10, ...
    'Axisymmetric reconstruction failed. Max error: %.3g', max_reconstruction_error)

selected_rows = [1; specified_modes + 1];
disp(table(mode_extractf(selected_rows), expected_amplitudes(selected_rows), ...
    amp_extractf(selected_rows), phase_extractf(selected_rows), ...
    'VariableNames', {'Mode', 'ExpectedAmplitude', 'ExtractedAmplitude', 'Phase'}))
fprintf('Axisymmetric extraction test passed. Max amplitude error = %.3g\n', ...
    max_amplitude_error)
