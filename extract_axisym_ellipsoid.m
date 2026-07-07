%% extract_axisym_ellipsoid.m
% Initialize a prolate ellipsoid and extract axisymmetric amplitudes.

clear
close all
clc

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'src', 'common'))
addpath src\common\

% Ellipsoid principal radii. The two minor principal radii are equal, and
% the major principal radius is aligned with the z axis.
minorRadius = 30;
majorRadius = 40;

maxMode = 20;
numTheta = 721;
theta = linspace(0, pi, numTheta).';

% Axisymmetric ellipsoid radius in spherical coordinates, with theta
% measured from the +z axis:
%   x^2 / minorRadius^2 + y^2 / minorRadius^2 + z^2 / majorRadius^2 = 1
radius = 1 ./ sqrt((sin(theta) .^ 2) ./ minorRadius .^ 2 + ...
    (cos(theta) .^ 2) ./ majorRadius .^ 2);

simdata = [theta, radius];

[mode_extractf, amp_extractf, phase_extractf] = ...
    fft_extract_axissym(simdata, maxMode);
epn_axisym = normalize_axisym_amplitudes(amp_extractf, amp_extractf(1));

radiusReconstructed = reconstructAxisymRadius(theta, amp_extractf);
maxRadiusError = max(abs(radiusReconstructed - radius));

resultsTable = table(mode_extractf, amp_extractf, phase_extractf, epn_axisym, ...
    'VariableNames', {'Mode', 'Amplitude', 'Phase', 'NormalizedAmplitude'});
disp(resultsTable)
fprintf('Max reconstructed radius error with modes 0:%d = %.6g\n', ...
    maxMode, maxRadiusError)

% Coordinates of the meridional cross-section, useful for quick visual checks.
x = radius .* sin(theta);
z = radius .* cos(theta);

figure('Color', 'w')
plot(x, z, 'k-', -x, z, 'k-', 'LineWidth', 1.5)
axis equal
box on
grid on
xlabel('x')
ylabel('z')
title('Prolate ellipsoid cross-section')

figure('Color', 'w')
stem(mode_extractf, amp_extractf, 'filled')
box on
grid on
xlabel('Mode n')
ylabel('Amplitude')
title('Extracted axisymmetric amplitudes')

function radiusOut = reconstructAxisymRadius(theta, ampExtract)
    theta = theta(:);
    ampExtract = ampExtract(:);

    maxModeLocal = numel(ampExtract) - 1;
    radiusOut = ampExtract(1) .* ones(size(theta));
    x = cos(theta);

    for n = 1:maxModeLocal
        radiusOut = radiusOut + ampExtract(n + 1) .* ...
            extract_real_spherical_harmonic(x, 0, n, 0);
    end
end
