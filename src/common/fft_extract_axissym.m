function [mode_extractf, amp_extractf, phase_extractf] = fft_extract_axissym(simdata, max_mode)
%FFT_EXTRACT_AXISYM Extract axisymmetric spherical-harmonic amplitudes.
%   simdata(:,1) is theta and simdata(:,2) is the measured radius. The
%   fitted signal is
%
%       r(theta) = a_0 + sum_{n=1}^{max_mode} a_n Y_n^0(cos(theta)).
%
%   a_0 is returned as the mean-radius offset, not as the coefficient of
%   Y_0^0. This preserves the convention used by the processing scripts.

    if nargin < 2 || isempty(max_mode)
        max_mode = 20;
    end

    validateattributes(max_mode, {'numeric'}, ...
        {'scalar', 'integer', 'nonnegative'}, mfilename, 'max_mode');

    [theta, radius] = clean_cross_section_data(simdata, max_mode + 1);

    basis = ones(numel(theta), max_mode + 1);
    x = cos(theta);
    for n = 1:max_mode
        basis(:, n + 1) = extract_real_spherical_harmonic(x, 0, n, 0);
    end

    coeff = basis \ radius;

    mode_extractf = (0:max_mode).';
    amp_extractf = coeff(:);
    phase_extractf = zeros(max_mode + 1, 1);
    phase_extractf(amp_extractf < 0) = pi;
end
