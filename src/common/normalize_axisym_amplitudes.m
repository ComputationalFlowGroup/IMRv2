function epn = normalize_axisym_amplitudes(amp, radius)
%NORMALIZE_AXISYM_AMPLITUDES Normalize signed Y_n^0 coefficients by radius.
%   The first row is kept as the dimensional radius. Rows 2:end are
%   returned as dimensionless signed perturbation amplitudes.

    if nargin < 2 || isempty(radius)
        radius = amp(1, :);
    end

    validateattributes(amp, {'numeric'}, {'2d'}, mfilename, 'amp');

    radius = radius(:).';
    if numel(radius) ~= size(amp, 2)
        error('%s:RadiusSizeMismatch', mfilename, ...
            'radius must have one value per frame.');
    end

    epn = nan(size(amp));
    epn(1, :) = radius;

    finite_radius = radius(isfinite(radius));
    if isempty(finite_radius)
        radius_scale = 1;
    else
        radius_scale = max(abs(finite_radius));
    end
    radius_tol = max(eps, eps(radius_scale));
    valid_radius = isfinite(radius) & abs(radius) > radius_tol;
    epn(2:end, valid_radius) = amp(2:end, valid_radius) ./ radius(valid_radius);
end
