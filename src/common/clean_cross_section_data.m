function [theta, radius] = clean_cross_section_data(simdata, min_points)
%CLEAN_CROSS_SECTION_DATA Validate and sort two-column cross-section data.

    if nargin < 2
        min_points = 2;
    end

    validateattributes(simdata, {'numeric'}, {'2d', 'ncols', 2}, ...
        mfilename, 'simdata');

    finite_rows = all(isfinite(simdata), 2);
    simdata = simdata(finite_rows, :);

    if size(simdata, 1) < min_points
        error('%s:NotEnoughData', mfilename, ...
            'Need at least %d finite points; received %d.', ...
            min_points, size(simdata, 1));
    end

    [theta, order] = sort(simdata(:, 1));
    radius = simdata(order, 2);

    theta = theta(:);
    radius = radius(:);
end
