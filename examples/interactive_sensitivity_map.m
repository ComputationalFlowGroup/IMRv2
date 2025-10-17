function interactive_sensitivity_map(metric_map, G_ratio_vals, Lambda_vals, l1_vals)
    % Create figure and axes
    fig = figure('Name', 'Interactive Sensitivity Map');
    ax = axes('Parent', fig);
    
    % Initial slice index for l1
    l1_index = 1;
    
    % Plot initial surface (G_ratio vs Lambda at fixed l1)
    [LambdaGrid, GratioGrid] = meshgrid(Lambda_vals, G_ratio_vals);
    Z = squeeze(metric_map(:, :, l1_index));
    hSurf = surf(ax, LambdaGrid, GratioGrid, Z);
    xlabel(ax, '\Lambda');
    ylabel(ax, 'G_1/G_0');
    zlabel(ax, 'Sensitivity');
    title(ax, ['Sensitivity at \ell_1/R_0 = ', num2str(l1_vals(l1_index))]);
    colorbar;
    shading interp;
    
    % Add slider to select l1 slice
    uicontrol('Style', 'text', 'Position', [150 10 200 20], 'String', 'Select \ell_1/R_0');
    slider = uicontrol('Style', 'slider', ...
        'Min', 1, 'Max', length(l1_vals), 'Value', l1_index, ...
        'SliderStep', [1/(length(l1_vals)-1) , 1/(length(l1_vals)-1)], ...
        'Position', [150 40 300 20], ...
        'Callback', @(src,~) slider_callback(src, ax, metric_map, G_ratio_vals, Lambda_vals, l1_vals, hSurf));
    
end

function slider_callback(src, ax, metric_map, G_ratio_vals, Lambda_vals, l1_vals, hSurf)
    % Get the slider value and round to nearest integer index
    idx = round(src.Value);
    
    % Update surface data
    Z = squeeze(metric_map(:, :, idx));
    [LambdaGrid, GratioGrid] = meshgrid(Lambda_vals, G_ratio_vals);
    
    set(hSurf, 'ZData', Z);
    set(hSurf, 'XData', LambdaGrid);
    set(hSurf, 'YData', GratioGrid);
    
    title(ax, ['Sensitivity at \ell_1/R_0 = ', num2str(l1_vals(idx))]);
    drawnow;
end