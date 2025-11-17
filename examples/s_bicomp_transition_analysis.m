%% Regime Transition Criteria for S_total
% Analytical identification of transition boundaries

clear; clc; close all;

%% Define symbolic variables
syms Lambda ell1_sym G1G0 G0_sym real positive

% Lambda1 as function of Lambda and ell1
Lambda1_sym = (1 + (Lambda^3 - 1)/ell1_sym^3)^(1/3);

% Stress components
S_near_sym = (G0_sym/2) * (1/Lambda^4 + 4/Lambda - 1/Lambda1_sym^4 - 4/Lambda1_sym);
S_far_sym = (G0_sym*G1G0/2) * (1/Lambda1_sym^4 + 4/Lambda1_sym - 5);
S_total_sym = S_near_sym + S_far_sym;

% Compute partial derivatives
dS_dLambda = diff(S_total_sym, Lambda);
dS_dell1 = diff(S_total_sym, ell1_sym);
dS_dG1G0 = diff(S_total_sym, G1G0);

% Simplify
dS_dLambda = simplify(dS_dLambda);
dS_dell1 = simplify(dS_dell1);
dS_dG1G0 = simplify(dS_dG1G0);

fprintf('=== ANALYTICAL DERIVATIVES ===\n\n');
fprintf('∂S/∂Λ:\n');
pretty(dS_dLambda);
fprintf('\n∂S/∂ℓ₁:\n');
pretty(dS_dell1);
fprintf('\n∂S/∂(G₁/G₀):\n');
pretty(dS_dG1G0);

%% Convert to numerical functions
dS_dLambda_func = matlabFunction(dS_dLambda, 'Vars', [G0_sym, G1G0, Lambda, ell1_sym]);
dS_dell1_func = matlabFunction(dS_dell1, 'Vars', [G0_sym, G1G0, Lambda, ell1_sym]);
dS_dG1G0_func = matlabFunction(dS_dG1G0, 'Vars', [G0_sym, G1G0, Lambda, ell1_sym]);
S_total_func = matlabFunction(S_total_sym, 'Vars', [G0_sym, G1G0, Lambda, ell1_sym]);

%% CRITERION 1: Where does ∂S/∂Λ = 0? (Critical stretch)
fprintf('\n\n=== CRITERION 1: Critical Stretch (∂S/∂Λ = 0) ===\n');

G0 = 1000; % Pa
G1G0_vals = [0.01, 0.1, .25, 0.5, .75, 1.0];
ell1_vals = logspace(0, 1, 20); % 1 to 10

Lambda_crit = zeros(length(ell1_vals), length(G1G0_vals));
Lambda_range = linspace(0.1, 10,100);

figure('Position', [100, 100, 1200, 800]);

for j = 1:length(G1G0_vals)
    subplot(3, 3, j);
    hold on;
    
    for i = 1:length(ell1_vals)
        % Evaluate derivative across Lambda range
        dS_vals = arrayfun(@(L) dS_dLambda_func(G0, G1G0_vals(j), L, ell1_vals(i)), Lambda_range);
        
        % Find zero crossings
        zero_crossings = find(diff(sign(dS_vals)));
        
        if ~isempty(zero_crossings)
            for zc = zero_crossings'
                % Refine with fzero
                try
                    L_crit = fzero(@(L) dS_dLambda_func(G0, G1G0_vals(j), L, ell1_vals(i)), ...
                                   [Lambda_range(zc), Lambda_range(zc+1)]);
                    Lambda_crit(i, j) = L_crit;
                    plot(ell1_vals(i), L_crit, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
                catch
                    Lambda_crit(i, j) = NaN;
                end
            end
        else
            Lambda_crit(i, j) = NaN;
        end
    end
    
    xlabel('$\ell_1$','Interpreter','latex');
    ylabel('\Lambda_{critical}');
    title(sprintf('G_1/G_0 = %.2f', G1G0_vals(j)));
    grid on;
    set(gca, 'XScale', 'log');
end
sgtitle('Critical Stretch: Where ∂S/∂Λ = 0');
%saveas(gcf,'./crit1_stretch_Sbi','png')
%% dry version of code above for dS/dLambda, dS/dell1, and dS/dalpha
G1G0_vals_sets = {
    [0.01, 0.1, 0.5, 1.0, 5, 10], ... % for dS/dLambda
    [0.01, 0.1, 0.5, 1.0, 5, 10], ...      % for dS/dell1
    [0.01, 0.1, 0.5, 1.0, 5, 10]           % for dS/dalpha
};
ell1_vals = logspace(0, 1, 20); % 1 to 10
Lambda_range = linspace(0.1, 10, 100);

% Function handles for derivatives
deriv_funcs = {dS_dLambda_func, dS_dell1_func, dS_dG1G0_func};
titles = {
    'Critical Stretch: Where \partialS/\partial\Lambda = 0', ...
    'Critical Stretch: Where \partialS/\partial\ell_1 = 0', ...
    'Critical Stretch: Where \partialS/\partial\alpha = 0'
};

%%Plotting loop
for k = 1:3
    G1G0_vals = G1G0_vals_sets{k};
    Lambda_crit = zeros(length(ell1_vals), length(G1G0_vals));
    
    figure('Position', [100, 100, 1200, 800]);
    
    for j = 1:length(G1G0_vals)
        subplot(3, 3, j); hold on;
        
        for i = 1:length(ell1_vals)
            % Evaluate derivative across Lambda range
            dS_vals = arrayfun(@(L) deriv_funcs{k}(G0, G1G0_vals(j), L, ell1_vals(i)), Lambda_range);
            
            % Find zero crossings
            zero_crossings = find(diff(sign(dS_vals)));
            
            if ~isempty(zero_crossings)
                for zc = zero_crossings'
                    try
                        L_crit = fzero(@(L) deriv_funcs{k}(G0, G1G0_vals(j), L, ell1_vals(i)), ...
                                       [Lambda_range(zc), Lambda_range(zc+1)]);
                        Lambda_crit(i, j) = L_crit;
                        plot(ell1_vals(i), L_crit, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
                    catch
                        Lambda_crit(i, j) = NaN;
                    end
                end
            else
                Lambda_crit(i, j) = NaN;
            end
        end
        
        xlabel('$\ell_1$', 'Interpreter', 'latex');
        ylabel('\Lambda_{critical}');
        title(sprintf('G_1/G_0 = %.2f', G1G0_vals(j)));
        grid on;
        set(gca, 'XScale', 'log');
    end
    
    sgtitle(titles{k});
end


%% CRITERION 2: Where does the sensitivity change most? (∂²S/∂Λ² = 0)
fprintf('\n=== CRITERION 2: Inflection Points (∂²S/∂Λ² = 0) ===\n');

d2S_dLambda2 = diff(dS_dLambda, Lambda);
d2S_dLambda2_func = matlabFunction(d2S_dLambda2, 'Vars', [G0_sym, G1G0, Lambda, ell1_sym]);

Lambda_inflection = zeros(length(ell1_vals), length(G1G0_vals));

figure('Position', [100, 100, 1200, 800]);

for j = 1:length(G1G0_vals)
    subplot(3, 3, j);
    hold on;
    
    for i = 1:length(ell1_vals)
        % Find inflection points
        d2S_vals = arrayfun(@(L) d2S_dLambda2_func(G0, G1G0_vals(j), L, ell1_vals(i)), Lambda_range);
        
        zero_crossings = find(diff(sign(d2S_vals)));
        
        if ~isempty(zero_crossings)
            for zc = zero_crossings'
                try
                    L_infl = fzero(@(L) d2S_dLambda2_func(G0, G1G0_vals(j), L, ell1_vals(i)), ...
                                   [Lambda_range(zc), Lambda_range(zc+1)]);
                    Lambda_inflection(i, j) = L_infl;
                    plot(ell1_vals(i), L_infl, 'bs', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
                catch
                    Lambda_inflection(i, j) = NaN;
                end
            end
        else
            Lambda_inflection(i, j) = NaN;
        end
    end
    
    xlabel('$\ell_1$','Interpreter','latex');
    ylabel('\Lambda_{inflection}');
    title(sprintf('G_1/G_0 = %.2f', G1G0_vals(j)));
    grid on;
    set(gca, 'XScale', 'log');
end
sgtitle('Inflection Points: Where ∂²S/∂Λ² = 0');
saveas(gcf,'./crit2_inflection','png')

%% CRITERION 3: Relative contribution analysis
% When does S_far dominate over S_near?
fprintf('\n=== CRITERION 3: Far-field Dominance (S_far/S_near = 1) ===\n');

S_near_func = matlabFunction(S_near_sym, 'Vars', [G0_sym, G1G0, Lambda, ell1_sym]);
S_far_func = matlabFunction(S_far_sym, 'Vars', [G0_sym, G1G0, Lambda, ell1_sym]);

figure('Position', [100, 100, 1400, 500]);

for j = 1:length(G1G0_vals)
    subplot(1, 6, j);
    
    % Create grid for contour plot
    [Lambda_grid, ell1_grid] = meshgrid(linspace(0.5, 10, 100), logspace(0, 1, 100));
    ratio_grid = zeros(size(Lambda_grid));
    
    for ii = 1:size(Lambda_grid, 1)
        for jj = 1:size(Lambda_grid, 2)
            S_near_val = S_near_func(G0, G1G0_vals(j), Lambda_grid(ii,jj), ell1_grid(ii,jj));
            S_far_val = S_far_func(G0, G1G0_vals(j), Lambda_grid(ii,jj), ell1_grid(ii,jj));
            
            if abs(S_near_val) > 1e-10
                ratio_grid(ii,jj) = abs(S_far_val / S_near_val);
            else
                ratio_grid(ii,jj) = NaN;
            end
        end
    end
    
    % Plot contours
    contourf(Lambda_grid, ell1_grid, log10(ratio_grid), 20);
    hold on;
    contour(Lambda_grid, ell1_grid, ratio_grid, [1, 1], 'r-', 'LineWidth', 3);
    
    xlabel('\Lambda');
    ylabel('$\ell_1$','Interpreter','latex');
    title(sprintf('G_1/G_0 = %.2f', G1G0_vals(j)));
    colorbar;
    set(gca, 'YScale', 'log');
    caxis([-2, 2]);
    title(colorbar, 'log_{10}(|S_{far}/S_{near}|)');
end
sgtitle('Far-field Dominance: Red line where S_{far} = S_{near}');
saveas(gcf,'./crit3_Sfar_dominate','png')

%% CRITERION 4: Dimensionless transition parameter
% Define: Π = (Λ₁/Λ) * (G₁/G₀)
fprintf('\n=== CRITERION 4: Dimensionless Transition Parameter ===\n');

figure('Position', [100, 100, 1200, 800]);

G1G0_range = logspace(-1, 0, 30);
Lambda_vals = [0.5, 1.0, 5.0, 10];
S_all = [];

for idx = 1:length(Lambda_vals)
    subplot(2, 2, idx);
    hold on;
    
    for j = 1:length(G1G0_range)
        for i = 1:length(ell1_vals)
            Lambda1 = (1 + (Lambda_vals(idx)^3 - 1)/ell1_vals(i)^3)^(1/3);
            Pi = (Lambda1/Lambda_vals(idx)) * G1G0_range(j);
            
            S_val = S_total_func(G0, G1G0_range(j), Lambda_vals(idx), ell1_vals(i));
            S_normalized = S_val / G0;
            S_all(end+1) = log10(abs(S_normalized) + eps);
            
            scatter(Pi, ell1_vals(i), 50, log10(abs(S_normalized)+eps), 'filled');
            caxis(clim);
        end
    end
    
    xlabel('\Pi = (\Lambda_1/\Lambda) \cdot (G_1/G_0)');
    ylabel('$\ell_1$','Interpreter','latex');
    title(sprintf('\\Lambda = %.2f', Lambda_vals(idx)));
    %colorbar;

    set(gca, 'YScale', 'log');
    grid on;
end
clim = [min(S_all),max(S_all)];
h = colorbar('Position',[0.92 01.1 0.02 0.8]);
title(h,'log_{10}(|S|/G_0)')
sgtitle('Collapse via Dimensionless Parameter Π');

%% CRITERION 5: Gradient magnitude map
fprintf('\n=== CRITERION 5: Sensitivity/Gradient Magnitude Map ===\n');

Lambda_test = logspace(log10(0.3), log10(10), 40);
ell1_test = logspace(0, 1, 40);
G1G0_test = logspace(-1, 0, 20);

[G1G0_grid, Lambda_grid, ell1_grid] = ndgrid(G1G0_test, Lambda_test, ell1_test);

% Compute gradient magnitude
grad_mag = zeros(size(G1G0_grid));

for i = 1:numel(G1G0_grid)
    dL = dS_dLambda_func(G0, G1G0_grid(i), Lambda_grid(i), ell1_grid(i));
    de = dS_dell1_func(G0, G1G0_grid(i), Lambda_grid(i), ell1_grid(i));
    dG = dS_dG1G0_func(G0, G1G0_grid(i), Lambda_grid(i), ell1_grid(i));
    grad_mag(i) = sqrt(dL^2 + de^2 + dG^2);
end

% Find top 5% highest gradients (transition regions)
threshold = prctile(grad_mag(:), 95);
high_grad_idx = grad_mag > threshold;

figure('Position', [100, 100, 800, 600]);
scatter3(G1G0_grid(high_grad_idx), Lambda_grid(high_grad_idx), ell1_grid(high_grad_idx), ...
         50, log10(grad_mag(high_grad_idx)), 'filled');
xlabel('G_1/G_0');
ylabel('\Lambda');
zlabel('\ell_1');
set(gca, 'YScale', 'log', 'ZScale', 'log');
colorbar;
title('High Sensitivity Regions (Top 5%)');
view(45, 30);

%% Summary of Transition Criteria
fprintf('\n\n========== SUMMARY OF TRANSITION CRITERIA ==========\n');
fprintf('1. CRITICAL STRETCH: ∂S/∂Λ = 0\n');
fprintf('   → Peak stress locations\n\n');
fprintf('2. INFLECTION POINTS: ∂²S/∂Λ² = 0\n');
fprintf('   → Maximum rate of change\n\n');
fprintf('3. FAR-FIELD DOMINANCE: |S_far/S_near| = 1\n');
fprintf('   → Crossover between near and far field contributions\n\n');
fprintf('4. DIMENSIONLESS PARAMETER: Π = (Λ₁/Λ)·(G₁/G₀)\n');
fprintf('   → Universal scaling parameter\n\n');
fprintf('5. GRADIENT MAGNITUDE: ||∇S|| > threshold\n');
fprintf('   → Regime boundaries (sharpest transitions)\n');
fprintf('====================================================\n');
