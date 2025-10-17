%% composite mat: 2 homogeneous domains
Rst = linspace(0.2,10,200); %R_max/R_0
G0 = 500;
G1 = 1000;
l1 = 1.5; %nondim
l2 = 3; %nondim
a = 2.5;
n = 0.3;
b = 0.5;
%%
% goal: what value of l1 are the total stress contributions from each region equal
G0 = 1000; G1=1000;
% Parameters
Lambda = linspace(1,1.5,10);      % Rmax / R0
%ell1 = logspace(1, 4, 500);  % Dimensionless inner boundary
ell1 = linspace(1,10,100000);
cmap=parula(round(1.25*length(Lambda)));
figure
hold on;
for i=1:length(Lambda)
%Lambda1 = (1 + ((Lambda(i)^3 - 1)./ell1.^3)).^(1/3);
Lambda1 = (ell1.^3 + Lambda(i)^3 - 1).^(1/3);
% Compute stress contributions
S_near = abs((G0/2) * (1./Lambda(i)^4 + 4/Lambda(i) - (1./Lambda1.^4 + 4./Lambda1)));
S_far  = abs((G1/2) * (1./Lambda1.^4 + 4./Lambda1 - 5));
S_total = S_near + S_far;

% Fractional contributions
%S_frac_near = S_near ./ S_total;
S_frac_far  = S_far  ./ S_total;

% % Find where they are approximately equal
% [~, idx_equal] = min(abs(S_frac_near - S_frac_far));
% ell1_equal = ell1(idx_equal);

% Plot
%plot(ell1, S_frac_near, 'b-', 'LineWidth', 2);
plot(ell1, S_frac_far, 'Color',cmap(i,:), 'LineWidth', 2);
yline(.5, 'k:', 'LineWidth', 2);
%legend({'Near-field', 'Far-field', 'Equal contribution'}, 'Location', 'best');
xlabel('$\ell_1 = l_1 / R_0$', 'Interpreter','latex','FontSize',20);
ylabel('$S_{\rm far} / S_{\rm total}$', 'Interpreter','latex','FontSize',20);
grid on; box on;
set(gcf,'color','w');
set(gca, 'FontName', 'Times', 'FontSize',20); 
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.03 .03];
xa.LineWidth = 1.5;
xticks(0:2:10);
end
hold off;
%saveas(gcf,'Sfracfar_G0eqG1','png')
%% for different Lambda values
% Parameters
Lambda_vals = linspace(1.01, 10, 100);
ell1_vals = linspace(1.01, 10, 200);

% Initialize matrix to store crossover ell1 for each Lambda
ell1_crossover = NaN(size(Lambda_vals));

for i = 1:length(Lambda_vals)
    Lambda = Lambda_vals(i);
    
    Lambda1 = nthroot(1 + (Lambda^3 - 1)./ell1_vals.^3, 3);
    
    % Stress contributions
    S_near = (G0/2) * (1/Lambda^4 + 4/Lambda - (1./Lambda1.^4 + 4./Lambda1));
    S_far  = (G1/2) * (1./Lambda1.^4 + 4./Lambda1 - 5);
    S_total = S_near + S_far;
    
    % Fractional contributions
    S_frac_near = S_near ./ S_total;
    S_frac_far  = S_far  ./ S_total;
    
    % Find crossover by minimal difference
    [~, idx_equal] = min(abs(S_frac_near - S_frac_far));
    ell1_crossover(i) = ell1_vals(idx_equal);
end

figure;
plot(Lambda_vals, ell1_crossover, 'LineWidth', 2)
xlabel('$\Lambda$', 'Interpreter','latex','FontSize',16)
ylabel('$\ell_1$ crossover', 'Interpreter','latex','FontSize',16)
title('Crossover point where Near-field = Far-field stress', 'Interpreter','latex')
grid on; box on;
%%
% 
% % Numerical differentiation - finite difference approx of d(ell1)/d(Lambda)
% dLambda = diff(Lambda_vals);
% dell1 = diff(ell1_crossover);
% 
% slope = dell1 ./ dLambda;  % slope at midpoints between Lambda_vals
% 
% % To plot slope vs Lambda, assign slope to midpoints of Lambda_vals
% Lambda_mid = (Lambda_vals(1:end-1) + Lambda_vals(2:end))/2;
% 
% % Plot slope
% figure;
% plot(Lambda_mid, slope, 'LineWidth', 2);
% xlabel('$\Lambda$', 'Interpreter','latex', 'FontSize', 16);
% ylabel('Slope $d\ell_1/d\Lambda$', 'Interpreter','latex', 'FontSize', 16);
% title('Slope of crossover curve $\ell_1(\Lambda)$', 'Interpreter','latex');
% grid on; box on;
% Assuming Lambda_vals and ell1_crossover are vectors of the same length

p = polyfit(Lambda_vals, ell1_crossover, 1);  % 1 = linear fit

slope_reg = p(1);  % slope of the fitted line
intercept_reg = p(2);

% Plot the data and the regression line
figure;
plot(Lambda_vals, ell1_crossover, 'bo', 'MarkerSize', 6, 'DisplayName', 'Data'); hold on;
plot(Lambda_vals, polyval(p, Lambda_vals), 'r-', 'LineWidth', 2, 'DisplayName', 'Linear fit');
xlabel('$\Lambda$', 'Interpreter','latex', 'FontSize', 16);
ylabel('$\ell_1$ crossover', 'Interpreter','latex', 'FontSize', 16);
title('Linear regression of crossover length vs Lambda', 'Interpreter','latex');
legend('Location', 'best');
grid on; box on;

fprintf('Regression slope = %.4f\n', slope_reg);

%%

% Precompute near/far field stress dominance
S_near_dominates = false(length(l1_range), length(Lambda_vals));

for i = 1:length(Lambda_vals)
    Rst = Lambda_vals(i);
    diffs = NaN(size(l1_range));
    
    for j = 1:length(l1_range)
        l1 = l1_range(j);
        Lambda1 = (1 + (Rst^3 - 1)/l1^3)^(1/3);
        
        % Skip nonphysical
        if ~isreal(Lambda1) || Lambda1 < 1
            continue
        end

        S_near = (G0/2) * (1/Rst^4 + 4/Rst - (1/Lambda1^4 + 4/Lambda1));
        S_far = (G1/2) * (1/Lambda1^4 + 4/Lambda1 - 5);

        diffs(j) = abs(S_near - S_far);

        % For shading: track if near > far
        S_near_dominates(j, i) = (S_near > S_far);
    end

    % Find crossover
    [~, min_idx] = min(diffs);
    crossover_l1(i) = l1_range(min_idx);
end

% Plot shaded dominance regions
figure;
hold on;

% Create shaded image
imagesc(Lambda_vals, l1_range, S_near_dominates);
colormap([0.8 0.9 1; 1 0.9 0.9]); % blue = far field dominates, red = near field
set(gca,'YDir','normal');

% Overlay crossover curve
plot(Lambda_vals, crossover_l1, 'k-', 'LineWidth', 2);

% Labels and formatting
xlabel('$\Lambda$ (Stretch)', 'Interpreter','latex', 'FontSize', 18);
ylabel('$\ell_1$', 'Interpreter','latex', 'FontSize', 18);
title('Stress Dominance Map (Near vs Far Field)', 'Interpreter','latex');

legend({'$S_{\mathrm{near}} = S_{\mathrm{far}}$'}, 'Interpreter','latex', 'Location','northeast');

% Custom colorbar
cb = colorbar('Ticks',[0.25 0.75], 'TickLabels',{'Far Field','Near Field'});
ylabel(cb, 'Dominant Stress Source', 'Interpreter','latex');
set(gca, 'FontSize', 14, 'TickLabelInterpreter','latex');
box on;

%%
% Parameters
G0 = 1;       % Shear modulus near field
G1 = 2;       % Shear modulus far field

% Domain for stretch and near-field thickness
Lambda = linspace(1.01, 10, 200);    % R_max / R0
l1 = linspace(0.1, 6, 200);          % l1 / R0

% Meshgrid
[LambdaGrid, l1Grid] = meshgrid(Lambda, l1);

% Compute Lambda1 from incompressibility
Lambda1 = (1 + (LambdaGrid.^3 - 1) ./ l1Grid.^3).^(1/3);
Lambda2 = LambdaGrid;

% Compute stress contributions
S_near = (G0/2) .* (1 ./ LambdaGrid.^4 + 4 ./ LambdaGrid - (1 ./ Lambda1.^4 + 4 ./ Lambda1));
S_far = (G1/2) .* (1 ./ Lambda2.^4 + 4 ./ Lambda2 - 5);
S_total = S_near + S_far;

% Compute imbalance ratio
imbalance_ratio = (S_near - S_far) ./ S_total;
imbalance_ratio(~isfinite(imbalance_ratio)) = NaN;


% Plotting the imbalance map
figure;
contourf(LambdaGrid, l1Grid, imbalance_ratio, 100, 'LineColor', 'none');
hold on;
contour(LambdaGrid, l1Grid, imbalance_ratio, [0 0], 'k', 'LineWidth', 2);  % Equal contribution curve
clim([-50 50]);
colorbar;
colormap(turbo(200));
xlabel('$\Lambda$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\ell_1$', 'Interpreter', 'latex', 'FontSize', 18);
title('Stress Imbalance Map: $(S_{\mathrm{near}} - S_{\mathrm{far}})/S_{\mathrm{total}}$', ...
    'Interpreter', 'latex', 'FontSize', 16);
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');

%% material contrast sensitivity map, different G1/G0 ratios
% Parameters
Lambda_vals = linspace(1.001, 10, 100);    % Rmax/R0 (Λ)
ell1_vals   = linspace(1.001, 10, 300);  % Dimensionless l1
G_ratios    = [0.01, 0.1, 0.5, 1, 2, 5, 10, 100];  % G1/G0 values

% Preallocate
L1_crossover = NaN(size(G_ratios));
all_crossovers = NaN(length(G_ratios), length(Lambda_vals));

figure; hold on;
colors = lines(length(G_ratios));

for k = 1:length(G_ratios)
    G0 = 1;
    G1 = G_ratios(k) * G0;
    
    crossover_l1 = NaN(size(Lambda_vals));
    
    for j = 1:length(Lambda_vals)
        Lambda = Lambda_vals(j);
        ell1 = ell1_vals;
        
        Lambda1 = nthroot(1 + (Lambda^3 - 1) ./ ell1.^3, 3);
        %Lambda1 = ell1 / Lambda; %linear approx

        % if G0 == G1
        %    S_total = (G0/2) *(1/Lambda^(4) + 4/Lambda -5);
        %    S_near = ((ell1.^3 - 1) ./ (Lambda^3 - 1)) .* S_total;
        %    S_far = ((Lambda^3 - ell1.^3) ./ (Lambda^3 - 1)) .* S_total;
        %    crossover_l1(j) = ((Lambda^3 + 1)/2)^(1/3); %equal volume (energy) nsplit in homogeneous case
        % else
        
        S_near = (G0/2) * (1/Lambda^4 + 4/Lambda - (1./Lambda1.^4 + 4./Lambda1));
        S_far  = (G1/2) * (1./Lambda1.^4 + 4./Lambda1 - 5);
        S_total = S_near + S_far;
        
        % Fractional contributions
        S_frac_near = S_near ./ S_total;
        S_frac_far = S_far ./ S_total;

        if G0 == G1
            S_homo = (G0/2) *(1/Lambda^(4) + 4/Lambda -5);
            remainder = S_homo - S_total;
            if remainder > 1e-14; stop; end
        end
        
        % Find crossover point (where fractions ~ equal)
        [~, idx_eq] = min(abs(S_frac_near - S_frac_far));
        crossover_l1(j) = ell1(idx_eq);
        % end
    end
    
    % Store and plot
    all_crossovers(k, :) = crossover_l1;
    plot(Lambda_vals, crossover_l1, 'LineWidth', 2, 'Color', colors(k,:), ...
        'DisplayName', ['G_1/G_0 = ' num2str(G_ratios(k))]);
end
% % Numerical slope (from your computed crossover)
%numerical_slope = diff(crossover_l1) ./ diff(Lambda_vals);

% Extract numerical crossover for G1 = G0
idx_equal = find(G_ratios == 1);
crossover_equal = all_crossovers(idx_equal, :);

% Plot
figure;
plot(Lambda_vals, crossover_equal, 'bo-', 'DisplayName', 'Numerical (G_0 = G_1)');
hold on;
plot(Lambda_vals, ell1_analytical, 'r--', 'LineWidth', 2, 'DisplayName', 'Analytical (Energy Theory)');
xlabel('$\Lambda = R_{\mathrm{max}} / R_0$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\ell_1^{\mathrm{crossover}}$', 'Interpreter', 'latex', 'FontSize', 18);
title('Crossover $\ell_1$ vs $\Lambda$ for Varying Material Contrast', 'Interpreter', 'latex');
legend('Location', 'northwestoutside');
grid on; box on;
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');

idx_equal = find(G_ratios == 1);
% Extract Lambda_vals and crossover_l1 for G0=G1 case
Lambda_equal = Lambda_vals;
crossover_equal = all_crossovers(idx_equal, :);

% Perform linear regression (polyfit) for the equal case
p = polyfit(Lambda_equal, crossover_equal, 1);  % Linear fit
slope_reg = p(1);
intercept_reg = p(2);

% Plot data and regression line for G0=G1 case only
figure;
plot(Lambda_equal, crossover_equal, 'bo', 'MarkerSize', 6, 'DisplayName', 'Data (G_0=G_1)');
hold on;
plot(Lambda_equal, polyval(p, Lambda_equal), 'r-', 'LineWidth', 2, 'DisplayName', 'Linear fit');
xlabel('$\Lambda$', 'Interpreter','latex', 'FontSize', 16);
ylabel('$\ell_1$ crossover', 'Interpreter','latex', 'FontSize', 16);
title('Linear regression of crossover length vs Lambda $(G_0=G_1)$', 'Interpreter','latex');
legend('Location', 'best');
grid on; box on;

% Print slope
fprintf('Regression slope $(G0=G1)$ = %.4f\n', slope_reg,'Interpreter','latex');

%
% Numerical slope from your data
numerical_slope = diff(crossover_equal) ./ diff(Lambda_vals);
Lambda_mid = (Lambda_vals(1:end-1) + Lambda_vals(2:end)) / 2;

% Plot
figure;
plot(Lambda_mid, numerical_slope, 'b-', 'LineWidth', 2, 'DisplayName', 'Numerical');
hold on;
plot(Lambda_vals, ell1_slope, 'r--', 'LineWidth', 2, 'DisplayName', 'Analytical (Implicit d\ell_1/d\Lambda)');
xlabel('$\Lambda$', 'Interpreter', 'latex');
ylabel('$d\ell_1/d\Lambda$', 'Interpreter', 'latex');
legend('Location', 'best');
title('Slope Comparison: Numerical vs Analytical');
grid on; box on;


%% stress field profile - 2composite
% Parameters
Lambda = 3;         % Outer boundary stretch (R_max/R_0)
G0 = 1;             % Inner modulus
G1 = 5;             % Outer modulus
ell1 = 1.1;           % Interface location in real space (dimensionless l1)

% Compute the switch point in stretch space
Lambda1 = nthroot(1 + (Lambda^3 - 1)/ell1^3, 3);

% Lambda grid
lambda_vals = linspace(1, Lambda, 500);
S_field = zeros(size(lambda_vals));

% Compute the stress field with modulus jump at Lambda1 (where switch happens)
for i = 1:length(lambda_vals)
    lam = lambda_vals(i);
    if lam <= Lambda1
        G = G0;
    else
        G = G1;
    end
    S_field(i) = (G/2) * (lam^-4 + 4*lam^-1 - 5); % instantaneous value of stress at each point
    S_near = (G0/2) * (1/Lambda^4 + 4/Lambda - (1/Lambda1.^4 + 4/Lambda1));
    S_far = (G1/2) * (1/Lambda1.^4 + 4/Lambda1 -5);
    S_total = S_near + S_far;
    S_frac_near = S_near ./ S_total;
    S_frac_far = S_far ./ S_total;
end

% Plot
figure;
plot(lambda_vals, S_field, 'k', 'LineWidth', 2); hold on;
xline(Lambda1, 'r--', 'LineWidth', 2, 'DisplayName', '$\ell_1$ interface');

xlabel('$\lambda = r/R_0$', 'Interpreter','latex', 'FontSize', 18);
ylabel('$S^{e}(\lambda)$', 'Interpreter','latex', 'FontSize', 18);
title('Radial Stress Field with Modulus Jump at $\ell_1$', 'Interpreter','latex');
legend('Location', 'best', 'Interpreter', 'latex');
grid on; box on;
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');

%% 2composite stress integral and taus
x1 = @(Rst) (1+(Rst.^3-1)./(l1).^3).^(1/3); %Lambda1
%x2 = @(Rst) (1+(Rst.^3-1)./(l2).^3).^(1/3); %Lambda2

S0 = (G0/2)*(1./Rst.^4 + 4./Rst - (1./x1(Rst).^4 + 4./x1(Rst)));
S1 = -(G1/2)*(5 - 4./x1(Rst) - 1./x1(Rst).^4);
S0_homo = -(G0/2)*(5 - 4./Rst - 1./Rst.^4);
S1_homo = -(G1/2)*(5 - 4./Rst - 1./Rst.^4);

tau = @(x) (2/3)*(1./x.^4 - x.^2);
gtau0 = zeros(size(Rst));
gtau1 = zeros(size(Rst));
for i = 1:length(Rst)
    rst = Rst(i);
    gtau0(i) = G0*tau(rst);
    gtau1(i) = G1*tau(rst);
end

figure
hold on;
plot(Rst,S0_homo/G0,'r','LineWidth',3)
plot(Rst,S1_homo/G0,'k--','LineWidth',3)
plot(Rst,(S0+S1)/G0,'-.g','LineWidth',3)
ylim([-5 5])
xlabel('$R_{\mathrm{max}}/R_{0}$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$S/G_0$','Interpreter','Latex','FontSize',24);
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.015 .015];
xa.LineWidth = 1.5;
tickrange= 0:2:10;
xticks(tickrange)
tickrange= -5:2:5;
yticks(tickrange)
box on;
hold off;

figure
hold on;
plot(Rst,gtau0,'r','LineWidth',3)
plot(Rst,gtau1,'k--','LineWidth',3)
plot(Rst,gtau0+gtau1,'-.g','LineWidth',3)
ylim([-5 5])
xlabel('$R_{\mathrm{max}}/R_{0}$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$\tau_{rr}/p_{\infty}$','Interpreter','Latex','FontSize',24);
set(gcf,'color','w');
set(gca,'FontName','Times','FontSize',20);
set(gca,'TickLabelInterpreter','latex')
xa = gca;
xa.TickLength = [.015 .015];
xa.LineWidth = 1.5;
tickrange= 0:2:10;
xticks(tickrange)
tickrange= -5:2:5;
yticks(tickrange)
box on;
hold off;

%%
%% Parameters
G0 = 1;           % Shear modulus near-field
G1 = 5;           % Shear modulus far-field
a  = 2.5;         % Shape param (for C-Y)
n  = 0.3;         % Shape param (for C-Y)
b = 0.5;          % Shape param (for tanh)
Lambda = 3;       % Bubble expansion ratio
R0 = 1;           % Initial radius (normalized)
Rmax = Lambda * R0;

% Graded geometry
ell1 = 1.5;              % Start of graded region (dimless)
ell2 = 3;                % End of graded region (dimless)
r0 = linspace(R0, Rmax, 500);  % Radial positions

% Compute local stretches at each r0
lambda = (1 + (Lambda^3 - 1) ./ (r0 .^3)).^(1/3);

% Define graded modulus profile using tanh
% m = @(r) 0.5 * (1 + tanh(a * ( (r - ell1) / (ell2 - ell1) - 0.5 )));
% G_r = @(r) G0 + (G1 - G0) * m(r);

% Define f(lambda) based on theory
f_tanh = @(lambda) (ell2 + ell1)/(ell2 - ell1) * ...
          ( ((Lambda^3 - 1) ./ (lambda.^3 - 1)).^(1/3) + (ell1 + ell2)/2 );

% Graded modulus function using your tanh form
m_tanh = @(lambda) 0.5 * (1 + tanh(b * f_tanh(lambda)));

G_tanh = @(lambda) G0 + (G1 - G0) * m(lambda);  % Graded modulus


% Define CY shape function
f_cy = @(lambda) (ell2 .* ((lambda.^3 - 1) ./ (Lambda^3 - 1)).^(1/3) - 1) ./ ...
       (1 - ell1 .* ((lambda.^3 - 1) ./ (Lambda^3 - 1)).^(1/3));

% Graded modulus using CY model
G_cy = @(lambda) G0 + (G1 - G0) .* (1 + f_cy(lambda).^a).^((n - 1)/a);


% Compute local stress density: neo-Hookean
S_field_cy = G_cy(lambda) .* (lambda.^(-5) + lambda.^(-2));
S_field_tanh = G_tanh(lambda) .* (lambda.^(-5) + lambda.^(-2));

% Plot stress field
figure;
plot(lambda, S_field_cy, 'LineWidth', 2);
hold on;
% Highlight graded region
yL = ylim;
patch([ell1 ell2 ell2 ell1], [yL(1) yL(1) yL(2) yL(2)], ...
      [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

xlabel('$r_0 / R_0$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Stress density $S(\lambda)_{cy}$', 'Interpreter', 'latex', 'FontSize', 16);
title(['Stress Field with Graded Region: $\ell_1 = ', num2str(ell1), ...
       '$, $\ell_2 = ', num2str(ell2), '$'], 'Interpreter', 'latex');
legend('Stress', 'Graded Region', 'Location', 'best');
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
grid on;
hold off;

figure;
plot(lambda, G_cy, 'LineWidth', 2);
xlabel('$r_0/R_0$', 'Interpreter','latex', 'FontSize', 16);
ylabel('Stress density $S(\lambda)_{tanh}$', 'Interpreter', 'latex', 'FontSize', 16);
% Highlight graded region
yL = ylim;
patch([ell1 ell2 ell2 ell1], [yL(1) yL(1) yL(2) yL(2)], ...
      [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

xlabel('$r_0 / R_0$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Stress density $S(r_0)$', 'Interpreter', 'latex', 'FontSize', 16);
title(['Stress Field with Graded Region: $\ell_1 = ', num2str(ell1), ...
       '$, $\ell_2 = ', num2str(ell2), '$'], 'Interpreter', 'latex');
legend('Stress', 'Graded Region', 'Location', 'best');
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
grid on;
hold off;
grid on; box on;
set(gca, 'FontSize', 14, 'TickLabelInterpreter','latex');


ell1_vals = linspace(1.1, 3, 20);
ell2_vals = linspace(3.1, 6, 20);

[ELL1, ELL2] = meshgrid(ell1_vals, ell2_vals);
MAX_STRESS = zeros(size(ELL1));

for i = 1:numel(ELL1)
    e1 = ELL1(i);
    e2 = ELL2(i);
    
    m = @(r) 0.5 * (1 + tanh(a * ((r - e1) / (e2 - e1) - 0.5)));
    G_r0 = @(r) G0 + (G1 - G0) * m(r);
    S_field = G_r0(r0) .* (lambda.^(-5) + lambda.^(-2));
    
    MAX_STRESS(i) = max(S_field);
end

% Heatmap of max stress vs graded geometry
figure;
contourf(ELL1, ELL2, MAX_STRESS, 20, 'LineColor', 'none');
xlabel('$\ell_1$', 'Interpreter', 'latex');
ylabel('$\ell_2$', 'Interpreter', 'latex');
title('Max Stress vs Graded Region Geometry', 'Interpreter', 'latex');
colorbar;
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');


%% 3d plots: Se vs l1/l2 vs Rst
% option A: graded width / extent of mat transition
% option B: graded location and extent: where and how much grading occurs

% explore A (l1,l2 dimless already)
Lambda = linspace(0.2, 10, 100);
ratios = linspace(1.1, 5, 50);  % l2/l1
[LAM, RATIO] = meshgrid(Lambda, ratios);

% EXtent = 1;
% L1 = linspace(1.05,5,100);
% [LAM,L1GRID] = meshgrid(Lambda,L1);
% L2GRID = L1GRID + EXtent;

S_tanh = zeros(size(LAM));
S_ycy = zeros(size(LAM));
S_frac_gtanh = NaN(size(LAM));
S_frac_ntanh = NaN(size(LAM));
S_frac_ftanh = NaN(size(LAM));
S_frac_gcy = NaN(size(LAM));
S_frac_ncy = NaN(size(LAM));
S_frac_fcy = NaN(size(LAM));

for i = 1:numel(LAM)
    Rst = LAM(i);
    l2 = RATIO(i) * l1;  % l2 = ratio * l1

    % l1 = L1GRID(i);
    % l2 = L2GRID(i);
    
    Lambda1 = (1 + (Rst^3 - 1) / l1^3)^(1/3);
    Lambda2 = (1 + (Rst^3 - 1) / l2^3)^(1/3);

    % tanh
    f_tanh = @(x) (l2 + l1)/(l2 - l1) * ...
             (( (Rst^3 - 1)./(x.^3 - 1) ).^(1/3) + (l1 + l2)/2);
    y_tanh = @(x) (G0 + (G1 - G0) * 0.5*(1 + tanh(0.5 * f_tanh(x)))) .* (1./x.^5 + 1./x.^2);

    % cy
    fcy = @(x) (l2*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3) - 1)./...
      (1-l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3));
    ycy = @(x) (G0+(G1-G0)*(1+( fcy(x) ).^a).^((n-1)/a)).*(1./x.^5+1./x.^2);

    try
        if Lambda1 < Lambda2 && isreal(Lambda1) && isreal(Lambda2); continue; end
        % Stress integral
        Sg_ycy = 2 * integral(ycy, Lambda1, Lambda2, 'AbsTol',1e-8,'RelTol',1e-8);
        Sg_ytanh = 2 * integral(y_tanh, Lambda1, Lambda2, 'AbsTol',1e-8,'RelTol',1e-8);
    catch
        Sg_tanh = NaN; Sg_ycy = NaN;
    end
    S0 = G0/2 * (1/Rst^4 + 4/Rst - (1/Lambda1^4 + 4/Lambda1));
    S1 = G1/2 * (1/Lambda2^4 + 4/Lambda2 - 5);

    S_tanh(i) = (S0 + Sg_ytanh + S1) / G0;
    S_ycy(i) = (S0 + Sg_ycy + S1) / G0;

    % fractional contribution
    S_frac_gtanh(i) = Sg_ytanh / S_tanh(i);
    S_frac_ntanh(i) = S0 / S_tanh(i);
    S_frac_ftanh(i) = S1 / S_tanh(i);
    S_frac_gcy(i) = Sg_ycy / S_ycy(i);
    S_frac_ncy(i) = S0 / S_ycy(i);
    S_frac_fcy(i) = S1 / S_ycy(i);
end

figure;
surf(LAM, RATIO, S_tanh, 'EdgeColor','none')
xlabel('$\Lambda$', 'Interpreter','latex','FontSize',20)
ylabel('$\ell_1 / \ell_2$', 'Interpreter','latex','FontSize',20)
zlabel('$S/G_0$', 'Interpreter','latex','FontSize',20)
colorbar
view(45,30)

figure;
surf(LAM, RATIO, S_ycy, 'EdgeColor','none')
xlabel('$\Lambda$', 'Interpreter','latex','FontSize',20)
ylabel('$\ell_1 / \ell_2$', 'Interpreter','latex','FontSize',20)
zlabel('$S/G_0$', 'Interpreter','latex','FontSize',20)
colorbar
view(45,30)

% fractional contributions - dominance map
figure;
contourf(LAM, RATIO, S_frac_gtanh,20,'LineColor','none')
xlabel('$\Lambda$', 'Interpreter','latex','FontSize',20)
ylabel('$\ell_1 / \ell_2$', 'Interpreter','latex','FontSize',20)
cb = colorbar;
ylabel(cb,'$S_g/S_{tanh}$','Interpreter','latex')
set(gca,'FontSize',16,'TickLabelInterpreter','latex')
% Overlay regime transition contours
% hold on;
% [~, h1] = contour(LAM, RATIO, S_frac_gtanh, [0.5 0.5], 'r', 'LineWidth', 2);
% [~, h2] = contour(LAM, RATIO, S_frac_ntanh, [0.5 0.5], 'b--', 'LineWidth', 2);
% [~, h3] = contour(LAM, RATIO, S_frac_ftanh, [0.5 0.5], 'k-.', 'LineWidth', 2);
% %legend([h1(1) h2(1) h3(1)], {'Graded = 50%', 'Near-field = 50%', 'Far-field = 50%'}, ...
% %    'Interpreter','latex', 'FontSize',12, 'Location','southwest');
% hold off;


figure;
contourf(LAM, RATIO, S_frac_gcy,20,'LineColor','none')
% xlabel('$\Lambda$', 'Interpreter','latex','FontSize',20)
% ylabel('$\ell_1 / \ell_2$', 'Interpreter','latex','FontSize',20)
% cb = colorbar;
% ylabel(cb,'$S_g/S_{cy}$','Interpreter','latex')
% set(gca,'FontSize',16,'TickLabelInterpreter','latex')

% Overlay regime transition contours
% hold on;
% [~, h1] = contour(LAM, RATIO, S_frac_gcy, [0.25,0.5, 0.5], 'r', 'LineWidth', 2);
% [~, h2] = contour(LAM, RATIO, S_frac_ncy, [0.5 0.5], 'b--', 'LineWidth', 2);
% [~, h3] = contour(LAM, RATIO, S_frac_fcy, [0.5 0.5], 'k-.', 'LineWidth', 2);
%legend([h1(1) h2(1) h3(1)], {'Graded = 50%', 'Near-field = 50%', 'Far-field = 50%'}, ...
%    'Interpreter','latex', 'FontSize',12, 'Location','southwest');
xlabel('$\Lambda$', 'Interpreter','latex','FontSize',20)
ylabel('$\ell_1 / \ell_2$', 'Interpreter','latex','FontSize',20)
cb = colorbar;
ylabel(cb,'$S_g/S_{cy}$','Interpreter','latex')
set(gca,'FontSize',16,'TickLabelInterpreter','latex')
hold off;

% figure;
% contourf(LAM, L1GRID, S_frac_gtanh,20,'LineColor','none');
% xlabel('$\Lambda$', 'Interpreter','latex','FontSize',20)
% ylabel('$\ell_1$', 'Interpreter','latex','FontSize',20)
% cb = colorbar;
% ylabel(cb,'$S_g/S_{tanh}$','Interpreter','latex')
% set(gca,'FontSize',16,'TickLabelInterpreter','latex')
% 
% figure;
% contourf(LAM, L1GRID, S_frac_gcy,20,'LineColor','none');
% xlabel('$\Lambda$', 'Interpreter','latex','FontSize',20)
% ylabel('$\ell_1$', 'Interpreter','latex','FontSize',20)
% cb = colorbar;
% ylabel(cb,'$S_g/S_{cy}$','Interpreter','latex')
% set(gca,'FontSize',16,'TickLabelInterpreter','latex')


%% option B
l1_range = linspace(1.1,4,100);
l2_range = linspace(1.1,7,100);
l_ratio = l1_range ./ l2_range;
[L1,L2] = meshgrid(l1_range,l2_range);
valid = L2 > L1;

L_ratio = L2 ./ L1;
L_ratio(~valid) = NaN;

figure;
contourf(LAM, L1_range, S_frac_gcy, 20, 'LineColor', 'none');
colorbar;
xlabel('\Lambda', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('\ell_1', 'Interpreter', 'latex', 'FontSize', 14);
title('Influence of Graded Region Location ($\ell_1$)', 'Interpreter', 'latex');

%%
function x_sol = solve_for_x(Lambda)
    % C(Lambda)
    C = 1/Lambda^4 + 4/Lambda + 5;

    % Define the polynomial: C*x^4 - 8*x^3 - 2 = 0
    coeffs = [C, 0, 0, -8, -2]; % corresponds to x^4, x^3, ..., constant

    % Find all real roots
    roots_all = roots(coeffs);
    roots_real = roots_all(imag(roots_all) == 0); % real roots only

    % Choose the positive real root
    x_sol = roots_real(roots_real > 0);
    
    % If multiple positive roots, choose the smallest (safest)
    if length(x_sol) > 1
        x_sol = min(x_sol);
    end
end
function dydx = implicit_slope(x, y)
    u = (1 + (x^3 - 1)/y^3)^(1/3);

    dB_du = -4/u^5 - 4/u^2;

    du_dx = x^2 / (y^3 * u^2);
    df_dx = -4/x^5 - 4/x^2 + 2 * dB_du * du_dx;

    du_dy = - (x^3 - 1) / (y^4 * u^2);
    df_dy = 2 * dB_du * du_dy;

    dydx = - df_dx / df_dy;
end
 
