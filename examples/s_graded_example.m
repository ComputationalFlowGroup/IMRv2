%clc; clear; close all;
% figure 1 JMPS

Rst = linspace(0.2,10,200); %R_max/R_0

G0 = 500;
G1 = 1000;
l1 = 1.5; %nondim
l2 = 3; %nondim
a = 2.5;
n = 0.3;
b = 0.5;

%%
x1 = @(Rst) (1+(Rst.^3-1)./(l1).^3).^(1/3); %Lambda1
x2 = @(Rst) (1+(Rst.^3-1)./(l2).^3).^(1/3); %Lambda2

% ycy = @(x,Rst) (1+(G1-1)*(1+((x-x1(Rst))./(x2(Rst)-x)).^a).^((n-1)/a)).*(1./x.^5+1./x.^2);
fcy = @(x,Rst) (l2*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3) - 1)./...
      (1-l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3));
ycy = @(x,Rst) (G0+(G1-G0)*(1+( fcy(x,Rst) ).^a).^((n-1)/a)).*(1./x.^5+1./x.^2);
% ype = @(x,Rst) (G3+(G1-G3)*asinh((x-x1(Rst))./(x2(Rst)-x))).*(1./x.^5+1./x.^2);
% ympe = @(x,Rst) (G3+(G1-G3)*log((x-x1(Rst))./(x2(Rst)-x)+1)./((x-x1(Rst))./(x2(Rst)-x)).^a).*(1./x.^5+1./x.^2);
% ycr = @(x,Rst) (G3+(G1-G3)*1./(1+(x-x1(Rst))./(x2(Rst)-x)).^n).*(1./x.^5+1./x.^2);
% yscr = @(x,Rst) (G3+(G1-G3)*1./(1+(x-x1(Rst))./(x2(Rst)-x))).*(1./x.^5+1./x.^2);
% ymcr = @(x,Rst) (G3+(G1-G3)*(1./(1+(x-x1(Rst))./(x2(Rst)-x)).^n).^a).*(1./x.^5+1./x.^2);

ftanh = @(x,Rst) ((l2 + l1)/(l2 - l1)) .* (((Rst.^3 - 1)./(x.^3 - 1)).^(1/3) + ((l1+l2)/2));
ytanh = @(x,Rst) (G0+(G1-G0)*(1/2)*(1+tanh(b*(2*ftanh(x,Rst) -1)))).*(1./x.^5+1./x.^2);

ytanh_fcy = @(x,Rst) (G0+(G1-G0)*(1/2)*(1+tanh(b*(2*fcy(x,Rst) -1)))).*(1./x.^5+1./x.^2);
ycy_ftanh= @(x,Rst) (G0+(G1-G0)*(1+( ftanh(x,Rst) ).^a).^((n-1)/a)).*(1./x.^5+1./x.^2);

% linear can be checked by just using fcy and ftanh

reltol = 1e-8;
abstol = 1e-8;
S2 = zeros(1,length(Rst));
for i = 1:length(Rst)
    rst = Rst(i);
    S2(i) = 2*integral(@(x) ftanh(x,rst),x1(rst),x2(rst),...
            'RelTol',reltol,'AbsTol',abstol);
end


S0 = (G0/2)*(1./Rst.^4 + 4./Rst - (1./x1(Rst).^4 + 4./x1(Rst)));
S1 = -(G1/2)*(5 - 4./x2(Rst) - 1./x2(Rst).^4);
SG0 = -(G0/2)*(5 - 4./Rst - 1./Rst.^4);
SG1 = -(G1/2)*(5 - 4./Rst - 1./Rst.^4);
%figure(1)
figure
hold on;
% plot(Rst,S0,'m','LineWidth',3)
% plot(Rst,S1,'c','LineWidth',3)
% plot(Rst,S2,'b','LineWidth',3)
plot(Rst,SG0/G0,'r','LineWidth',3)
plot(Rst,SG1/G0,'k--','LineWidth',3)
plot(Rst,(S0+S1+S2)/G0,'-.g','LineWidth',3)
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
%saveas(gcf,'./fig_graded_stress_integral','png')


%% composite mat: 2 homogeneous domains
% goal: what value of l1 are the total stress contributions from each region equal

% Parameters
Lambda = 10;                % Rmax / R0
ell1 = linspace(1.01, 10, 500);  % Dimensionless inner boundary
Lambda1 = nthroot(1 + (Lambda^3 - 1)./ell1.^3, 3);

% Compute stress contributions
S_near = (G0/2) * (1/Lambda^4 + 4/Lambda - (1./Lambda1.^4 + 4./Lambda1));
S_far  = (G1/2) * (1./Lambda1.^4 + 4./Lambda1 - 5);
S_total = S_near + S_far;

% Fractional contributions
S_frac_near = S_near ./ S_total;
S_frac_far  = S_far  ./ S_total;

% Find where they are approximately equal
[~, idx_equal] = min(abs(S_frac_near - S_frac_far));
ell1_equal = ell1(idx_equal);

% Plot
figure;
plot(ell1, S_frac_near, 'b-', 'LineWidth', 2); hold on;
plot(ell1, S_frac_far, 'r--', 'LineWidth', 2);
xline(ell1_equal, 'k:', 'LineWidth', 2);
legend({'Near-field', 'Far-field', 'Equal contribution'}, 'Location', 'best');
xlabel('$\ell_1 = l_1 / R_0$', 'Interpreter','latex','FontSize',16);
ylabel('Fraction of total stress', 'Interpreter','latex','FontSize',16);
title(['Stress Contributions at $\Lambda = ', num2str(Lambda), '$'], 'Interpreter','latex');
grid on; box on;
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');

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

%% material contrast sensitivity map
% Parameters
Lambda_vals = linspace(1.1, 10, 100);    % Rmax/R0 (Λ)
ell1_vals   = linspace(1.001, 10, 300);  % Dimensionless l1
%G_ratios    = [0.01, 0.1, 0.5, 1, 2, 5, 10, 100];  % G1/G0 values
G_ratios = 1;

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
        
        %Lambda1 = nthroot(1 + (Lambda^3 - 1) ./ ell1.^3, 3);
        Lambda1 = ell1 / Lambda; %linear approx
        
        S_near = (G0/2) * (1/Lambda^4 + 4/Lambda - (1./Lambda1.^4 + 4./Lambda1));
        S_far  = (G1/2) * (1./Lambda1.^4 + 4./Lambda1 - 5);
        S_total = S_near + S_far;
        
        % Fractional contributions
        S_frac_near = S_near ./ S_total;
        S_frac_far = S_far ./ S_total;
        
        % Find crossover point (where fractions ~ equal)
        [~, idx_eq] = min(abs(S_frac_near - S_frac_far));
        crossover_l1(j) = ell1(idx_eq);
    end
    
    % Store and plot
    all_crossovers(k, :) = crossover_l1;
    plot(Lambda_vals, crossover_l1, 'LineWidth', 2, 'Color', colors(k,:), ...
        'DisplayName', ['G_1/G_0 = ' num2str(G_ratios(k))]);
    plot(Lambda_vals,Lambda_vals, 'r--','LineWidth',2)
end

xlabel('$\Lambda = R_{\mathrm{max}} / R_0$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\ell_1^{\mathrm{crossover}}$', 'Interpreter', 'latex', 'FontSize', 18);
title('Crossover $\ell_1$ vs $\Lambda$ for Varying Material Contrast', 'Interpreter', 'latex');
legend('Location', 'northwestoutside');
grid on; box on;
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');

p = polyfit(Lambda_vals, crossover_l1, 1);  % 1 = linear fit
slope_reg = p(1);  % slope of the fitted line
intercept_reg = p(2);
% Plot the data and the regression line
figure;
plot(Lambda_vals, crossover_l1, 'bo', 'MarkerSize', 6, 'DisplayName', 'Data'); hold on;
plot(Lambda_vals, polyval(p, Lambda_vals), 'r-', 'LineWidth', 2, 'DisplayName', 'Linear fit');
xlabel('$\Lambda$', 'Interpreter','latex', 'FontSize', 16);
ylabel('$\ell_1$ crossover', 'Interpreter','latex', 'FontSize', 16);
title('Linear regression of crossover length vs Lambda', 'Interpreter','latex');
legend('Location', 'best');
grid on; box on;
fprintf('Regression slope = %.4f\n', slope_reg);

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

%%

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
% figure 2 JMPS
% dimensional, in terms of lambda
tau = @(x) (2/3)*(1./x.^4 - x.^2);
%gtau = @(x,Rst) G3+((G1-G3)*(1+( (l2*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3) - 1)./...
%        (1-l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3)) ).^a).^((n-1)/a));
gtau = @(x,Rst) (G1+(G3-G1)*(1/2)*(1+tanh(a*(2*f(x,Rst) -1))));

gtau1 = zeros(size(Rst));
gtau2 = zeros(size(Rst));
gtau3 = zeros(size(Rst));
for i = 1:length(Rst)
    rst = Rst(i);
    gtau1(i) = G1*tau(rst);
    gtau2(i) = gtau(rst,rst)*tau(rst);
    gtau3(i) = G3*tau(rst);
end

figure(2)
hold on;
plot(Rst,gtau1,'r','LineWidth',3)
plot(Rst,gtau3,'k--','LineWidth',3)
plot(Rst,gtau2,'c--','LineWidth',3)
%plot(Rst,gtau1+gtau2+gtau3,'-.g','LineWidth',3)
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
%saveas(gcf,'./fig_graded_stress_Rst','png')

%% stress integral as a function of Lambda
% Lambda = linspace(0.1,1,200); %R(t)/R0
% Rt = Lambda * Req;
% Rst = R0/Req; %Rmax/Req
R0=1; %initial radius Rmax
Req=0.125; % equilibrium radius R0
Rst = R0/Req; % Rmax/Req (constant)
Lambda = linspace(1,Rst,200); % R(t)/Req, decreasing from Rmax to Req

G1 = 0.5;
G3 = 1;
l1 = 1.5;
l2 = 3;
a = 2.5;
n = 0.3;
x1 = @(Rst) (1+(Rst.^3-1)./(l1).^3).^(1/3); %Lambda1
x2 = @(Rst) (1+(Rst.^3-1)./(l2).^3).^(1/3); %Lambda2


% ycy = @(x,Rst) (G1+(G3-G1)*(1+( (l2*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3) - 1)./...
%       (1-l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3)) ).^a).^((n-1)/a)).*(1./x.^5+1./x.^2);
% 
% f = @(x,Rst) (l2*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3) - 1)./...
%       (1-l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3));
% ytanh = @(x,Rst) (G1+(G3-G1)*(1/2)*(1+tanh(a*(2*f(x,Rst) -1)))).*(1./x.^5+1./x.^2);

%
reltol = 1e-8;
abstol = 1e-8;
S2 = zeros(1,length(Lambda));
for i = 1:length(Lambda)
    lambda= Rst(i);
    S2(i) = 2*integral(@(x) ytanh(x,rst),x1(rst),x2(rst),...
            'RelTol',reltol,'AbsTol',abstol);
end


S1 = (G1/2)*(1./Rst.^4 + 4./Rst - (1./x1(Rst).^4 + 4./x1(Rst)));
S3 = -(G3/2)*(5 - 4./x2(Rst) - 1./x2(Rst).^4);
SG1 = -(G1/2)*(5 - 4./Rst - 1./Rst.^4);
SG3 = -(G3/2)*(5 - 4./Rst - 1./Rst.^4);
figure(1)
hold on;
% plot(Rst,S1,'m')
% plot(Rst,S2,'k')
% plot(Rst,S3,'b')
plot(Rst,SG1/G1,'r','LineWidth',3)
plot(Rst,SG3/G1,'k--','LineWidth',3)
plot(Rst,(S1+S2+S3)/G1,'-.g','LineWidth',3)
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
%saveas(gcf,'./fig_graded_stress_integral','png')
%% stress as a function of ratio l1/l2
% at an instant, how does graded material's stress in graded region change with ratio
G0 = 1000; G1 = 5000;
v_nc = 0.3; v_a = 2;
R0 = 300e-5; %Rmax
Req = R0/3; %R_0

l1_l2 = linspace(0.1,1,100);
r_coord = 1; %linspace(0.1,3,500); %fixed Eulerian grid
Rnow = R0;
r0coord = (r_coord.^3 - Rnow.^3 + Req^3);
% valid = r0coord > 0; %where r0_coord values are negative (inside bubble)
% r0_coord = (r0coord.*valid).^(1/3);
comp = zeros(size(l1_l2));
for i = 1:length(l1_l2)
    % r0_l2 = r0_coord / l2; %this still relies on a FIXED value of l2
    % fcy = (1 - r0_l2) ./ (r0_l2 - l1_l2);
    % mcy = @(fcy) (1 + fcy.^v_a).^((v_nc-1)/v_a);

    % at a fixed location - midpoint of graded region
    gamma = l1_l2(i);
    beta = 0.5*(gamma + 1);
    fcy_mid = ((1 - beta) ./ (beta - gamma));
    %mcy_mid = mcy(fcy_mid);
    mcy_mid = (1 + fcy_mid.^v_a).^((v_nc-1)/v_a);

   tau = 2/3 * ((r0_coord.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord.^2));
   taurr1 = G0*tau;
   taurr2 = (G0 + (G1 - G0)*mcy_mid).*tau;
  % taurr3 = G1*tau;

  % taurr = taurr1 + taurr2 + taurr3;

   %still have to fix l_2 to use l_1 or vice versa
   l2 = 1;
   l1 = gamma*l2;
   comp(i) = 1./ ((r0_coord./l1).* (((1/(G1-G0)).*((3/2).*(taurr2./taurr1) -1)).^(v_a/(v_nc-1)) -1).^(1/v_a) + (r0_coord/l1) - 1);
end
figure;
plot(l1_l2,comp,'k','LineWidth',3)
xlabel('$l_1 / l_2$', 'Interpreter', 'Latex', 'FontSize', 20);
%ylabel('$\tau_{rr} / \mathrm{max}(\tau_{rr})$','Interpreter','Latex','FontSize',20);
ylabel('$\tau_{rr}$','Interpreter','Latex','FontSize',20);

%% compare other non-Newtonian graded functions

Rst = linspace(0.2,10,200); %R_max/R_0

G0 = 0.5;
G1 = 1;
%should G1 and G0 be switched since G0 < G1 for positive gradient
l1 = 1.5;
l2 = 3;
a = 2;
n = 0.3;
x1 = @(Rst) (1+(Rst.^3-1)./(l1).^3).^(1/3); %Lambda1
x2 = @(Rst) (1+(Rst.^3-1)./(l2).^3).^(1/3); %Lambda2

G = G1 + (G1 - G0);

% ycy = @(x,Rst) (1+(G1-1)*(1+((x-x1(Rst))./(x2(Rst)-x)).^a).^((n-1)/a)).*(1./x.^5+1./x.^2);
ycy = @(x,Rst) (G*(1+( (l2*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3) - 1)./...
      (1-l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3)) ).^a).^((n-1)/a)).*(1./x.^5+1./x.^2);
% NEED TO CONSIDER HOW THESE CHANGE FOR EACH GRADED FUNCTION WITH l1 AND l2
ype = @(x,Rst) (G*asinh((x-x1(Rst))./(x2(Rst)-x))).*(1./x.^5+1./x.^2);
ympe = @(x,Rst) (G*log((x-x1(Rst))./(x2(Rst)-x)+1)./((x-x1(Rst))./(x2(Rst)-x)).^a).*(1./x.^5+1./x.^2);
ycr = @(x,Rst) (G*1./(1+(x-x1(Rst))./(x2(Rst)-x)).^n).*(1./x.^5+1./x.^2);
yscr = @(x,Rst) (G*1./(1+(x-x1(Rst))./(x2(Rst)-x))).*(1./x.^5+1./x.^2);
ymcr = @(x,Rst) (G*(1./(1+(x-x1(Rst))./(x2(Rst)-x)).^n).^a).*(1./x.^5+1./x.^2);
functions = {@ycy, @ympe, @ycr, @yscr, @ymcr};


reltol = 1e-8;
abstol = 1e-8;
%S2 = zeros(1,length(Rst));
S2 = zeros(length(Rst),length(functions));
for i = 1:length(Rst)
    rst = Rst(i);
    for f = 1:length(functions)
        func = functions{f};
        S2(i,f) = 2*integral(@(x) func(x,rst),x1(rst),x2(rst),...
            'RelTol',reltol,'AbsTol',abstol);
    end
end

%%
S1 = (G1/2)*(1./Rst.^4 + 4./Rst - (1./x1(Rst).^4 + 4./x1(Rst)));
S3 = -(G3/2)*(5 - 4./x2(Rst) - 1./x2(Rst).^4);
SG1 = -(G1/2)*(5 - 4./Rst - 1./Rst.^4);
SG3 = -(G3/2)*(5 - 4./Rst - 1./Rst.^4);
figure(1)
hold on;
% plot(Rst,S1,'m')
% plot(Rst,S2,'k')
% plot(Rst,S3,'b')
plot(f, ype, 'r', 'LineWidth', 2); hold on;
plot(f, ympe, 'g', 'LineWidth', 2);
plot(f, ycr, 'b', 'LineWidth', 2);
plot(f, yscr, 'm', 'LineWidth', 2);
plot(f, ymcr, 'c', 'LineWidth', 2);
plot(Rst,SG1/G1,'r','LineWidth',3)
plot(Rst,SG3/G1,'k--','LineWidth',3)
plot(Rst,(S1+S2+S3)/G1,'-.g','LineWidth',3)

% Plotting the results of different functions (S2)
figure;
hold on;
for f = 1:length(functions)
    plot(Rst, S2(:, f), 'DisplayName', func2str(functions{f}));
end
xlabel('Rst');
ylabel('Integral Result');
legend show;
title('Integral Results for Different Functions');
hold off;

% Plotting the functions themselves for a specific rst value
rst_example = Rst(1); % Example rst to plot for
x_vals = linspace(x1(rst_example), x2(rst_example), 100); % x-values for plotting

figure;
hold on;
for f = 1:length(functions)
    y_vals = arrayfun(@(x) functions{f}(x, rst_example), x_vals);
    plot(x_vals, y_vals, 'DisplayName', func2str(functions{f}));
end
xlabel('x');
ylabel('Function Value');
legend show;
title('Functions for a Specific Rst');
hold off;
