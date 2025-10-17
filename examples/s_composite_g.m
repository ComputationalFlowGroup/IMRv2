clear; clc;

% Parameters
Lambda = 2;                   % R_max / R0
G0 = 1;                       % Shear modulus near the bubble
G1 = 5;                      % Shear modulus far away
a = 2.5; n = 0.3;             % C-Y parameters
b = 0.5;                      % tanh steepness parameter

% Domain for grading width
Delta_vals = linspace(0.01, 3, 100);  % Graded region half-widths

% ---- Step 1: Find l1_crossover for this Lambda ----
ell1 = linspace(1.01, 10, 1000);
Lambda1 = nthroot(1 + (Lambda^3 - 1)./ell1.^3, 3);
S_near = (G0/2) * (1/Lambda^4 + 4/Lambda - (1./Lambda1.^4 + 4./Lambda1));
S_far  = (G1/2) * (1./Lambda1.^4 + 4./Lambda1 - 5);
S_total = S_near + S_far;
S_frac_near = S_near ./ S_total;
S_frac_far  = S_far  ./ S_total;
[~, idx_cross] = min(abs(S_frac_near - S_frac_far));
l1_cross = ell1(idx_cross);

% ---- Step 2: Loop over grading widths Δ ----
S_frac_gtanh = zeros(size(Delta_vals));
S_frac_gcy   = zeros(size(Delta_vals));
S_total_tanh = zeros(size(Delta_vals));
S_total_cy   = zeros(size(Delta_vals));

for j = 1:length(Delta_vals)
    Delta = Delta_vals(j);
    l1 = l1_cross - Delta;
    l2 = l1_cross + Delta;

    % Integration bounds
    Lambda1 = nthroot(1 + (Lambda^3 - 1)/l1^3, 3);
    Lambda2 = nthroot(1 + (Lambda^3 - 1)/l2^3, 3);

    % Stress expressions
    x = @(r) r;

    % tanh grading
    f_tanh = @(x) (l2 + l1)/(l2 - l1) * (( (Lambda^3 - 1)./(x.^3 - 1) ).^(1/3) + (l1 + l2)/2);
    m_tanh = @(x) 0.5 * (1 + tanh(b * f_tanh(x)));
    Gtanh  = @(x) G0 + (G1 - G0) .* m_tanh(x);
    ytanh  = @(x) Gtanh(x) .* (1./x.^5 + 1./x.^2);
    
    % cy grading
    f_cy = @(x) (l2*((x.^3 - 1)/(Lambda^3 - 1)).^(1/3) - 1) ./ ...
                 (1 - l1*((x.^3 - 1)/(Lambda^3 - 1)).^(1/3));
    m_cy = @(x) (1 + (f_cy(x)).^a).^((n - 1)/a);
    Gcy  = @(x) G0 + (G1 - G0).*m_cy(x);
    ycy  = @(x) Gcy(x) .* (1./x.^5 + 1./x.^2);

    % Integrate
    try
        Sg_tanh = 2 * integral(ytanh, Lambda1, Lambda2, 'AbsTol', 1e-8, 'RelTol', 1e-8);
        Sg_cy   = 2 * integral(ycy, Lambda1, Lambda2, 'AbsTol', 1e-8, 'RelTol', 1e-8);
    catch
        Sg_tanh = NaN; Sg_cy = NaN;
    end

    % Boundary stress terms
    S0 = G0/2 * (1/Lambda^4 + 4/Lambda - (1/Lambda1^4 + 4/Lambda1));
    S1 = G1/2 * (1/Lambda2^4 + 4/Lambda2 - 5);

    % Total stress
    S_total_tanh(j) = S0 + Sg_tanh + S1;
    S_total_cy(j)   = S0 + Sg_cy + S1;

    % Fractional contributions
    S_frac_gtanh(j) = Sg_tanh / S_total_tanh(j);
    S_frac_gcy(j)   = Sg_cy / S_total_cy(j);
end

% ---- Step 3: Plot results ----

figure;
plot(Delta_vals, S_frac_gtanh, 'r-', 'LineWidth', 2); hold on;
plot(Delta_vals, S_frac_gcy, 'b--', 'LineWidth', 2);
xlabel('Graded half-width $\Delta$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Graded region contribution (fraction)', 'Interpreter', 'latex', 'FontSize', 16);
legend({'Tanh grading', 'C-Y grading'}, 'Interpreter','latex');
title(['Graded stress contribution vs $\Delta$ at $\Lambda = ', num2str(Lambda), '$'], 'Interpreter','latex');
grid on; box on;

figure;
plot(Delta_vals, S_total_tanh, 'r-', 'LineWidth', 2); hold on;
plot(Delta_vals, S_total_cy, 'b--', 'LineWidth', 2);
xlabel('Graded half-width $\Delta$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Total Stress $S$', 'Interpreter', 'latex', 'FontSize', 16);
legend({'Tanh grading', 'C-Y grading'}, 'Interpreter','latex');
title('Total stress vs graded width', 'Interpreter','latex');
grid on; box on;
 %%

% Domain
Delta_vals = linspace(0.01, 3, 60);    % Graded region half-width
Lambda_vals = linspace(1.1, 10, 80);   % Stretch ratios

[DELTA, LAMBDA] = meshgrid(Delta_vals, Lambda_vals);
S_frac_gtanh = NaN(size(DELTA));
S_frac_gcy   = NaN(size(DELTA));

for i = 1:numel(DELTA)
    Lambda = LAMBDA(i);
    Delta  = DELTA(i);

    % ----- Step 1: Compute crossover l1 for current Lambda -----
    ell1_range = linspace(1.01, 15, 1000);
    Lambda1 = nthroot(1 + (Lambda^3 - 1)./ell1_range.^3, 3);
    S_near = (G0/2) * (1/Lambda^4 + 4/Lambda - (1./Lambda1.^4 + 4./Lambda1));
    S_far  = (G1/2) * (1./Lambda1.^4 + 4./Lambda1 - 5);
    S_total = S_near + S_far;
    [~, idx_cross] = min(abs(S_near - S_far));
    l1_cross = ell1_range(idx_cross);

    % ----- Step 2: Define grading bounds -----
    l1 = l1_cross - Delta;
    l2 = l1_cross + Delta;
    if l1 <= 0 || l2 <= 0
        continue;  % Skip unphysical cases
    end

    % Integration bounds
    Lambda1 = nthroot(1 + (Lambda^3 - 1)/l1^3, 3);
    Lambda2 = nthroot(1 + (Lambda^3 - 1)/l2^3, 3);
    if ~isreal(Lambda1) || ~isreal(Lambda2) || Lambda1 >= Lambda2
        continue;  % Skip invalid domains
    end

    % ---- tanh grading ----
    f_tanh = @(x) (l2 + l1)/(l2 - l1) * (( (Lambda^3 - 1)./(x.^3 - 1) ).^(1/3) + (l1 + l2)/2);
    m_tanh = @(x) 0.5 * (1 + tanh(b * f_tanh(x)));
    Gtanh  = @(x) G0 + (G1 - G0) .* m_tanh(x);
    ytanh  = @(x) Gtanh(x) .* (1./x.^5 + 1./x.^2);

    % ---- C-Y grading ----
    f_cy = @(x) (l2*((x.^3 - 1)/(Lambda^3 - 1)).^(1/3) - 1) ./ ...
                 (1 - l1*((x.^3 - 1)/(Lambda^3 - 1)).^(1/3));
    m_cy = @(x) (1 + (f_cy(x)).^a).^((n - 1)/a);
    Gcy  = @(x) G0 + (G1 - G0) .* m_cy(x);
    ycy  = @(x) Gcy(x) .* (1./x.^5 + 1./x.^2);

    try
        Sg_tanh = 2 * integral(ytanh, Lambda1, Lambda2, 'AbsTol', 1e-8, 'RelTol', 1e-8);
        Sg_cy   = 2 * integral(ycy, Lambda1, Lambda2, 'AbsTol', 1e-8, 'RelTol', 1e-8);
    catch
        Sg_tanh = NaN; Sg_cy = NaN;
    end

    % Boundary stresses
    S0 = G0/2 * (1/Lambda^4 + 4/Lambda - (1/Lambda1^4 + 4/Lambda1));
    S1 = G1/2 * (1/Lambda2^4 + 4/Lambda2 - 5);

    % Total
    S_total_tanh = S0 + Sg_tanh + S1;
    S_total_cy   = S0 + Sg_cy + S1;

    % Fractional contributions
    S_frac_gtanh(i) = Sg_tanh / S_total_tanh;
    S_frac_gcy(i)   = Sg_cy   / S_total_cy;
end
Z_gcy = S_frac_gcy;
Z_gcy(~isfinite(Z_gcy)) = NaN;
Z_tanh = S_frac_gtanh;
Z_tanh(~isfinite(Z_tanh)) = NaN;

% ---- PLOT: Graded Region Contribution ----
figure;
contourf(DELTA, LAMBDA, Z_tanh, 30, 'LineColor', 'none');
xlabel('Graded half-width $\Delta$', 'Interpreter','latex','FontSize',14);
ylabel('Stretch ratio $\Lambda$', 'Interpreter','latex','FontSize',14);
title('Tanh: Graded region stress fracstion', 'Interpreter','latex');
colorbar; colormap parula;
set(gca, 'FontSize', 12, 'TickLabelInterpreter','latex');

figure;
contourf(DELTA, LAMBDA, Z_gcy, 30, 'LineColor', 'none');
xlabel('Graded half-width $\Delta$', 'Interpreter','latex','FontSize',14);
ylabel('Stretch ratio $\Lambda$', 'Interpreter','latex','FontSize',14);
title('C-Y: Graded region stress fraction', 'Interpreter','latex');
colorbar; colormap turbo;
set(gca, 'FontSize', 12, 'TickLabelInterpreter','latex');

