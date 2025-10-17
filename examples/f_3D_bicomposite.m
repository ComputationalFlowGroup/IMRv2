%% Sensitivity Map: Bi-Composite vs Mono-Composite Bubble Dynamics

% Parameter ranges
Lambda_vals = linspace(1.1, 3, 20);         % Bubble expansion ratios
G_ratio_vals = logspace(-1, 1, 20);         % G1/G0 from 0.1 to 10 (log-spaced)
l1_vals = [0.2, 0.5, 0.8];                  % Interface positions (r0/R0)

% Preallocate result array: [G_ratio x Lambda x l1]
metric_map = zeros(length(G_ratio_vals), length(Lambda_vals), length(l1_vals));

for k = 1:length(l1_vals)
    l1 = l1_vals(k);

    for i = 1:length(G_ratio_vals)
        G_ratio = G_ratio_vals(i);

        for j = 1:length(Lambda_vals)
            Lambda = Lambda_vals(j);

            % Compute sensitivity metric (RMSE)
            try
                metric = simulate_bubble_diff(Lambda, G_ratio, l1);
            catch ME
                warning('Simulation failed at Lambda=%.2f, G_ratio=%.2f, l1=%.2f', Lambda, G_ratio, l1);
                metric = NaN;
            end

            metric_map(i, j, k) = metric;
        end
    end
end

% Plot contour maps for each l1
for k = 1:length(l1_vals)
    figure;
    contourf(Lambda_vals, log10(G_ratio_vals), metric_map(:,:,k), 20);
    xlabel('\Lambda = R_{max}/R_{eq}');
    ylabel('log_{10}(G_1 / G_0)');
    title(['Sensitivity Map at \ell_1/R_0 = ', num2str(l1_vals(k))]);
    colorbar;
    colormap('turbo');
end

surf(Lambda_vals, log10(G_ratio_vals), metric_map(:,:,k));
xlabel('\Lambda'); ylabel('log_{10}(G_1/G_0)'); zlabel('RMSE');

figure;
plot(G_ratio_vals, squeeze(metric_map(:,j,k)));  % fix Lambda, l1


%% SIMULATION FUNCTIONS

%function delta = simulate_bubble_diff(Lambda, G_ratio, l1)
    % Run bi-composite simulation
    [t, R_bi] = run_graded(Lambda, G_ratio, 'bi', l1);

    % Run mono-composite simulation
    [t2, R_mono] = run_graded(Lambda, G_ratio, 'mono', l1);  % l1 is unused for mono

    % Interpolate if needed
    if length(t2) ~= length(t)
        R_mono = interp1(t2, R_mono, t, 'linear', 'extrap');
    end

    % Compute RMSE between R_bi and R_mono
    delta = sqrt(mean((R_bi - R_mono).^2));

    % Optional: normalize
    % delta = delta / max(R_mono);
%end

%function [t, R] = run_graded(Lambda, G_ratio, model_type, l1)
    addpath('src/forward_solver/');

    % Simulation constants
    kappa = 1.4;
    T8 = 298.15;
    rho8 = 1064;
    mu = 0;
    lambda1 = 0;
    lambda2 = 0;
    alphax = 0;
    Pref = 101325;
    R0 = 100e-6;
    Req = R0 / Lambda;
    tfin = 75E-6;
    tvector = linspace(0, tfin, 1000);

    collapse = 0;
    radial = 2;
    vapor = 1;
    bubtherm = 0;
    medtherm = 0;
    masstrans = 0;
    stress = 1;

    % Graded material properties
    G0 = 1E3;
    G1 = G_ratio * G0;
    l2 = 2.2;  % Outer radius
    gfun = 4;
    v_nc = 0.3;
    v_a = 2;

    % Shared options
    varin = {'progdisplay',0,...
        'radial',radial,...
        'bubtherm',bubtherm,...
        'tvector',tvector,...
        'vapor',vapor,...
        'medtherm',medtherm,...
        'masstrans',masstrans,...
        'method',23,...
        'stress',stress,...
        'collapse',collapse,...
        'mu',mu,...
        'lambda1',lambda1,...
        'lambda2',lambda2,...
        'alphax',alphax,...
        'r0',R0,...
        'req',Req,...
        'kappa',kappa,...
        't8',T8,...
        'rho8',rho8,...
        'v_nc',v_nc,...
        'v_a',v_a};

    switch model_type
        case 'bi'
            varin = [varin,...
                {'graded',1,...
                'g',G0,...
                'g1',G1,...
                'l1',l1,...
                'l2',l2,...
                'gfun',gfun}];
        case 'mono'
            varin = [varin,...
                {'graded',0,...
                'g',G0}];
        otherwise
            error('Invalid model_type. Use ''bi'' or ''mono''.');
    end

    % Run the simulation
    [t, R, ~] = f_imr_fd(varin{:}, 'Nt', 16, 'Mt', 64);
%end



%%
Lambda_vals = linspace(1.1, 10, 20);
%G_ratio_vals = logspace(-1, 1, 20);
G_ratio_vals    = [0.01, 0.1, 0.5, 1, 2, 5, 10];%, 100];  % G1/G0 values
l1_vals = linspace(0.1, 0.9, 15);
%ell1_vals   = linspace(1.001, 10, 300);  % Dimensionless l1 - but not wtr R0

metric_map = zeros(length(G_ratio_vals), length(Lambda_vals), length(l1_vals));

for k = 1:length(l1_vals)
    l1 = l1_vals(k);
    for i = 1:length(G_ratio_vals)
        G_ratio = G_ratio_vals(i);
        for j = 1:length(Lambda_vals)
            Lambda = Lambda_vals(j);
            try
                metric = simulate_bubble_diff(Lambda, G_ratio, l1);
            catch
                metric = NaN;
            end
            metric_map(i, j, k) = metric;
        end
    end
end

% save('sensitivity_map.mat','metric_map','Lambda_vals','G_ratio_vals','l1_vals');
% g_index = 10;  % or use find(G_ratio_vals == desired_value)
% imagesc(Lambda_vals, l1_vals, squeeze(metric_map(g_index,:,:))');
% xlabel('\Lambda = R_{max}/R_0'); ylabel('\ell_1 / R_0');
% title(['Sensitivity at G_1/G_0 = ', num2str(G_ratio_vals(g_index))]);
% colorbar;
l_index = 4;  % some middle value
%surf(Lambda_vals, log10(G_ratio_vals), squeeze(metric_map(:,:,l_index)));
surf(Lambda_vals, G_ratio_vals, metric_map(:,:,:));
xlabel('\Lambda'); ylabel('log_{10}(G_1/G_0)'); zlabel('Sensitivity');
title(['Sensitivity at %\ell_1/R_0 = %', num2str(l1_vals(l_index))],'Interpreter','latex');
% for k = 1:length(l1_vals)
%     figure(k); clf;
%     contourf(Lambda_vals, G_ratio_vals, metric_map(:,:,k), 20);
%     title(['\ell_1/R_0 = ', num2str(l1_vals(k))]);
%     xlabel('\Lambda'); ylabel('log_{10}(G_1/G_0)');
%     colorbar; drawnow;
%     pause(0.2);
% end
%%
[X, Y, Z] = meshgrid(G_ratio_vals, Lambda_vals, l1_vals);
scatter3(X(:), Y(:), Z(:), 36, metric_map(:), 'filled');
xlabel('G_1/G_0');
ylabel('\Lambda');
zlabel('$\ell_1/R_0$','Interpreter','latex');
colorbar;
title('Sensitivity metric');
grid on;
%%
%metric_map(i, j, k); 
% i = G_ratio index
% j = Lambda index
% k = l1 index
[X, Y, Z] = ndgrid(G_ratio_vals, Lambda_vals, l1_vals);
scatter3(X(:), Y(:), Z(:), 36, metric_map(:), 'filled');
xlabel('G_1/G_0');
ylabel('\Lambda');
zlabel('$\ell_1/R_0$','Interpreter','latex');
colorbar;
title('Sensitivity metric');
grid on;
