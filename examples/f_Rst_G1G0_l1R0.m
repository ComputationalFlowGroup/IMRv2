%% composite mat: 2 homogeneous domains
Rst = linspace(0.2,10,200); %R_max/R_0 = Lambda
G_ratio = [0.1, 0.5, 1, 2, 5, 10]; %G1/G0
l1s = [0.1,0.3,0.5,0.7,0.9]; %l1/R0

% 3D regime map
% Define parameter ranges
Lambda_vals = linspace(1.1, 3, 20);         % Bubble expansion ratios
G_ratio_vals = logspace(-1, 1, 20);         % G1/G0 from 0.1 to 10 (log-spaced)
l1_vals = [0.2, 0.5, 0.8];                  % Interface positions (r0/R0)

% call simulation function
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

% Preallocate result arrays
metric_map = zeros(length(G_ratio_vals), length(Lambda_vals), length(l1_vals));

for k = 1:length(l1_vals)
    l1 = l1_vals(k);
    
    for i = 1:length(G_ratio_vals)
        G_ratio = G_ratio_vals(i);
        
        for j = 1:length(Lambda_vals)
            Lambda = Lambda_vals(j);
            
            % Compute sensitivity metric
            metric = simulate_bubble_diff(Lambda, G_ratio, l1);
            
            % Store
            metric_map(i, j, k) = metric;
        end
    end
end
% Plot for each l1 slice
for k = 1:length(l1_vals)
    figure;
    contourf(Lambda_vals, log10(G_ratio_vals), metric_map(:,:,k), 20);
    xlabel('\Lambda = R_{max}/R_{eq}');
    ylabel('log_{10}(G_1 / G_0)');
    title(['Sensitivity Map at \ell_1/R_0 = ', num2str(l1_vals(k))]);
    colorbar;
    colormap('turbo');  % or 'parula', 'hot', etc.
end

function [t,R] = run_graded(Lambda,G_ratio,model_type,l1)
addpath('src/forward_solver/');
% options
kappa = 1.4;
T8 = 298.15;
rho8 = 1064;
mu = 0; %5E-2;
lambda1 = 0; %1e-7;
lambda2 = 0;
alphax = 0; %1e-3;
Pref = 101325;
R0 = 100e-6;
Req = R0 / Lambda; 
tfin = 75E-6;
tvector = linspace(0,tfin,1000);
collapse = 0;
radial = 2;
vapor = 1;
bubtherm = 0;
medtherm = 0;
masstrans = 0;
stress = 1;
% graded parameters
graded = 1;
v_nc = 0.3; 
v_a = 2; 
%l2 = 0;
G0 = 1E3;
G1 = G_ratio * G0;
gfun = 4;

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
    'g',G0,...
    'graded',graded,...
    'g1',G1,...
    'l1',l1,...
    'v_a',v_a,...
    'v_nc',v_nc,...
    'gfun',gfun,...
    'lambda1',lambda1,...
    'lambda2',lambda2,...
    'alphax',alphax,...
    'r0',R0,...
    'req',Req,...
    'kappa',kappa,...
    't8',T8,...
    'rho8',rho8};
switch model_type
   case 'bi'
            varin = [varin,...
                {'graded',1,...
                'g',G0,...
                'g1',G1,...
                'l1',l1,...
                'l2',l2,...
                'gfun',4}];
        case 'mono'
            varin = [varin,...
                {'graded',0,...
                'g',G0}];
end 
% generate R v t data
[t,R,~] = f_imr_fd(varin{:},'Nt',16,'Mt',64);

% figure
% hold on;
% plot(t,R,'bx')
% ylim([0 1.2])
% hold off;
end
