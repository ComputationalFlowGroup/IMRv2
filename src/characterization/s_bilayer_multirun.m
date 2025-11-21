%% Bilayered multirun for synthetic experiments

% user-choice of params
eps = .02; % 2% variation
Lambda = 2;
Lambda_perturb = [Lambda-eps, Lambda, Lambda+eps];
l1 = 1.4e-4;
l2 = 1;
G0 = 100;
G1 = 1000;

% preallocating output file
format long;
nexp = length(Lambda_perturb);
data = zeros(nexp,3); %[Rmax, Lambdamax, tc]
%Rdata = cell(3,nexp); % 1 column for R, 2nd for Rdot, 3rd for time
Rdata = cell(1,nexp);

% call IMR, with KM, keep vapor on but turn off medtherm and masstherm
% turn off mu, lambda1, lambda2, alphax
addpath('../forward_solver/');
for i = 1:nexp
    % options
    R0 = 100e-6;
    Req = R0/Lambda_perturb(i);
    
    kappa = 1.4;
    T8 = 298.15;
    rho8 = 1064;
    mu = 0;
    lambda1 = 0;
    lambda2 = 0;
    alphax = 0;
    Pref = 101325;
    
    % compute tfin
    tfin = 1.25*R0*sqrt(rho8/Pref);
    tvector = linspace(0,tfin,1000);
    
    collapse = 0;
    radial = 1; % 1 RP, 2 KM
    vapor = 1;
    bubtherm = 0;
    %1; can't directly compare bubble pressure/temp in energy balance
    medtherm = 0;
    masstrans = 0;
    stress = 1;
    
    % graded parameters (nondim)
    el1 = l1/Req;
    el2 = l2/Req;
    
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
        'graded',1,...
        'gfun',4,...
        'g1',G1,...
        'l1',el1,...
        'l2',el2,...
        'lambda1',lambda1,...
        'lambda2',lambda2,...
        'alphax',alphax,...
        'r0',R0,...
        'req',Req,...
        'kappa',kappa,...
        't8',T8,...
        'rho8',rho8};
    
    % generate R v t data
    [t,R,U,~] = f_imr_fd(varin{:},'Nt',16,'Mt',64);
    [~,idx_minR] = min(R);
    % if want nondim tc
    %tc = t(idx_minR);
    %need to dimensionalize time
    tchar = sqrt(rho8/Pref)*R0;
    tc = t(idx_minR)*tchar;
    data(i,:) = [R0, Lambda_perturb(i), tc];

    % Rdata{1,i} = R * R0; %dim bubble radius
    % Rdata{2,i} = U * (R0/tchar); %dim velocity
    % Rdata{3,i} = t*tchar; %dimensional time
    Rdata{i} = [t(:)*tchar,R(:)*R0,U(:)*(R0/tchar)];
end
% save the output
save('data.mat','data');
save('Rdata.mat','Rdata');
% Outputs
% nX = size(data,1);   % # of experiments
% RX = data(:,1);      % All Rmax
% LX = data(:,2);      % All amplification Lmax
% T1X = data(:,3);     % All collapse time t1
% ti = Rdata{i}(:,1);     % All time
% Ri = Rdata{i}(:,2);     % spatial info of bubble wall
% Ui = Rdata{i}(:,3);     % All Rdot



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extracting stress integral from experiments

% load data file containing t, R, U
%load("Rdata.mat")

% fluid properties
rho = 1064;       % kg/m^3
c = 1500;         % m/s, sound speed in water
gamma = 0.072;    % N/m
p_inf = 101325;   % Pa

%nexp = length(Rdata);
% preallocate output
Sdata = cell(1,nexp);

% read data files
% assume you have: Rdata{i} = [t, R, Rdot]
for i = 1:nexp
    t = Rdata{i}(:,1);
    R = Rdata{i}(:,2);
    Rdot = Rdata{i}(:,3);

% bubble pressure (adiabatic ideal gas)
R0 = max(R);
kappa = 1.4;
pb = 101325*(R0./R).^(3*kappa);  % bubble internal pressure

% uniform time spacing outputted from IMR
dt = t(2)-t(1);
% compute accelerations
Rddot = compute_Rddot(dt,R);

Sdata{i} = rayleighplesset(rho,R,Rdot,Rddot,p_inf,pb,gamma);
%S_KM = kellermiksis(t,R,Rdot,Rddot,rho,c,pb,p_inf,gamma);

%plot
figure(i);
plot(t,Sdata{i},'LineWidth',1.3);
xlabel('t [s]');
ylabel('Stress integral S(t) [Pa]');
title('Extracted stress integral using Rayleigh Plesset');
end
save("Sdata_RP.mat","Sdata")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inversion problem - step 2 comparing stress integrals with effective monomaterial stress integral
% % compute effective monomaterial shear moduli Geffs for each experiment
% Geff = zeros(1,nexp);
% % add options to set optimization tolerances
% opts = optimset('Display','iter','TolX',1e-4,'TolFun',1e-4,'MaxIter',500,'MaxFunEvals',50,'PlotFcns',@optimplotfval); %'OutputFcn',@output_log
% for i = 1:nexp
%     t = Rdata{i}(:,1);
%     R = Rdata{i}(:,2);
%     Rdot = Rdata{i}(:,3);
%     S_target = Sdata{i};
%     % L2 error between target stress and mono stress
%     objfun = @(G) stress_mismatch(G,R0,R,S_target);
%     % initial guess for Geff_try = initial guess, lower and upper bounds
%     Geff_try = 500; Geff_lb = 10; Geff_ub = 100000;
%     % search for G in a constrained way with informed initial guess
%     Geff(i) = fminsearch(objfun,Geff_try,opts);
%     %Geff(i) = fminsearchbnd(objfun,Geff_try,Geff_lb,Geff_ub,opts);
%     %Geff(i) = lsqnonlin(objfun,Geff_try,Geff_lb,Geff_ub,opts);
%     %Geff(i) = lsqnonlin(objfun,Geff_try,opts);
% end

% % using Rdata, Lambda_perturb, Sdata, recover params - fitting
% obj_bilayer = @(x) bilayer_mismatch(x,R0,Lambda_perturb,Geff);
% G0_guess = 200; G0_lb = 10; G0_ub = 100000;
% G1_guess = 2000; G1_lb = 10; G1_ub = 100000;
% l1_guess = 5e-6; l1_lb = R0; l1_ub = 1e-3; %dim
% x0 = [G0_guess,G1_guess,l1_guess];
% lb = [G0_lb,G1_lb,l1_lb];
% ub = [G0_ub,G1_ub,l1_ub];
% %params = fminsearch(obj_bilayer,x0,opts);
% %params = fminsearchbnd(obj_bilayer,x0,lb,ub,opts);
% %params = lsqnonlin(obj_bilayer,x0,lb,ub,opts);
% params = lsqnonlin(obj_bilayer,x0,opts);
% 
% 
% % save
% inversion_results.Sdata = Sdata;
% inversion_results.Geff = Geff;
% inversion_results.params = params;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inversion Part 2: Sensitivity analysis and weighting functions

% Choose sensitivity mode: 'full' = compute at all time steps
%                          'adaptive' = focus around collapse or high S change
sensitivity_mode = 'adaptive';

% Preallocate sensitivity storage
param_list = {'Lambda','alpha','l1'}; % parameters to test
nparams = length(param_list);
S_sens = cell(1,nexp); % store sensitivities per experiment

for i = 1:nexp
    fprintf('Computing sensitivities for experiment %d/%d...\n',i,nexp);
    t = Rdata{i}(:,1);
    R = Rdata{i}(:,2);
    U = Rdata{i}(:,3);
    S_target = Sdata{i};
    % Characteristic collapse time
    [~,idx_collapse] = min(R);
    t_collapse = t(idx_collapse);

    % Compute derivatives using finite differences
    S_sens{i} = compute_sensitivities(R0,Lambda_perturb(i),G0,G1,l1,t,S_target,param_list,sensitivity_mode,t_collapse);
end

% Plot sensitivities with optional weighting
figure;
for i = 1:nparams
    subplot(nparams,1,i);
    for j = 1:nexp
        plot(Rdata{j}(:,1), S_sens{j}(:,i), 'LineWidth',1.2); hold on;
        legend_string = arrayfun(@(x) sprintf('Exp %d, \\Lambda=%.2f',x,Lambda_perturb(x)),1:nexp,'UniformOutput',false);
        legend(legend_strings);
    end
    xlabel('t [s]');
    ylabel(['\partial S / \partial ', param_list{i}]);
    title(['Sensitivity of S to ', param_list{i}]);
end

% Define weighting functions (collapse-centered or stress-change)
weights_collapse = cell(1,nexp);
weights_stress = cell(1,nexp);
for i = 1:nexp
    t = Rdata{i}(:,1);
    S_target = Sdata{i};
    [~,idx_collapse] = min(Rdata{i}(:,2));
    t_collapse = t(idx_collapse);
    weights_collapse{i} = generate_weighting(t,S_target,weight_type,t_collapse);
    weights_stress{i} = generate_weighting(t,S_target,'stress',t_collapse);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3: Inversion with weighted cost functions and error landscape

% Loop over experiments or perform global fit
fprintf('Starting bilayer inversion...\n');

% Initial guesses and bounds
G0_guess = 200; G0_lb = 10; G0_ub = 100000;
G1_guess = 2000; G1_lb = 10; G1_ub = 100000;
l1_guess = 5e-6; l1_lb = R0; l1_ub = 1e-3;
x0 = [G0_guess, G1_guess, l1_guess];
lb = [G0_lb, G1_lb, l1_lb];
ub = [G0_ub, G1_ub, l1_ub];

% set options for optimizer
opts_lsq = optimoptions('lsqnonlin','Display','iter','MaxIterations',500,'MaxFunctionEvaluations',1000);
opts_fmin = optimset('Display','iter','TolX',1e-4,'TolFun',1e-4,'MaxIter',500,'MaxFunEvals',50,'PlotFcns',@optimplotfval); %'OutputFcn',@output_log

% Weighted cost function
weighted_obj = @(x) weighted_bilayer_mismatch(x,R0,Lambda_perturb,Geff,weights);

% Inversion
params_opt = lsqnonlin(weighted_obj,x0,lb,ub,opts_lsq);
%params_opt = fminsearch(weighted_obj,x0,opts_fmin);
G0_opt = params_opt(1);
G1_opt = params_opt(2);
l1_opt = params_opt(3);

fprintf('Optimal parameters:\n G0 = %.2f, G1 = %.2f, l1 = %.2e\n',G0_opt,G1_opt,l1_opt);

% Generate error landscape for visualization (G0 vs G1)
G0_range = linspace(G0_lb,G0_ub,30);
G1_range = linspace(G1_lb,G1_ub,30);
[GG0,GG1] = meshgrid(G0_range,G1_range);
cost_grid = zeros(size(GG0));

for i = 1:numel(GG0)
    cost_grid(i) = norm(bilayer_Geff(GG0(i),GG1(i),l1_opt./R0,Lambda_perturb) - Geff);
end

figure;
surf(G0_range,G1_range,cost_grid);
xlabel('G0 [Pa]'); ylabel('G1 [Pa]'); zlabel('Error norm');
title('Error landscape: G0 vs G1');
shading interp; colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper functions for extracting the stress integral from experiments
function Sdata = rayleighplesset(rho,R,Rdot,Rddot,p_inf,pb,gamma)
    % store
    Sdata = rho*(R .* Rddot + 1.5*Rdot.^2) + p_inf - pb + 2*gamma./R;
end

function Sdata = kellermiksis(t,R,Rdot,Rddot,rho,c,pb,p_inf,gamma)
    % solve Q
    Q = solve_Q_KM(t, R, Rdot, Rddot, rho, c);

    % compute stress integral
    Sdata = Q - (pb - p_inf - 2*gamma./R);
end


function Rddot = compute_Rddot(dt,R)
    N = length(R);

    % Second finite difference on R(t)
        Rddot = zeros(size(R));
    for i = 2:N-1
        Rddot(i) = (R(i+1) - 2*R(i) + R(i-1))/(dt^2);
    end
    % endpoints (forward/backward difference)
    Rddot(1) = (R(3) - 2*R(2) + R(1))/(dt^2);
    Rddot(N) = (R(N) - 2*R(N-1) + R(N-2))/(dt^2);

    % % First finite difference on U(t)
    % Rddot = zeros(size(Rdot));
    % % interior points
    % for i = 2:N-1
    %     Rddot(i) = (Rdot(i+1) - Rdot(i-1))/(2*dt);
    % end
    % 
    % % endpoints (forward/backward difference)
    % Rddot(1) = (Rdot(2) - Rdot(1))/dt;
    % Rddot(N) = (Rdot(N) - Rdot(N-1))/dt;
end

function Q = solve_Q_KM(t, R, Rdot, Rddot, rho, c)
    % Compute A(t)
    A = (1 - Rdot./c).*R.*Rddot + 1.5*(1 - Rdot/(3*c)).*Rdot.^2;
    
    % Interpolation for ode solver
    R_fun = @(tt) interp1(t,R,tt,'linear','extrap');
    Rdot_fun = @(tt) interp1(t,Rdot,tt,'linear','extrap');
    A_fun = @(tt) interp1(t,A,tt,'linear','extrap');
    
    % Define RHS of linear ODE Q' = f(t,Q)
    odefun = @(tt,Q) (c/rho)/R_fun(tt) * rho*A_fun(tt) - ((c + Rdot_fun(tt))/R_fun(tt))*Q;
    
    % Solve with ode45
    [~, Qsol] = ode45(odefun, t, 0); % assume Q(0)=0
    Q = Qsol;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper functions for computing effective monolayered material
function F = stress_mismatch(G,R0,R,S_target)
Smono = -(G/2).*(5 -(R0./R).^4 - 4*R0./R);
F = norm(Smono - S_target);
end

function F = bilayer_mismatch(x,R0,Lambda_vec,Geff_data)
G0_try = x(1);
G1_try = x(2);
l1_try = x(3); 
ell1_try = l1_try / R0;
Geff_model = zeros(size(Lambda_vec));
for j = 1:length(Lambda_vec)
    Geff_model(j) = bilayer_Geff(G0_try,G1_try,ell1_try,Lambda_vec(j));
end
F = norm(Geff_model - Geff_data);
end

function Geff = bilayer_Geff(G0,G1,ell1,Lambda_vec_i)
Lambda1 = (1+ (Lambda_vec_i^3 - 1)/ell1^3)^(1/3);
A = ((1/Lambda_vec_i^4 + 4/Lambda_vec_i) - 1/Lambda1^4 - 4/Lambda1) / ...
    (1/Lambda_vec_i^4 + 4/Lambda_vec_i -5);
B = (1/Lambda1^4 + 4/Lambda1 - 5)/(1/Lambda_vec_i^4 + 4/Lambda_vec_i -5);
Geff = G0*A + G1*B;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper functions for part 2 sensitivity analysis

function S_sens = compute_sensitivities(R0,Lambda,G0,G1,l1,t,S_target,param_list,mode,t_collapse)
    nparams = length(param_list);
    S_sens = zeros(length(t),nparams);
    rel_eps = 1e-3; % relative perturbation
    for k = 1:nparams
        % Compute perturbed parameter
        x_perturb = struct('Lambda',Lambda,'alpha',G1/G0,'l1',l1);
        pname = param_list{k};
        x_val = x_perturb.(pname);
        dx = rel_eps * x_val;
        x_perturb.(pname) = x_val + dx;
        
        % Compute S with perturbed parameter
        switch pname
            case 'Lambda'
                Geff_pert = bilayer_Geff(G0,G1,l1/R0,x_perturb.Lambda);
            case 'alpha'
                Geff_pert = bilayer_Geff(G0, x_perturb.alpha*G0,l1/R0,Lambda);
            case 'l1'
                Geff_pert = bilayer_Geff(G0,G1,x_perturb.l1/R0,Lambda);
        end
        S_sens(:,k) = (Geff_pert - S_target)/dx;
    end
    
    % Apply weighting if adaptive mode
    if strcmp(mode,'adaptive')
        w = generate_weighting(t,S_target,'collapse',t_collapse);
        S_sens = S_sens .* w;
    end
end

function w = generate_weighting(t,S,weight_type,t_collapse)
    switch weight_type
        case 'collapse'
            sigma = 0.1*(max(t)-min(t));
            w = exp(-(t-t_collapse).^2/(2*sigma^2));
        case 'stress'
            dS = [0; abs(diff(S))];
            w = dS/max(dS);
        otherwise
            w = ones(size(t));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper function for part 3

function F = weighted_bilayer_mismatch(x,R0,Lambda_vec,Geff_data,weights)
    G0_try = x(1);
    G1_try = x(2);
    l1_try = x(3);
    ell1_try = l1_try/R0;
    Geff_model = zeros(size(Lambda_vec));
    for j = 1:length(Lambda_vec)
        Geff_model(j) = bilayer_Geff(G0_try,G1_try,ell1_try,Lambda_vec(j));
    end
    F = zeros(size(Geff_model));
    for j = 1:length(Geff_model)
        F(j) = weights{j}(:)' * (Geff_model(j) - Geff_data(j));
    end
end
