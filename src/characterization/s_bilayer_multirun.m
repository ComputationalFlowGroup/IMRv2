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
pb = p_inf*(R0./R).^(3*kappa);  % bubble internal pressure

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inversion Part 2: Sensitivity analysis and weighting functions

% Preallocate sensitivity storage
param_list = {'Lambda','alpha','l1'}; % parameters to test
nparams = length(param_list);

S_sens = cell(1,nexp); % store sensitivities per experiment
for i = 2
    fprintf('Computing sensitivities for experiment %d/%d...\n',i,nexp);
    t = Rdata{i}(:,1);
    R = Rdata{i}(:,2);
    U = Rdata{i}(:,3);
    S_target = Sdata{i};
    % Characteristic collapse time
    % [~,idx_collapse] = min(R);
    % t_collapse = t(idx_collapse);
    t_collapse = data(i,3);

    % Compute derivatives using finite differences
    S_sens{i} = compute_sensitivities(R0,Lambda_perturb(i),G0,G1,l1,t,S_target,param_list,rho,p_inf,gamma,kappa);
end
%%
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

weights = cell(nparams,1);  % one cell per parameter
for i_param = 1:nparams
    weights{i_param} = cell(1,nexp); % one weighting vector per experiment
    
    for i_exp = 1:nexp
        t = Rdata{i_exp}(:,1);
        S_target = Sdata{i_exp};
        S_sens_param = S_sens{i_exp}(:,i_param);  % sensitivity for this param
        tc = data(i_exp,3); %tc = data(:,i_exp);
        % Generate weight automatically based on sensitivity
        weights{i_param}{i_exp} = generate_weighting(t, tc, S_target, S_sens_param, 'auto');
    end
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
weighted_obj = @(x) weighted_bilayer_mismatch(x,R0,Lambda_perturb,Sdata,S_sens,param_list,weights);

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
    % for i = 2:N-1
    %     Rddot(i) = (R(i+1) - 2*R(i) + R(i-1))/(dt^2);
    % end
    Rddot(2:end-1) = (R(3:end) - 2*R(2:end-1) + R(1:end-2)) /(dt^2);
    % endpoints (forward/backward difference)
    Rddot(1) = (R(3) - 2*R(2) + R(1))/(dt^2);
    Rddot(end) = (R(end) - 2*R(end-1) + R(end-2))/(dt^2);

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

function Rddot = tikhonov_Rddot(dt, R, lambda)

N = length(R);

% Construct second derivative finite difference matrix
e = ones(N,1);
D2 = spdiags([e -2*e e], -1:1, N, N) / dt^2;

% Tikhonov regularization matrix (identity)
L = speye(N);

% Solve the Tikhonov problem:
% (D2' D2 + lambda L' L) u = D2' R
A = (D2' * D2 + lambda * (L' * L));
b = D2' * R;

Rddot = A \ b;   % Solve linear system

end

function Rddot = sgolay_Rddot(dt, R)

% Choose polynomial order and frame size
order = 5;      % Typical: 3–5
frame = 11;     % Must be odd (7, 9, 11, 13)

% Get SG differentiation coefficients
[b,g] = sgolay(order, frame);

N = length(R);
half = (frame-1)/2;

Rddot = zeros(size(R));

% Apply SG filter for second derivative
for n = half+1 : N-half
    % Convolve with 2nd derivative coefficient row
    Rddot(n) = dot(g(:,3), R(n-half:n+half)) / dt^2;
end

% Pad endpoints by copying nearest valid value
Rddot(1:half) = Rddot(half+1);
Rddot(end-half+1:end) = Rddot(end-half);

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
function S_sens = compute_sensitivities(R0, Lambda, G0, G1, l1, t, S_target, param_list, rho, p_inf, gamma, kappa)
    nparams = length(param_list);
    nt = length(t);
    S_sens = zeros(nt,nparams);

    % Relative perturbation
    rel_eps = 1e-3; % can adjust or make adaptive

    alpha = G1/G0;
    dt = t(2) - t(1);
    
    parfor i = 1:nparams
        param = param_list{i};
        % Set perturbation
        % switch param
        %     case 'Lambda', delta = Lambda*rel_eps;
        %     case 'alpha',  delta = (G1/G0)*rel_eps;
        %     case 'l1',     delta = l1*rel_eps;
        %     otherwise, error('Unknown parameter');
        % end
        % 
        % % Perturb parameter and run forward solver
        % params_plus = struct('Lambda',Lambda,'G0',G0,'G1',G1,'l1',l1);
        % if strcmp(param,'alpha')
        %     params_plus.G1 = G0*(G1/G0 + delta); % perturb alpha = G1/G0
        % else
        %     params_plus.(param) = params_plus.(param) + delta;
        % end
        switch param
            case 'Lambda'
                fprintf(' Perturbing Lambda...\n')
                delta = Lambda*rel_eps;
                params_plus = struct('Lambda',Lambda+delta,'G0',G0,'G1',G1,'l1',l1);
            case 'alpha' 
                fprintf(' Perturbing alpha...\n')
                delta = alpha*rel_eps;
                alpha_p = alpha + delta;
                G1_p = G0*alpha_p;
                params_plus = struct('Lambda',Lambda,'G0',G0,'G1',G1_p,'l1',l1);
            case 'l1' 
                fprintf(' Perturbing l1...\n')
                delta = l1*rel_eps;
                params_plus = struct('Lambda',Lambda,'G0',G0,'G1',G1,'l1',l1+delta);
            otherwise, error('Unknown parameter');
        end

        % Forward solver
        [t_sim, R_sim, U_sim] = f_imr_fd_wrapper(params_plus,R0);

        % ensure time alignment
        if length(t_sim) ~= length(t)
            R_sim = interp1(t_sim,R_sim,t,'pchip');
            U_sim = interp1(t_sim,U_sim,t,'pchip');
        end
        % compute Rddot
        Rddot_sim = compute_Rddot(dt,R_sim);
        % compute internal gas bubble pressure
        %Rmax = max(R_sim);
        pbp = p_inf*(R0 ./ R_sim).^(3*kappa);
        S_sim = rayleighplesset(rho, R_sim, U_sim, Rddot_sim, p_inf, pbp, gamma);

        % Sensitivity by finite difference
        S_sens(:,i) = (S_sim - S_target)/delta;
    end
end

%function w = generate_weighting(t, S, weight_type, tc)
function w = generate_weighting(t, tc, S_target, S_sens_param, weight_type);
%GENERATE_WEIGHTING Generate time-dependent weighting for inversion
% Inputs:
%   t            - time vector
%   tc           - collapse time (scalar)
%   S_target     - stress integral over time per experiment
%   S_sens_param - sensitivity vector for given param
%   weight_type  - 'collapse', 'stress', or 'auto'
%
% Output:
%   w           - weighting vector same size as t

col_window = 0.05; %fraction of tc for collapse width
% Gaussian centered at collapse
sigma = col_window*tc; % width = col_window % of collapse time (adjustable)
w_col = exp(-0.5*((t - tc)/sigma).^2);
% Weight proportional to stress evolution magnitude
deltaS = diff([0;S_target]); %[0; diff(S)];   % simple forward difference to compute delta S
%w = abs(dS) + 1e-6;  % avoid zero weight
w_evo = 1 + (deltaS.^2) / max(deltaS.^2);

switch lower(weight_type)
    case 'collapse'
        w = w_col;
        
    case 'stress'
        w = w_evo;
           
     case 'auto'
        
        % Fraction of total sensitivity magnitude near collapse (±20% tc)
        dt = t(2)-t(1);
        window = round(sigma/dt); 
        idx_window = max(find(t<=tc,1)-window,1):min(find(t>=tc,1)+window,length(t)); %(t >= tc - 0.2*tc) & (t <= tc + 0.2*tc));
        col_influence = sum(abs(S_sens_param(idx_window))) / sum(abs(S_sens_param)); 
        %close to 1, collapse dominated; close to 0, not-collapse dominated

        % relative importance of sensitivity at non-first-collapse time regions
        % num   = sum(abs(S_sens_param) .* w_evo);
        % denom = sum(abs(S_sens_param) .* (w_evo + w_col));
        % beta = num / denom;
        % 
        % w = beta*w_evo + (1-beta)*w_col;
        
        % hybrid
        w = col_influence * w_col + (1-col_influence)*w_evo; %
    otherwise
        error('Unknown weight type. Choose collapse, stress, or auto.');
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper function for part 3

function F = weighted_bilayer_mismatch(x, R0, Lambda_vec, Sdata, param_name, param_list, weights)
    nexp = length(Sdata);
    nparams = length(param_list);
    F = [];
    
    for j = 1:nexp
        S_target = Sdata{j};

	% model prediction for this experiment
        Geff_model = bilayer_Geff(x(1), x(2), x(3)/R0, Lambda_vec(j));
        
        % Compute weight vector
       % t = Rdata{j}(:,1);
       % switch weight_type
       %     case 'collapse'
       %         [~, idx_minR] = min(S_target); % approximate collapse time
       %         tc = t(idx_minR);
       %         w = collapse_weight(t, tc, 0.2);
       %     case 'stress_change'
       %         w = stress_change_weight(S_target, 1e-3);
       % end
        for k = 1:nparams
	w = weights{k}{j} %weighting factor for this param and experiment
        % Weighted residual (time-resolved)
        F_j = w .* (Geff_model - S_target); 
        F = [F; F_j]; % stack all experiments and time points
	end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper function for sensitivities, which requires running IMR multiple times

function [t_dim, R_dim, U_dim] = f_imr_fd_wrapper(params, R0)
% ---------------------------------------------------------------
% A clean wrapper for running f_imr_fd for inversion & sensitivity
%
% INPUT:
%   params.Lambda   – amplification factor
%   params.G0       – shear modulus layer 0
%   params.G1       – shear modulus layer 1
%   params.l1       – physical thickness (m)
%   R0              – maximum bubble radius (m)
%
% OUTPUT (DIMENSIONAL):
%   t_dim  – time [s]
%   R_dim  – radius [m]
%   U_dim  – Rdot [m/s]
%
% ---------------------------------------------------------------

    % ------------------------------
    % 1. Unpack parameters
    % ------------------------------
    Lambda = params.Lambda;
    G0     = params.G0;
    G1     = params.G1;
    l1     = params.l1;
    %radial_eq = params.radial;

    % ------------------------------
    % 2. Physical constants
    % ------------------------------
    rho8   = 1064;      % water density
    Pref   = 101325;    % ambient pressure
    kappa  = 1.4;
    T8     = 298.15;
    mu     = 0;
    lambda1 = 0;
    lambda2 = 0;
    alphax  = 0;
    vapor   = 1;
    medtherm = 0;
    masstrans = 0;
    bubtherm = 0;
    stress = 1;
    radial = 1; %radial_eq;

    % ------------------------------
    % 3. Compute equilibrium radius
    % ------------------------------
    Req = R0 / Lambda;

    % ------------------------------
    % 4. Characteristic collapse time
    % ------------------------------
    tchar = sqrt(rho8 / Pref) * R0;
    % total simulation time (same as your code)
    tfin = 1.25 * R0 * sqrt(rho8 / Pref);
    tvector = linspace(0, tfin / tchar, 1000); % f_imr_fd expects nondimensional time

    % ------------------------------
    % 5. Grade the thickness (nondimensional)
    % ------------------------------
    el1 = l1 / Req;
    el2 = 1;    % you had l2 = 1, so keep for now

    % ------------------------------
    % 6. Construct varin input list
    % ------------------------------
    varin = {'progdisplay',0,...
             'radial',radial,...
             'bubtherm',bubtherm,...
             'tvector',tvector,...
             'vapor',vapor,...
             'medtherm',medtherm,...
             'masstrans',masstrans,...
             'method',23,...
             'stress',stress,...
             'collapse',0,...
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

    % ------------------------------
    % 7. Call IMR forward solver
    % ------------------------------
    [t_ndim, R_ndim, U_ndim, ~] = f_imr_fd(varin{:}, 'Nt', 16, 'Mt', 64);

    % ------------------------------
    % 8. Convert to dimensional output
    % ------------------------------
    t_dim = t_ndim * tchar;
    R_dim = R_ndim * R0;
    U_dim = U_ndim * (R0 / tchar);
end
