% file s_NSLIC_inverse_char.m
% brief thile script contains the code to load in non-=spherical laser
% induced cavitation radius versus time and perturbation versus time data.
% The code then uses a hybrid optimization method, BO for course graining
% the parameter space and then lsqcurvefit for fine grained refining
% comparing the experimental data to the IMR forward solver.


clear all
close all
clc

%% general optimizer setup/ options + global parameters

% parameter bounds
lb = log10([1e-3, 1e-4]);
ub = log10([1e1, 1e-1]);

% optimization options
lsqopts = optimoptions('lsqcurvefit', ...
    'Display', 'iter-detailed', ...
    'OptimalityTolerance', 1e-12, ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, ...
    'MaxFunctionEvaluations', 75, ...
    'MaxIterations', 25, ...
    'UseParallel', false, ...
    'FiniteDifferenceType', 'central', ...
    'Algorithm', 'interior-point'); 
bayesoptVars = [
    optimizableVariable('Alpha', [lb(1) ub(1)])
    optimizableVariable('Mu', [lb(2), ub(2)])
    ];

%% Load In data

addpath Synthetic_data\ ../forward_solver/ ../common/

load("Synthetic_data\synthetic_data_NSLIC.mat")

datatype = 'synthetic';

perturbed = 1; % turns on characterization for surface perturbations, 1 is on


%% Loop over all datasets, extract material properties

% preallocate space for goodness of fit and extracted material properties
xsole = zeros(length(kindata), length(lb));
R2 = zeros(length(kindata),1);

tic
for s = 1:length(kindata)
    %extract experiment data
    exp = kindata{s};

    %assign data from experiment to variables
    texp = exp.time;
    R = exp.R;
    Req = exp.Req;
    epnm = exp.epnm;
    Gqs = exp.G;
    ST = exp.ST;
    fps = exp.fps;
    rho = exp.rho;
    n = exp.n; m = exp.m;

    % characteristic values
    [Rmax, idx_max] = max(R);
    Lc = Rmax;
    rhoc = rho;
    pc = 101325;
    tc = Lc*sqrt(rhoc/pc);

    % nondim values
    t_nondim = texp./tc;
    R_nondim = R./Lc;

    % finite difference stencils for computing perturbation initial
    % velcoty
    fdstenc_for = [-49/20 6 -15/2 20/3 -15/4 6/5 -1/6]; % forward finite difference
    fdstenc_cent = [-1/60 3/20 -3/4 0 3/4 -3/20 1/60]; % central finite difference

    % smooth perturbation data from experiment
    if isequal(datatype, 'exp')
        for i = 1:length(n)
            epnm(:,i) = sgolayfilt(fillmissing(epnm(:,i), 'spline'), 3, 15);
        end
    end

    % Center the data around t=0 at Rmax if needed, compute initial
    % perturbation velocities, and chinm values
    chinm = zeros(length(n), 1);
    epnmd0 = chinm;
    if idx_max == 1 || isequal(datatype, 'synthetic')
        for i = 1:length(n)
            chinm(i) = compute_chi(n(i), m(i));
            epnmd0(i) = dot(epnm(1:7, i),fdstenc_for)/mean(diff(t_nondim));
        end
    else
        texp(1:idx_max) = [];
        R(1:idx_max) = [];
        for i = 1:length(n)
            chinm(i) = compute_chi(n(i), m(i));
            epnmd0(i) = dot(epnm(idx_max-3:idx_max+3, i),fdstenc_cent)/mean(diff(t_nondim));
        end
        epnm(1:idx_max, :) = [];
        texp = texp-texp(1);
    end

    % define perturbation initial condition
    epnm0 = epnm(idx_max,:);

    % Find idx at collapse, only fit perturbations to collapse
    if isequal(datatype, 'synthetic')
        [~, idx_col] = min(R(t_nondim < 1));
    end
    epnm(idx_col+1:end, :) = [];

    % compute norms of the experimental data
    sR = norm(R_nondim);
    sEP = vecnorm(epnm);

    % equally weight the radial data with the sum of the perturbation
    % data
    wR = length(n);
    wEP = ones(1, length(n));

    % scaling factors
    aR = sqrt(wR)/sR;
    aEP = sqrt(wEP)./sEP;

    % create y_data vector which contains the radial data then all of
    % the perturbation data, both are non-dimensional
    y_data = [aR.*R_nondim; reshape(epnm.*aEP, [], 1)];

    if isequal(datatype, 'exp')
        texp(R==0) = [];
        epnm(R==0, :) = [];
        R(R==0) = [];
    end

    % create x_data struct needed for the optimization
    x_data = struct('time', texp, 'n', n, 'Req', Req, 'chi', chinm, 'ST', ST, ...
        'rho', rho, 'Gqs', Gqs, 'epnmd0', epnmd0, 'epnm0', epnm0, 'R0', Rmax, ...
        'idx_col', idx_col, 'perturbed', perturbed, 'aR', aR, 'aEP', aEP);

    % ------------- Perform Optimization --------------- %

    % define objective function that outputs the y_data_sim
    objfun = @(params, x_data) f_run_fd_IMR(10.^params, x_data);

    % create scalar function that computes sum or squared error
    objfun_scalar = @(params, x_data, y_data) sum((y_data - objfun(params, x_data)).^2);

    % create loss function that can be read by bayesopt
    lossfun = @(T) objfun_scalar([T.Alpha, T.Mu], x_data, y_data);

    % Start with bayesopt to quickly narrow the parameter space to a few
    % candidates
    results = bayesopt(lossfun, bayesoptVars, ...
        'MaxObjectiveEvaluations', 250, ...
        'UseParallel', true, ...
        'IsObjectiveDeterministic', true, ...
        'AcquisitionFunctionName', 'lower-confidence-bound', ...
        'NumSeedPoints', 50, ...
        'ExplorationRatio', 0.75, ...
        'Verbose', 1, ...
        'PlotFcn', {@plotObjectiveModel, @plotMinObjective});

    % results is the bayesopt output
    OT = results.ObjectiveTrace;      % numeric vector
    XT = results.XTrace;              % table of tried points (log-scale if you used log variables)

    % number of restarts for lsqcurvefit, extract top n_starts performers
    % from bayes opt to feed into lsqcurvefit
    n_starts = 5;
    [~,I] = sort(OT, 'ascend');
    k = min(n_starts, numel(I));           
    topk_table = XT(I(1:k),:);
    topk_mat = table2array(topk_table);
    starts = topk_mat;

    % create start points for the optimizer
    custom_starts = CustomStartPointSet(starts);

        % set up lsqcurvefit multirun optimization problem
    problem = createOptimProblem('lsqcurvefit', ...
        'x0', starts(1,:), ...
        'objective', objfun, ...
        'xdata', x_data, ...
        'ydata', y_data, ...
        'lb', lb, ...
        'ub', ub, ...
        'options', lsqopts);
    ms = MultiStart('UseParallel', true, 'StartPointsToRun','all');

    % run the optimizer, extract the best result
    [xsol, fval, exitflag, output, solutions] = run(ms, problem, custom_starts);

    % evaluate goodness of fit
    sol = objfun(xsol, x_data);
    R2(s) = 1-sum((sol-y_data).^2)/sum((y_data-mean(y_data)).^2);
    xsole(s, :) = 10.^xsol;
end
toc

%%
addpath ../../../cmap/
cmap = viridis(4); 

alph_e = xsole(:,1);
mu_e = xsole(:,2);


figure
countnew = length(alph_e);
% ===== Left subplot (alph) =====
subplot(1,2,1)
edgesalph = logspace(log10(min(alph_e)), 1.21*log10(max(alph_e)), ceil(sqrt(countnew))+1);
histogram(alph_e, edgesalph, 'FaceColor', cmap(1,:), 'EdgeColor', 'k')
hold on
xline(mean(alph_e), '--k', 'LineWidth', 2)
xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', 18)
ylabel('Count', 'FontSize', 18, 'Interpreter', 'latex')
set(gca, ...
    'TickLabelInterpreter', 'latex', ...
    'FontSize', 13, ...
    'XScale', 'log', ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off')

% Force evenly spaced ticks (log spacing but numeric labels)
xticksalph = logspace(log10(0.75*min(edgesalph)), log10(1.25*max(edgesalph)), 4);

grid on
xlim([10^(floor(log10(min(edgesalph)))) 10^(ceil(log10(max(edgesalph))))])
ylim([0 10])

% ===== Right subplot (mu) =====
subplot(1,2,2)
edgesMu = logspace(1.1*log10(min(mu_e.*1e3)), 1.1*log10(max(mu_e.*1e3)), ceil(sqrt(countnew))+1);
histogram(mu_e.*1e3, edgesMu, 'FaceColor', cmap(3,:), 'EdgeColor', 'k')
hold on
xline(mean(mu_e).*1e3, '--k', 'LineWidth', 2)
xlabel('$\mu$ [mPa$\cdot$s]', 'Interpreter', 'latex', 'FontSize', 18)
set(gca, ...
    'TickLabelInterpreter', 'latex', ...
    'FontSize', 13, ...
    'XScale', 'log', ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off')
grid on
xlim([10^(floor(log10(min(edgesMu)))) 10^(ceil(log10(max(edgesMu))))])
ylim([0 10])


%%
figure
plot(y_data, 'o')
hold on
plot(sol, '-')