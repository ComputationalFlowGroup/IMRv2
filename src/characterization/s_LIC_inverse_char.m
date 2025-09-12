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
lb = log10([1e-3, 1e-3]);
ub = log10([1e1, 1e-1]);

% optimization options
lsqopts = optimoptions('lsqcurvefit', ...
    'Display', 'iter-detailed', ...
    'OptimalityTolerance', 1e-12, ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, ...
    'MaxFunctionEvaluations', 25, ...
    'MaxIterations', 10, ...
    'UseParallel', false, ...
    'FiniteDifferenceType', 'central', ...
    'Algorithm', 'interior-point'); 
bayesoptVars = [
    optimizableVariable('Alpha', [lb(1) ub(1)])
    optimizableVariable('Mu', [lb(2), ub(2)])
    ];

%% Load In data

addpath Synthetic_data\ ../forward_solver/ ../common/

%load("Synthetic_data\synthetic_data_NSLIC.mat")
% load("Synthetic_data\synthetic_data_NSLIC_extracted_params.mat")
load("../../../Experimental_data/Processed_data/LIC/ns_Jin_polyacr.mat")
% datatype = 'synthetic';
datatype = 'exp';


G = 2.77e3; sigma = 0.056; rho = 1048;
%G = 2e3;
perturbed = 1; % turns on characterization for surface perturbations, 1 is on


%% Loop over all datasets, extract material properties

% preallocate space for goodness of fit and extracted material properties
%xsole = zeros(length(kindata), length(lb));
R2 = zeros(length(kindata),1);

tic
for s = 1:length(kindata)
    clear y_data x_data
    %extract experiment data
    exp = kindata{s};

    %assign data from experiment to variables
    texp = exp.time;
    R = exp.R;
    Req = exp.Req;
    epnm = exp.epnm;
    if isequal(datatype, 'synthetic')
        Gqs = exp.G;
        ST = exp.ST;
        rho = exp.rho;
        fps = exp.fps;
    else
        Gqs = G; ST = sigma;
        epnm = epnm';
        fps = 2e6;
        texp = 0:length(R)-1;
        texp = texp./fps;
    end
    
    n = exp.n; m = exp.m;

    % characteristic values
    [Rmax, idx_max] = max(R);
    Lc = Rmax;
    rhoc = rho;
    pc = 101325;
    tc = Lc*sqrt(rhoc/pc);

    % nondim values
    t_nondim = texp./tc;
    if size(R,1) == 1
        R = R';
    end
    

    % finite difference stencils for computing perturbation initial
    % velcoty
    fdstenc_for = [-49/20 6 -15/2 20/3 -15/4 6/5 -1/6]; % forward finite difference
    fdstenc_cent = [-1/60 3/20 -3/4 0 3/4 -3/20 1/60]; % central finite difference
    fdstencfor4 = [-25/12 4 -3 4/3 -1/4];

    % smooth perturbation data from experiment
    if isequal(datatype, 'exp')
        for i = 1:length(n)
            epnm(:,i) = sgolayfilt((fillmissing(epnm(:,i), 'spline')), 2, 7);
        end
    end

    % Center the data around t=0 at Rmax if needed, compute initial
    % perturbation velocities, and chinm values
    chinm = zeros(length(n), 1);
    epnmd0 = chinm;
    if idx_max == 1 || isequal(datatype, 'synthetic')
        idx_col1 = find(((texp-texp(idx_max))./tc < 1.25) == 1, 1, 'last');
        for i = 1:size(epnm,2)
            window = 1:1+round((idx_col1-idx_max)/8);
            twin = t_nondim(window);
            ywin = epnm(window,i);
            p = polyfit(twin, ywin, 1);   % quadratic
            epnmd0(i) = polyval(polyder(p), twin(1));
            epnm0(i) = epnm(1,i);
        end
    else
        idx_col1 = find(((texp-texp(idx_max))./tc < 1.25) == 1, 1, 'last');
        texp(1:idx_max-1) = [];
        R(1:idx_max-1) = [];

        for i = 1:15
            chinm(i) = compute_chi(n(i), m(i));
            window = idx_max:idx_max+round((idx_col1-idx_max)/5);
            twin = t_nondim(window);
            ywin = epnm(window,i);
            p = polyfit(twin, ywin, 1);   % quadratic
            epnmd0(i) = polyval(polyder(p), twin(1));
        end
        for i = 1:length(epnmd0)
            epnm0(i) = mean(epnm(idx_max:idx_max+5,i));
        end
        epnm(1:idx_max-1, :) = [];
        texp = texp-texp(1);
    end
    R_nondim = R./Lc;
    % Find idx at collapse, only fit perturbations to collapse
    % if isequal(datatype, 'synthetic')
        [~, idx_col] = min(R(t_nondim < 1));
    % end
    epnm(idx_col+1:end, :) = [];

    % compute norms of the experimental data
    sR = norm(R_nondim);
    sEP = vecnorm(epnm);

    % equally weight the radial data with the sum of the perturbation
    % data
    wR = length(n);
    wEP = ones(1, length(n));

    % scaling factors
    aR = wR/sR;
    aEP = 1./sEP;
    
    count = 0;
    idx = [];
    for i = 1:length(n)
        if ((sEP(i)/size(epnm,1) >= 1e-3) && i < 15)
            count = count + 1;
            epnm(:,count) = epnm(:,i);
        else
            idx = [idx, i];
        end
    end
    if count == 0
        continue
    end
    epnm(:,count+1:end) = [];
    n(idx) = [];
    m(idx) = [];
    aEP(idx) = [];
    epnm0(idx) = [];
    epnmd0(idx) = [];
    chinm(idx) = [];
    sEP(idx) = [];
    if isequal(datatype, 'exp')
        texp(R==0) = [];
        epnm(R(1:size(epnm,1))==0, :) = [];
        R_nondim(R==0) = [];
    end

    % create y_data vector which contains the radial data then all of
    % the perturbation data, both are non-dimensional
    y_data = [aR.*R_nondim; reshape(epnm.*aEP, [], 1)];

    % create x_data struct needed for the optimization
    x_data = struct('time', texp, 'n', n, 'Req', Req, 'chi', chinm, 'ST', ST, ...
        'rho', rho, 'Gqs', Gqs, 'epnmd0', epnmd0, 'epnm0', epnm0, 'R0', Rmax, ...
        'idx_col', idx_col, 'perturbed', perturbed, 'aR', aR, 'aEP', aEP);
    figure
    plot(texp(1:size(epnm,1))./tc, epnm, 'o--')
    hold on
    for i = 1: length(epnmd0)
        % yline(epnm0(i))
        plot(texp(1:size(epnm,1))./tc,epnmd0(i).*texp(1:size(epnm,1))./tc+epnm0(i), '-')
    end


     % ------------- Perform Optimization --------------- %

     %define objective function that outputs the y_data_sim
     objfun = @(params, x_data) f_run_fd_IMR(10.^params, x_data);

    % create scalar function that computes sum or squared error
    objfun_scalar = @(params, x_data, y_data) sqrt(sum((y_data - objfun(params, x_data)).^2))/norm(y_data);

    % create loss function that can be read by bayesopt
    lossfun = @(T) objfun_scalar([T.Alpha, T.Mu], x_data, y_data);

    % Start with bayesopt to quickly narrow the parameter space to a few
    % candidates
    results = bayesopt(lossfun, bayesoptVars, ...
        'MaxObjectiveEvaluations', 50, ...
        'UseParallel', true, ...
        'IsObjectiveDeterministic', true, ...
        'AcquisitionFunctionName', 'lower-confidence-bound', ...
        'NumSeedPoints', 5, ...
        'ExplorationRatio', 0.75, ...
        'Verbose', 1, ...
        'PlotFcn', {@plotObjectiveModel, @plotMinObjective});

    % results is the bayesopt output
    OT = results.ObjectiveTrace;      % numeric vector
    XT = results.XTrace;              % table of tried points (log-scale if you used log variables)

    % number of restarts for lsqcurvefit, extract top n_starts performers
    % from bayes opt to feed into lsqcurvefit
    n_starts = 3;
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
    figure 
    plot(sol, '-')
    hold on
    plot(y_data,'o')
end
toc

%% Histograms for extracted material properties from multiple datasets
 clear all
 clc
% 
load("Synthetic_data\synthetic_data_NSLIC_extracted_params.mat")

% load("../../../Experimental_data/Processed_data/LIC/ns_Jin_polyacr_optimized_newbds_fixedICs_new.mat")

addpath ../../../cmap/

cmap = viridis(4); 

alph_e = xsole(R2>= 0.9 ,1);
mu_e = xsole(R2>= 0.9 ,2);

figure
countnew = length(alph_e);
% ===== Left subplot (alph) =====
subplot(1,2,1)
edgesalph = logspace(0.9*log10(min(alph_e)), 1.1*log10(max(alph_e)), ceil(sqrt(countnew))+1);
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
edgesMu = logspace(0.9*log10(min(mu_e.*1e3)), 1.1*log10(max(mu_e.*1e3)), ceil(sqrt(countnew))+1);
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


%% Plot results
clear all
clc

load("../../../Experimental_data/Processed_data/LIC/ns_Jin_polyacr_optimized.mat")

% preallocate space for goodness of fit and extracted material properties
%xsole = zeros(length(kindata), length(lb));
%R2 = zeros(length(kindata),1);

tic
for s = 1:length(kindata)
    if R2(s) < 0.8
        continue
    end
 clear y_data x_data
    %extract experiment data
    exp = kindata{s};

    %assign data from experiment to variables
    texp = exp.time;
    R = exp.R;
    Req = exp.Req;
    epnm = exp.epnm;
    if isequal(datatype, 'synthetic')
        Gqs = exp.G;
        ST = exp.ST;
        rho = exp.rho;
        fps = exp.fps;
    else
        Gqs = G; ST = sigma;
        epnm = epnm';
        fps = 2e6;
        texp = 0:length(R)-1;
        texp = texp./fps;
    end
    
    n = exp.n; m = exp.m;

    % characteristic values
    [Rmax, idx_max] = max(R);
    Lc = Rmax;
    rhoc = rho;
    pc = 101325;
    tc = Lc*sqrt(rhoc/pc);

    % nondim values
    t_nondim = texp./tc;
    if size(R,1) == 1
        R = R';
    end
    

    % finite difference stencils for computing perturbation initial
    % velcoty
    fdstenc_for = [-49/20 6 -15/2 20/3 -15/4 6/5 -1/6]; % forward finite difference
    fdstenc_cent = [-1/60 3/20 -3/4 0 3/4 -3/20 1/60]; % central finite difference
    fdstencfor4 = [-25/12 4 -3 4/3 -1/4];

    % smooth perturbation data from experiment
    if isequal(datatype, 'exp')
        for i = 1:length(n)
            epnm(:,i) = sgolayfilt((fillmissing(epnm(:,i), 'spline')), 2, 7);
        end
    end

    % Center the data around t=0 at Rmax if needed, compute initial
    % perturbation velocities, and chinm values
    chinm = zeros(length(n), 1);
    epnmd0 = chinm;
    if idx_max == 1 || isequal(datatype, 'synthetic')
        idx_col1 = find(((texp-texp(idx_max))./tc < 1.25) == 1, 1, 'last');
        for i = 1:size(epnm,2)
            window = 1:1+round((idx_col1-idx_max)/8);
            twin = t_nondim(window);
            ywin = epnm(window,i);
            p = polyfit(twin, ywin, 1);   % quadratic
            epnmd0(i) = polyval(polyder(p), twin(1));
            epnm0(i) = epnm(1,i);
        end
    else
        idx_col1 = find(((texp-texp(idx_max))./tc < 1.25) == 1, 1, 'last');
        texp(1:idx_max-1) = [];
        R(1:idx_max-1) = [];

        for i = 1:15
            window = idx_max:idx_max+round((idx_col1-idx_max)/5);
            twin = t_nondim(window);
            ywin = epnm(window,i);
            p = polyfit(twin, ywin, 1);   % quadratic
            epnmd0(i) = polyval(polyder(p), twin(1));
        end
        for i = 1:length(epnmd0)
            epnm0(i) = mean(epnm(idx_max:idx_max+5,i));
        end
        epnm(1:idx_max-1, :) = [];
        texp = texp-texp(1);
    end
    R_nondim = R./Lc;
    % Find idx at collapse, only fit perturbations to collapse
    % if isequal(datatype, 'synthetic')
        [~, idx_col] = min(R(t_nondim < 1));
    % end
    epnm(idx_col+1:end, :) = [];

    % compute norms of the experimental data
    sR = norm(R_nondim);
    sEP = vecnorm(epnm);

    % equally weight the radial data with the sum of the perturbation
    % data
    wR = length(n);
    wEP = ones(1, length(n));

    % scaling factors
    aR = wR/sR;
    aEP = 1./sEP;
    
    count = 0;
    idx = [];
    for i = 1:length(n)
        if ((sEP(i)/size(epnm,1) >= 1e-3 || mean(abs(epnm(:,i))) >= 0.005) && i < 15)
            count = count + 1;
            epnm(:,count) = epnm(:,i);
        else
            idx = [idx, i];
        end
    end
    if count == 0
        continue
    end
    epnm(:,count+1:end) = [];
    n(idx) = [];
    m(idx) = [];
    aEP(idx) = [];
    epnm0(idx) = [];
    epnmd0(idx) = [];
    chinm(idx) = [];
    sEP(idx) = [];
    if isequal(datatype, 'exp')
        texp(R==0) = [];
        epnm(R(1:size(epnm,1))==0, :) = [];
        R_nondim(R==0) = [];
    end

    % create y_data vector which contains the radial data then all of
    % the perturbation data, both are non-dimensional
    y_data = [aR.*R_nondim; reshape(epnm.*aEP, [], 1)];

    % create x_data struct needed for the optimization
    x_data = struct('time', texp, 'n', n, 'Req', Req, 'chi', chinm, 'ST', ST, ...
        'rho', rho, 'Gqs', Gqs, 'epnmd0', epnmd0, 'epnm0', epnm0, 'R0', Rmax, ...
        'idx_col', idx_col, 'perturbed', perturbed, 'aR', aR, 'aEP', aEP);
    xsol = log10(xsole(s, :));
    sol = objfun(xsol, x_data);
    Rsol = sol(1:length(R_nondim))./aR;
    epnmsol = reshape(sol(length(R_nondim)+1:end), [size(epnm)])./aEP;
    figure
    plot((y_data), 'o')
    hold on 
    plot((sol), '-')
    % figure
    % plot(R_nondim, 'o')
    % hold on
    % plot(Rsol, '-')
    % 
    % figure
    % plot(epnm, 'o')
    % hold on
    % plot(epnmsol, '-')
    % yline(epnm0)
end
toc