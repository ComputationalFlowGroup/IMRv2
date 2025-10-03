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
ub = log10([1e1, 1e0]);

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
%load("Synthetic_data\synthetic_data_NSLIC_extracted_params.mat")
%load("../../../Experimental_data/Processed_data/LIC/ns_Jin_polyacr.mat")
 %datatype = 'synthetic';
%datatype = 'exp';


% Load In data

addpath Synthetic_data\ ../forward_solver/ ../common/
T1 = readtable('R1_anastas.csv');
T2 = readtable('R2_anastas.csv');

Gqs = 16.09e3; sigma = 0.00; rho = 1000;


texp = T2{:,1}.*1e-6;
R1exp = T1{:,2}.*1e-6;
R2exp = T2{:,2}.*1e-6;
Req = 50e-6;
Rmax = max(R1exp);
Rbar = (R1exp.*R2exp.^2).^(1/3);
ep = (R1exp-Rbar)./(Rbar*spherical_harmonic(cos(0),0,2,0));
tc = max(Rbar)*sqrt(1000/101325);

figure
plot(texp./tc, R1exp./max(Rbar), 'o')
hold on
plot(texp./tc, R2exp./max(Rbar), 'o')
plot(texp./tc, ep, '^')


perturbed = 1; % turns on characterization for surface perturbations, 1 is on


%% Loop over all datasets, extract material properties

% preallocate space for goodness of fit and extracted material properties
%xsole = zeros(length(kindata), length(lb));
%R2 = zeros(length(kindata),1);

tic
for s = 1%:length(kindata)
    % clear y_data x_data
    % %extract experiment data
    % exp = kindata{s};
    % 
    % %assign data from experiment to variables
    % texp = exp.time;
    % R = exp.R;
    % Req = exp.Req;
    % epnm = exp.epnm;
    % if isequal(datatype, 'synthetic')
    %     Gqs = exp.G;
    %     ST = exp.ST;
    %     rho = exp.rho;
    %     fps = exp.fps;
    % else
    %     Gqs = G; ST = sigma;
    %     epnm = epnm';
    %     fps = 2e6;
    %     texp = 0:length(R)-1;
    %     texp = texp./fps;
    % end
    % 
    % n = exp.n; m = exp.m;
    % 
    % % characteristic values
    % [Rmax, idx_max] = max(R);
    % Lc = Rmax;
    % rhoc = rho;
    % pc = 101325;
    % tc = Lc*sqrt(rhoc/pc);
    % 
    % % nondim values
    % t_nondim = texp./tc;
    % if size(R,1) == 1
    %     R = R';
    % end
    % 
    % 
    % % finite difference stencils for computing perturbation initial
    % % velcoty
    % fdstenc_for = [-49/20 6 -15/2 20/3 -15/4 6/5 -1/6]; % forward finite difference
    % fdstenc_cent = [-1/60 3/20 -3/4 0 3/4 -3/20 1/60]; % central finite difference
    % fdstencfor4 = [-25/12 4 -3 4/3 -1/4];
    % 
    % % smooth perturbation data from experiment
    % if isequal(datatype, 'exp')
    %     for i = 1:length(n)
    %         epnm(:,i) = sgolayfilt((fillmissing(epnm(:,i), 'spline')), 2, 7);
    %     end
    % end
    % 
    % % Center the data around t=0 at Rmax if needed, compute initial
    % % perturbation velocities, and chinm values
    % chinm = zeros(length(n), 1);
    % epnmd0 = chinm;
    % if idx_max == 1 || isequal(datatype, 'synthetic')
    %     idx_col1 = find(((texp-texp(idx_max))./tc < 1.25) == 1, 1, 'last');
    %     for i = 1:size(epnm,2)
    %         window = 1:1+round((idx_col1-idx_max)/8);
    %         twin = t_nondim(window);
    %         ywin = epnm(window,i);
    %         p = polyfit(twin, ywin, 1);   % quadratic
    %         epnmd0(i) = polyval(polyder(p), twin(1));
    %         epnm0(i) = epnm(1,i);
    %     end
    % else
    %     idx_col1 = find(((texp-texp(idx_max))./tc < 1.25) == 1, 1, 'last');
    %     texp(1:idx_max-1) = [];
    %     R(1:idx_max-1) = [];
    % 
    %     for i = 1:15
    %         chinm(i) = compute_chi(n(i), m(i));
    %         window = idx_max:idx_max+round((idx_col1-idx_max)/5);
    %         twin = t_nondim(window);
    %         ywin = epnm(window,i);
    %         p = polyfit(twin, ywin, 1);   % quadratic
    %         epnmd0(i) = polyval(polyder(p), twin(1));
    %     end
    %     for i = 1:length(epnmd0)
    %         epnm0(i) = mean(epnm(idx_max:idx_max+5,i));
    %     end
    %     epnm(1:idx_max-1, :) = [];
    %     texp = texp-texp(1);
    % end
    % R_nondim = R./Lc;
    % % Find idx at collapse, only fit perturbations to collapse
    % % if isequal(datatype, 'synthetic')
    %     [~, idx_col] = min(R(t_nondim < 1));
    % % end
    % epnm(idx_col+1:end, :) = [];
    % 
    % % compute norms of the experimental data
    % sR = norm(R_nondim);
    % sEP = vecnorm(epnm);
    % 
    % % equally weight the radial data with the sum of the perturbation
    % % data
    % wR = length(n);
    % wEP = ones(1, length(n));
    % 
    % % scaling factors
    % aR = wR/sR;
    % aEP = 1./sEP;
    % 
    % count = 0;
    % idx = [];
    % for i = 1:length(n)
    %     if ((sEP(i)/size(epnm,1) >= 1e-3) && i < 15 && n(i) > 1)
    %         count = count + 1;
    %         epnm(:,count) = epnm(:,i);
    %     else
    %         idx = [idx, i];
    %     end
    % end
    % if count == 0
    %     continue
    % end
    % epnm(:,count+1:end) = [];
    % n(idx) = [];
    % m(idx) = [];
    % aEP(idx) = [];
    % epnm0(idx) = [];
    % epnmd0(idx) = [];
    % chinm(idx) = [];
    % sEP(idx) = [];
    % if isequal(datatype, 'exp')
    %     texp(R==0) = [];
    %     epnm(R(1:size(epnm,1))==0, :) = [];
    %     R_nondim(R==0) = [];
    % end
    % 
    % % create y_data vector which contains the radial data then all of
    % % the perturbation data, both are non-dimensional
    % y_data = [aR.*R_nondim; reshape(epnm.*aEP, [], 1)];
    % 
    % idx_col = size(epnm,1);
    n = 2;
    Req = 50e-6; ST = 0; epnmd0 = 0;
    epnm0 = ep(1); Rmax = max(Rbar); idx_col = length(R1exp); perturbed  = 1;
    aR = 1; aEP = 1;
    y_data = [Rbar./max(Rbar); ep];
    % create x_data struct needed for the optimization
    x_data = struct('time', texp, 'n', n, 'Req', Req,'ST', sigma, ...
        'rho', rho, 'Gqs', Gqs, 'epnmd0', epnmd0, 'epnm0', epnm0, 'R0', Rmax, ...
        'idx_col', idx_col, 'perturbed', perturbed, 'aR', aR, 'aEP', aEP);
    % figure
    % plot(texp(1:size(epnm,1))./tc, epnm, 'o--')
    % hold on
    % for i = 1: length(epnmd0)
    %     % yline(epnm0(i))
    %     plot(texp(1:size(epnm,1))./tc,epnmd0(i).*texp(1:size(epnm,1))./tc+epnm0(i), '-')
    % end
    % 
    % 
    %  % ------------- Perform Optimization --------------- %
    % 
     %define objective function that outputs the y_data_sim
     objfun = @(params, x_data) f_run_fd_IMR(10.^params, x_data);

    % create scalar function that computes sum or squared error
    objfun_scalar = @(params, x_data, y_data) sqrt(sum((y_data - objfun(params, x_data)).^2))/norm(y_data);
    % 
    % % create loss function that can be read by bayesopt
    % lossfun = @(T) objfun_scalar([T.Alpha, T.Mu], x_data, y_data);
    % 
    % % Start with bayesopt to quickly narrow the parameter space to a few
    % % candidates
    % results = bayesopt(lossfun, bayesoptVars, ...
    %     'MaxObjectiveEvaluations', 50, ...
    %     'UseParallel', true, ...
    %     'IsObjectiveDeterministic', true, ...
    %     'AcquisitionFunctionName', 'lower-confidence-bound', ...
    %     'NumSeedPoints', 5, ...
    %     'ExplorationRatio', 0.75, ...
    %     'Verbose', 1, ...
    %     'PlotFcn', {@plotObjectiveModel, @plotMinObjective});
    % 
    % % results is the bayesopt output
    % OT = results.ObjectiveTrace;      % numeric vector
    % XT = results.XTrace;              % table of tried points (log-scale if you used log variables)
    % 
    % % number of restarts for lsqcurvefit, extract top n_starts performers
    % % from bayes opt to feed into lsqcurvefit
    % n_starts = 3;
    % [~,I] = sort(OT, 'ascend');
    % k = min(n_starts, numel(I));           
    % topk_table = XT(I(1:k),:);
    % topk_mat = table2array(topk_table);
    % starts = topk_mat;
    % 
    % % create start points for the optimizer
    % custom_starts = CustomStartPointSet(starts);
    % 
    %     % set up lsqcurvefit multirun optimization problem
    % problem = createOptimProblem('lsqcurvefit', ...
    %     'x0', starts(1,:), ...
    %     'objective', objfun, ...
    %     'xdata', x_data, ...
    %     'ydata', y_data, ...
    %     'lb', lb, ...
    %     'ub', ub, ...
    %     'options', lsqopts);
    % ms = MultiStart('UseParallel', true, 'StartPointsToRun','all');
    % 
    % % run the optimizer, extract the best result
    % [xsol, fval, exitflag, output, solutions] = run(ms, problem, custom_starts);

    % evaluate goodness of fit
    xsol = log10([1e-12 0.1]);

    sol = objfun(xsol, x_data);
    R2(s) = 1-sum((sol-y_data).^2)/sum((y_data-mean(y_data)).^2);
    xsole(s, :) = 10.^xsol;
    sol = reshape(sol, [length(R1exp),2]);
    tc = max(Rbar)*sqrt(1000/101325);
    Rbarsim = sol(:,1);
    epsim = sol(:,2);
    R1sim = Rbarsim.*(1+epsim*spherical_harmonic(cos(0),0, 2,0));
    R2sim = Rbarsim.*(1+epsim*spherical_harmonic(cos(pi/2),0, 2,0));
    hold on
    plot(texp./tc, R1sim, '-')
    plot(texp./tc, R2sim, '--')
    plot(texp./tc, epsim, '.-')
    xlabel('t/tc ', 'Interpreter','latex', FontSize=20)
    ylabel('$R(t)/R_{max}$', 'Interpreter','latex', FontSize=20)
    set(gca, 'TickLabelInterpreter', 'latex')
    % figure 
    % plot(sol, '-')
    % hold on
    % plot(y_data,'o')
end
toc

%% Histograms for extracted material properties from multiple datasets
 %clear all
 %clc
% 
load("Synthetic_data\synthetic_data_NSLIC_extracted_params.mat")
%%
% load("../../../Experimental_data/Processed_data/LIC/ns_Jin_polyacr_optimized_newbds_fixedICs_new.mat")
load("NSLIC_processed.mat")

exps =[3 5 9 11];

addpath ../../../cmap/

cmap = parula(4); 

alph_e = xsole(exps ,1);
mu_e = xsole(exps ,2);

figure
countnew = length(alph_e);
% ===== Left subplot (alph) =====
subplot(1,2,1)
edgesalph = logspace(1.1*log10(min(alph_e)), 0.9*log10(max(alph_e)), ceil(sqrt(countnew))+1);
histogram(alph_e, edgesalph, 'FaceColor', cmap(1,:), 'EdgeColor', 'k')
hold on
xline(mean(alph_e), '--k', 'LineWidth', 2)
xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', 24)
ylabel('Count', 'FontSize', 24, 'Interpreter', 'latex')
set(gca, ...
    'TickLabelInterpreter', 'latex', ...
    'FontSize', 16, ...
    'XScale', 'log', ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off')

% Force evenly spaced ticks (log spacing but numeric labels)
xticksalph = logspace(log10(0.75*min(edgesalph)), log10(1.25*max(edgesalph)), 4);

grid on
xlim([10^(floor(log10(min(edgesalph)))) 10^(ceil(log10(max(edgesalph))))])
ylim([0 4])

% ===== Right subplot (mu) =====
subplot(1,2,2)
edgesMu = logspace(0.9*log10(min(mu_e.*1e3)), 1.1*log10(max(mu_e.*1e3)), ceil(sqrt(countnew))+1);
histogram(mu_e.*1e3, edgesMu, 'FaceColor', cmap(3,:), 'EdgeColor', 'k')
hold on
xline(mean(mu_e).*1e3, '--k', 'LineWidth', 2)
xlabel('$\mu$ [mPa$\cdot$s]', 'Interpreter', 'latex', 'FontSize', 24)
set(gca, ...
    'TickLabelInterpreter', 'latex', ...
    'FontSize', 16, ...
    'XScale', 'log', ...
    'XMinorTick', 'off', ...
    'YMinorTick', 'off')
grid on
xlim([10^(floor(log10(min(edgesMu)))) 10^(ceil(log10(max(edgesMu))))])
ylim([0 4])


%% Plot results
clear all
clc

% load("../../../Experimental_data/Processed_data/LIC/ns_Jin_polyacr_optimized.mat")
load("Synthetic_data\synthetic_data_NSLIC_extracted_params.mat")

addpath Synthetic_data\ ../forward_solver/ ../common/


%%
addpath ../../../cmap/

%exps_extracted_whole = [2 3 4 5 6 7 8 9 10 11];
exps_extracted_whole = [2 3 4 5 8 11];
exps_extracted = [4, 11];
exps_modes{1} = [3 4 6 7 9 10];
exps_modes{2} = [2 3];
% exps_extracted = [2 5 9 11];
% exps_modes{1} = [3 4 6 7 9 10];
% exps_modes{2} = [3 4 6 7 9];
% exps_modes{3} = [3 4 6 7 10 9 12];
% exps_modes{4} = [2 3];


% preallocate space for goodness of fit and extracted material properties
%xsole = zeros(length(kindata), length(lb));
%R2 = zeros(length(kindata),1);
%close all
tic
for s = 3;%[2 5 9]% 9 11]%exps_extracted%1:length(kindata)
    % if R2(s) < 0.96
    %     continue
    % end
    s
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

    count = 0;
    idx = [];
    for i = 1:length(n)
        if ((sEP(i)/size(epnm,1) >= 1e-6) && i < 15 && n(i) > 1)
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
    idx_col = size(epnm,1);
    
    

    % 
    % if count == 0
    %     continue
    % end
    % epnm(:,count+1:end) = [];
    % n(idx) = [];
    % m(idx) = [];
    % aEP(idx) = [];
    % epnm0(idx) = [];
    % epnmd0(idx) = [];
    % chinm(idx) = [];
    % sEP(idx) = [];

%     [hasSpec, jSpec] = ismember(s, exps_extracted);
% if hasSpec
%     allowed_n = exps_modes{jSpec};
% else
%        for i = 1:length(n)
%         if ((sEP(i)/size(epnm,1) >= 1e-3) && i < 15 && n(i) > 1)
%             count = count + 1;
%             epnm(:,count) = epnm(:,i);
%         else
%             idx = [idx, i];
%         end
%        end
%     if count == 0
%         continue
%     end
%     epnm(:,count+1:end) = [];
%     n(idx) = [];
%     m(idx) = [];
%     aEP(idx) = [];
%     epnm0(idx) = [];
%     epnmd0(idx) = [];
%     chinm(idx) = [];
%     sEP(idx) = [];
% end
% 
% % Build a single mask over columns of epnm / entries of n
% keepMask = ismember(n, allowed_n);
% 
% % If nothing to keep, skip this experiment cleanly
% if ~any(keepMask)
%     fprintf('Experiment %d: no allowed modes -> skipping.\n', s);
%     continue
% end
% 
% % Apply mask to ALL per-mode quantities (columns)
% n      = n(keepMask);
% m      = m(keepMask);
% epnm   = epnm(:, keepMask);
% if exist('epnm0','var'),   epnm0   = epnm0(keepMask);   end
% if exist('epnmd0','var'),  epnmd0  = epnmd0(keepMask);  end
% if exist('chinm','var'),   chinm   = chinm(keepMask);   end
% 
% % Rows count for later reshapes (unchanged by column selection)
% idx_col = size(epnm,1);
% 
% % === NOW compute norms/weights/scales using the filtered columns ===
% sR  = norm(R_nondim);           % unchanged (radial series)
% sEP = vecnorm(epnm);            % 1 x numModes
% wR  = numel(n);
% wEP = ones(1, numel(n));
% aR  = wR / sR;
% aEP = 1 ./ sEP;                 % per-column scale
    n


    % create y_data vector which contains the radial data then all of
    % the perturbation data, both are non-dimensional
    y_data = [aR.*R_nondim; reshape(epnm.*aEP, [], 1)];

    idxcol1 = size(epnm,1);
    tcol = texp(idx_col);
    texp1 = texp;
    texp = linspace(0, texp(end), 1000);
    [~, idx_col] = min(abs(texp-tcol));

    % create x_data struct needed for the optimization
    x_data = struct('time', texp, 'n', n, 'Req', Req, 'chi', chinm, 'ST', ST, ...
        'rho', rho, 'Gqs', Gqs, 'epnmd0', epnmd0, 'epnm0', epnm0, 'R0', Rmax, ...
        'idx_col', idx_col, 'perturbed', perturbed, 'aR', aR, 'aEP', aEP);
    xsol = log10(xsole(s, :));
    xsol(1) = log10(0.1);
    sol = objfun(xsol, x_data);
    Rsol = sol(1:length(texp))./aR;
    epnmsol = reshape(sol(length(texp)+1:end), [idx_col, length(n)])./aEP;
%     % figure
%     % plot(abs(y_data), 'o')
%     % hold on 
%     % plot(abs(sol), '-')
% 
    cmap = viridis(length(n)+4);


    rows = 1 + ceil(length(n)/3);
    columns = 3;
    if length(n) < 3
        columns = length(n);
    end

    figure
    a = tiledlayout(rows,columns,'TileSpacing','compact','Padding','compact');

    % --- Top row: first plot spans all columns ---
    ax = nexttile([1 columns]);      % span 1 row by all columns
    plot(ax, texp1(1:2:end)./tc, R_nondim(1:2:end), 'o', 'Color', cmap(1,:), 'MarkerFaceColor',cmap(1,:), 'MarkerSize',8)
    hold on
    plot(ax, texp./tc, Rsol, '-', 'LineWidth',3, 'Color', cmap(5,:))
    set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 18)
    xlabel('$t/t_c$', 'Interpreter','latex', 'FontSize',20)
    ylabel('$R/R_{\textrm{max}}$', 'Interpreter', 'latex','FontSize',20)

for i = 1:length(n)
    ax = nexttile;
     plot(ax,texp1(1:2:size(epnm,1))./tc, epnm(1:2:end,i), 'o', 'Color', cmap(i+1,:), 'MarkerFaceColor',cmap(i+1,:), 'MarkerSize',8);
     hold on
     plot(ax, texp(1:idx_col)./tc, epnmsol(:,i), '-', 'LineWidth',3, 'Color', cmap(i+3,:))
     set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 18)
     xlabel('$t/t_c$', 'Interpreter','latex', 'FontSize',20)
     ylabel(['$\epsilon_{' num2str(n(i)) '}$'], 'Interpreter', 'latex','FontSize',20)
end





% 
%      %figure
%      %plot(R_nondim, 'o')
%      %hold on
%      %plot(Rsol, '-')
%     % 
%      %figure
%      %plot(abs(epnm), 'o')
%      %hold on
%      %plot(abs(epnmsol), '-')
%     % yline(epnm0)
end
toc
%%

    
