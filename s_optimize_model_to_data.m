%% Optimize anisotropic IMR model to experimental data
clear
clc
close all

[scriptDir, projectDir] = locateProjectPaths();
addpath(fullfile(scriptDir, 'src', 'common'));
addpath(fullfile(scriptDir, 'src', 'forward_solver'));
addpath(fullfile(scriptDir, 'src', 'characterization'));

%% User settings
dataFile = fullfile(projectDir, 'data', 'SicongJinChicken', ...
    'PVA_Rt_data', 'processed_11.mat');

material = "PVA";
maxmode = 22;
polyOrder = 3;
windowPts = 15;   % odd integer: 3, 5, 7, ...
icVelocityWindowPts = 15;  % forward polynomial derivative window from Rmax
icVelocityPolyOrder = 3;
opt.fit.NumPerturbationModes = 7;

% Loss priority weights after per-trace normalization. Radial receives
% RadialWeight, and mode n receives
% ModeBaseWeight*(ModeReference/(n + ModeWeightOffset))^ModeWeightPower.
opt.loss.RadialWeight = 1;
opt.loss.ModeBaseWeight = 1;
opt.loss.ModeReference = 2;
opt.loss.ModeWeightPower = 1;
opt.loss.ModeWeightOffset = 0;

% Parameter bounds. G and mu are optimized in log10-space by default.
opt.bounds.G = [5e4, 5e6];
opt.bounds.alph = [1e-3, 5];
opt.bounds.mu = [5e-2, 5e-1];
opt.bounds.ani = [0, 5; ...
                 0, 5];

% Set any entry to a finite value to remove that parameter from the
% optimizer and hold it fixed. Leave entries as NaN to optimize them.
% Examples:
%   opt.fixed.alph = 0;
%   opt.fixed.ani = [0, 0];
%   opt.fixed.ani = [NaN, 0];  % optimize ani1, fix ani2 = 0
opt.fixed.G = NaN;
opt.fixed.alph = NaN;
opt.fixed.mu = NaN;
opt.fixed.ani = [NaN, NaN];

opt.logScale.G = true;
opt.logScale.alph = true;
opt.logScale.mu = true;
opt.logScale.ani = [false, false];

% Forward-solver controls.
opt.sim.tsteps = 3000;          % fallback only; optimizer feeds experimental times
opt.sim.RelTol = 1e-4;
opt.sim.AbsTol = 1e-5;
opt.sim.Nt = 75;
opt.sim.Method = 23;
opt.sim.MaxWallTime = 120;       % seconds per simulation
opt.sim.UseHardTimeout = false; % true can fail if background workers do not inherit paths
opt.sim.TimeoutPollInterval = 1;
opt.sim.FailurePenalty = 1e5;   % scalar BO loss for failed/timeout runs
opt.sim.UseLogLoss = true;      % Bayes opt minimizes log10(normalized loss)
opt.sim.LogLossFloor = 1e-12;   % prevents log10(0)
opt.sim.TimeMatchTolerance = 1e-8; % dimensionless tolerance for t == t_exp
opt.sim.AllowNearestTimeExtraction = true; % no interpolation; nearest returned samples only
opt.sim.NearestTimeTolerance = 1e-6; % dimensionless tolerance for nearest returned samples
opt.sim.VerifyTimeExtraction = false; % final best fit is always verified
opt.sim.PrintFailures = true;
opt.sim.PrintSuccess = false;
opt.sim.DiagnosticLogFile = fullfile(scriptDir, ...
    'model_to_data_eval_diagnostics.tsv');
opt.sim.DiagnosticLogAppend = false;

% Bayes-opt controls.
opt.bayes.MaxObjectiveEvaluations = 500;
opt.bayes.NumSeedPoints = 50;
opt.bayes.UseLatinHypercubeInitialX = true;
opt.bayes.InitialRegionFraction = 1; % lower half of each optimizer bound range
opt.bayes.UseParallel = true;  % use true for larger budgets or an open pool
opt.bayes.NumWorkers = [];     % [] uses current/default pool size
opt.bayes.ParallelPoolProfile = ''; % '' uses MATLAB's default profile
opt.bayes.RestartParallelPool = false;
opt.bayes.MinWorkerUtilization = []; % [] keeps all pool workers busy
opt.bayes.ParallelMethod = 'clipped-model-prediction';
opt.bayes.IsObjectiveDeterministic = true;
opt.bayes.AcquisitionFunctionName = 'lower-confidence-bound';
opt.bayes.ExplorationRatio = 0.65;
opt.bayes.Verbose = 1;
opt.bayes.PlotFcn = {};         % plots add overhead during fast searches
opt.bayes.RunClientPreflight = true;
opt.bayes.RunWorkerPreflight = true;
opt.bayes.FallbackToSerialOnWorkerPreflightFailure = true;

% Optional local refinement from the best Bayes-opt points.
opt.refine.Enabled = true;     % can be expensive: finite differences call many simulations
opt.refine.NumStarts = 15;
opt.refine.Display = 'iter-detailed';
opt.refine.UseParallel = true;
opt.refine.MaxFunctionEvaluations = 100;
opt.refine.MaxIterations = 25;
opt.refine.OptimalityTolerance = 1e-4;
opt.refine.FunctionTolerance = 1e-4;
opt.refine.StepTolerance = 1e-5;
opt.refine.FiniteDifferenceType = 'forward';
opt.refine.Algorithm = 'interior-point';

opt.randomSeed = 1;
opt.outputFile = fullfile(scriptDir, 'optimized_model_to_data.mat');
opt.plotInitialConditionCheck = true;

%% Load and process data
if ~isempty(opt.randomSeed)
    rng(opt.randomSeed, 'twister');
end

if ~isfile(dataFile)
    error('Could not find the data file: %s', dataFile);
end
load(dataFile)

pxpermicron = 3.2;

if material == "PVA"
    tstepdt = 5e-7;
elseif material == "chicken"
    tstepdt = 1e-6;
end

expR = amp_extract_fft(1, :) .* 1e-6 .* pxpermicron;
texp = (0:numel(expR)-1) .* tstepdt;

amps_og = amp_extract_fft(3:end, :) ./ expR .* 1e-6 .* pxpermicron;
amps_og = fillNonfiniteTimeRows(amps_og);
Req = expR(end);

maxmode = min(maxmode, size(amps_og, 1) + 1);
windowPts = min(windowPts, 2 * floor((size(amps_og, 2) - 1) / 2) + 1);
if mod(windowPts, 2) == 0
    windowPts = windowPts - 1;
end

amps = sgolayfilt(amps_og(1:maxmode-1, :), polyOrder, windowPts, [], 2);

% Extract initial conditions at maximum radius.
[Rmax, maxidx] = max(expR);
tc = Rmax * sqrt(1000 / 101325);

epnm0 = amps(:, maxidx);
epnmd0 = computeInitialModeVelocities(amps, texp, maxidx, tc, ...
    icVelocityWindowPts, icVelocityPolyOrder);
eqWindow = max(1, size(amps, 2)-20):size(amps, 2);
epnmeq = mean(amps(:, eqWindow), 2);

% Fit only the data after Rmax because the forward simulation starts there.
fitIdx = maxidx:numel(texp);
tfit_nd = (texp(fitIdx) - texp(maxidx)) ./ tc;
R_data = expR(fitIdx).' ./ Rmax;
ep_data = amps(:, fitIdx).';
tf_nd = max(tfit_nd);
[firstCollapseIdx, collapseInfo] = findFirstCollapseIndex(R_data);
epFitIdx = 1:2*firstCollapseIdx;
firstCollapseTimeNd = tfit_nd(firstCollapseIdx);
firstCollapseTimeSeconds = firstCollapseTimeNd * tc;

modeRows = 3:maxmode+1;
n = mode_extract_fft(modeRows, 10);
n = n(:).';
m = zeros(size(n));

% Train only against the largest-energy perturbation modes before first
% collapse, while still simulating every retained mode for held-out testing.
nmodes = size(ep_data, 2);
modeEnergy = sum(ep_data(epFitIdx, :).^2, 1);
[~, modeEnergyOrder] = sort(modeEnergy, 'descend');
nFitModes = min(opt.fit.NumPerturbationModes, nmodes);
fitModeIdx = sort(modeEnergyOrder(1:nFitModes));
testModeIdx = setdiff(1:nmodes, fitModeIdx, 'stable');

% Normalize each trace, then apply priority weights: radial first, then
% lower-order perturbation modes before higher-order modes.
lossWeights = makePriorityLossWeights(n, opt.loss);
aR = lossWeights.radial / norm(R_data);
sEP = vecnorm(ep_data(epFitIdx, :), 2, 1);
sEP(sEP < eps) = 1;
aEP = lossWeights.modes ./ sEP;
yDataAllFitModes = [aR .* R_data; reshape(ep_data(epFitIdx, fitModeIdx) .* ...
    aEP(fitModeIdx), [], 1)];

xDataAll = struct( ...
    'Rmax', Rmax, ...
    'Req', Req, ...
    'epnm0', epnm0, ...
    'epnmd0', epnmd0, ...
    'epnmeq', epnmeq, ...
    'n', n, ...
    'm', m, ...
    'sig', 0.056, ...
    'p_a', -0 * 101325, ...
    'f_a', 50e3, ...
    'rho', 1000, ...
    'p8', 101325, ...
    'tf_nd', tf_nd, ...
    'tfit_nd', tfit_nd(:), ...
    'R_data', R_data, ...
    'ep_data', ep_data, ...
    'epFitIdx', epFitIdx, ...
    'firstCollapseIdx', firstCollapseIdx, ...
    'firstCollapseTimeNd', firstCollapseTimeNd, ...
    'firstCollapseTimeSeconds', firstCollapseTimeSeconds, ...
    'collapseInfo', collapseInfo, ...
    'fitModeIdx', fitModeIdx, ...
    'testModeIdx', testModeIdx, ...
    'modeEnergy', modeEnergy, ...
    'lossWeights', lossWeights, ...
    'aR', aR, ...
    'aEP', aEP, ...
    'y_data', yDataAllFitModes, ...
    'ultra', false);

xDataOpt = xDataAll;
yData = xDataOpt.y_data;

paramSpec = buildParameterSpec(opt);
bayesoptVars = buildBayesoptVariables(paramSpec);
[lbOpt, ubOpt] = optimizerBounds(paramSpec);
if ~isempty(bayesoptVars)
    opt.bayes = prepareBayesoptParallelPool(opt.bayes, scriptDir);
end
if opt.bayes.UseLatinHypercubeInitialX
    opt.bayes.InitialX = makeLowerRegionLatinHypercubeInitialX( ...
        paramSpec, opt.bayes.NumSeedPoints, opt.bayes.InitialRegionFraction);
end
initializeDiagnosticLog(opt.sim);

fprintf('Training perturbation modes: %s\n', num2str(n(fitModeIdx)));
fprintf('Held-out perturbation modes: %s\n', num2str(n(testModeIdx)));
printLossWeights(n, fitModeIdx, lossWeights);
printParameterSpec(paramSpec);
fprintf(['First collapse (raw radius minimum) at post-Rmax sample ', ...
    '%d/%d: t* = %.6g, dt = %.6g us, R/Rmax = %.6g\n'], ...
    firstCollapseIdx, numel(tfit_nd), firstCollapseTimeNd, ...
    firstCollapseTimeSeconds * 1e6, R_data(firstCollapseIdx));
fprintf(['  Sustained rebound confirmed through sample %d: t* = %.6g, ', ...
    'rise in R/Rmax = %.6g\n'], collapseInfo.confirmationIdx, ...
    tfit_nd(collapseInfo.confirmationIdx), collapseInfo.reboundRise);
fprintf('Perturbation loss uses %d samples from t* = %.6g to %.6g.\n', ...
    numel(epFitIdx), tfit_nd(epFitIdx(1)), tfit_nd(epFitIdx(end)));
if isfield(opt.bayes, 'InitialX') && ~isempty(opt.bayes.InitialX)
    disp('Bayes-opt initial Latin-hypercube points:')
    disp(opt.bayes.InitialX)
end

if opt.plotInitialConditionCheck
    plotInitialConditionCheck(texp, maxidx, amps_og, amps, epnm0, epnmd0, tc);
    plotCollapseDetection(tfit_nd, R_data, collapseInfo);
end

%% Bayesian optimization
tic
lossfun = @(T) f_optimize_model_to_data_loss(T, xDataOpt, paramSpec, opt.sim);
opt.bayes = runOptimizerPreflight(opt.bayes, xDataOpt, paramSpec, opt.sim);

if isempty(bayesoptVars)
    results = [];
    bestZFromBayes = [];
    fvalFromBayes = lossfun([]);
else
    bayesArgs = makeBayesoptArgs(opt.bayes);
    results = bayesopt(lossfun, bayesoptVars, bayesArgs{:});
    bestZFromBayes = tableToOptimizerMatrix(results.XAtMinObjective, ...
        paramSpec);
    fvalFromBayes = results.MinObjective;
end

% Optional local refinement
refineOutput = [];
solutions = [];
if opt.refine.Enabled && ~isempty(bayesoptVars)
    starts = bestBayesStarts(results, paramSpec, opt.refine.NumStarts);
    customStarts = CustomStartPointSet(starts);
    lsqopts = makeLsqOptions(opt.refine);
    objfun = @(z, xDataIn) f_optimize_model_to_data_predict( ...
        z, xDataIn, paramSpec, opt.sim);

    problem = createOptimProblem('lsqcurvefit', ...
        'x0', starts(1, :), ...
        'objective', objfun, ...
        'xdata', xDataOpt, ...
        'ydata', yData, ...
        'lb', lbOpt, ...
        'ub', ubOpt, ...
        'options', lsqopts);

    ms = MultiStart('UseParallel', opt.refine.UseParallel, ...
        'StartPointsToRun', 'all');
    [bestZ, fval, exitflag, refineOutput, solutions] = run(ms, problem, ...
        customStarts);
else
    bestZ = bestZFromBayes;
    fval = fvalFromBayes;
    exitflag = NaN;
end
toc

save('../optimized_data/Jin_data/PVA/optimized_11.mat')

%% Evaluate and plot best fit

% load('../optimized_data/Jin_data/PVA/optimized_22.mat')

bestParams = f_unpack_model_to_data_params(bestZ, paramSpec);
[bestY, bestRunInfo, bestSimOpt] = f_optimize_model_to_data_predict(bestZ, ...
    xDataOpt, paramSpec, opt.sim);
fullRunInfo = bestRunInfo;
bestSim = bestSimOpt;
xData = xDataAll;  % saved alias for compatibility with older analysis code
timeVerification = makeTimeExtractionVerification(bestSim, xDataAll, ...
    bestRunInfo);

bestLoss = sqrt(sum((yData - bestY).^2)) / norm(yData);
fitR2 = 1 - sum((bestY - yData).^2) / sum((yData - mean(yData)).^2);
trainModeLoss = computePerturbationSubsetLoss(bestSim, xDataAll, ...
    xDataAll.fitModeIdx);
heldoutModeLoss = computePerturbationSubsetLoss(bestSim, xDataAll, ...
    xDataAll.testModeIdx);
allModeLoss = computePerturbationSubsetLoss(bestSim, xDataAll, 1:nmodes);

fprintf('\nBest fit\n');
fprintf('  G     = %.6g Pa\n', bestParams.G);
fprintf('  alph  = %.6g\n', bestParams.alph);
fprintf('  mu    = %.6g Pa*s\n', bestParams.mu);
fprintf('  ani   = [%.6g %.6g]\n', bestParams.ani(1), bestParams.ani(2));
fprintf('  loss  = %.6g\n', bestLoss);
fprintf('  R2    = %.6g\n', fitR2);
fprintf('  train mode loss   = %.6g\n', trainModeLoss);
fprintf('  held-out mode loss = %.6g\n', heldoutModeLoss);
fprintf('  all mode loss     = %.6g\n', allModeLoss);
fprintf('  radial loss samples       = %d\n', bestRunInfo.nRadialLossTimes);
fprintf('  perturbation loss samples = %d (through t* = %.6g)\n', ...
    bestRunInfo.nPerturbationLossTimes, bestRunInfo.firstCollapseTimeNd);
fprintf('  max requested/returned simulation time mismatch = %.3g\n', ...
    bestRunInfo.maxRequestedReturnedTimeMismatch);
fprintf('  max radial time extraction mismatch = %.3g\n', ...
    timeVerification.maxRadialAbsDt);
fprintf('  max perturbation time extraction mismatch = %.3g\n', ...
    timeVerification.maxPerturbationAbsDt);
fprintf('  perturbation rows used    = %d:%d of %d\n', ...
    timeVerification.firstPerturbationRow, ...
    timeVerification.lastPerturbationRow, numel(xDataAll.tfit_nd));

plotOptimizedFit(bestSim, xDataAll, bestParams);

save(opt.outputFile, 'opt', 'paramSpec', 'xData', 'xDataAll', ...
    'xDataOpt', 'results', 'bestZ', 'bestParams', 'bestLoss', ...
    'fitR2', 'bestRunInfo', 'fullRunInfo', 'bestSimOpt', 'bestSim', ...
    'timeVerification', 'trainModeLoss', 'heldoutModeLoss', ...
    'allModeLoss', 'fval', ...
    'exitflag', 'refineOutput', 'solutions');

%% Local helper functions
function [scriptDir, projectDir] = locateProjectPaths()
    candidates = {};
    candidates = addCandidate(candidates, fileparts(mfilename('fullpath')));
    candidates = addCandidate(candidates, fileparts(which('s_optimize_model_to_data')));
    candidates = addCandidate(candidates, fileparts(which('f_call_IMRv2_exp')));
    candidates = addCandidate(candidates, pwd);
    candidates = addCandidate(candidates, fullfile(pwd, 'IMRv2'));
    candidates = addCandidate(candidates, fileparts(pwd));
    candidates = addCandidate(candidates, fullfile(fileparts(pwd), 'IMRv2'));

    for ii = 1:numel(candidates)
        candidate = candidates{ii};
        if isImrScriptDir(candidate)
            thisProjectDir = fileparts(candidate);
            if isfolder(fullfile(thisProjectDir, 'data'))
                scriptDir = candidate;
                projectDir = thisProjectDir;
                return
            end
        end
    end

    for ii = 1:numel(candidates)
        candidate = candidates{ii};
        if isfolder(fullfile(candidate, 'IMRv2')) && ...
                isImrScriptDir(fullfile(candidate, 'IMRv2')) && ...
                isfolder(fullfile(candidate, 'data'))
            scriptDir = fullfile(candidate, 'IMRv2');
            projectDir = candidate;
            return
        end
    end

    error(['Could not locate the Anisotropic_material_IMR project root. ', ...
        'Run this script from the project root or the IMRv2 folder.']);
end

function candidates = addCandidate(candidates, candidate)
    if isempty(candidate) || ~ischar(candidate)
        return
    end
    candidate = char(candidate);
    if isempty(candidate) || ~isfolder(candidate)
        return
    end
    candidate = char(strrep(candidate, '/', filesep));
    if ~any(strcmpi(candidates, candidate))
        candidates{end+1} = candidate;
    end
end

function tf = isImrScriptDir(candidate)
    tf = ~isempty(candidate) && isfolder(candidate) && ...
        isfile(fullfile(candidate, 'f_call_IMRv2_exp.m')) && ...
        isfolder(fullfile(candidate, 'src', 'forward_solver'));
end

function epnmd0 = computeInitialModeVelocities(amps, texp, maxidx, tc, ...
    windowPts, polyOrder)
    epnmd0 = zeros(size(amps, 1), 1);
    if nargin < 5 || isempty(windowPts)
        windowPts = 15;
    end
    if nargin < 6 || isempty(polyOrder)
        polyOrder = 3;
    end

    windowIdx = localForwardWindow(size(amps, 2), maxidx, windowPts);
    tstar = (texp(windowIdx) - texp(maxidx)) ./ tc;
    fitOrder = min(polyOrder, numel(windowIdx) - 1);
    if fitOrder < 1
        return
    end

    tstar = tstar(:);
    finiteTime = isfinite(tstar);
    for ii = 1:size(amps, 1)
        y = amps(ii, windowIdx).';
        if all(finiteTime) && all(isfinite(y))
            p = polyfit(tstar, y, fitOrder);
        else
            fitRows = finiteTime & isfinite(y);
            localFitOrder = min(polyOrder, nnz(fitRows) - 1);
            if localFitOrder < 1
                continue
            end
            p = polyfit(tstar(fitRows), y(fitRows), localFitOrder);
        end
        epnmd0(ii) = polyval(polyder(p), 0);
    end
end

function values = fillNonfiniteTimeRows(values)
    if all(isfinite(values(:)))
        return
    end

    sampleIdx = 1:size(values, 2);
    for ii = 1:size(values, 1)
        row = values(ii, :);
        good = isfinite(row);
        if all(good)
            continue
        end
        if nnz(good) < 2
            error(['Cannot repair nonfinite amplitude row %d; at least ', ...
                'two finite samples are required.'], ii);
        end
        row(~good) = interp1(sampleIdx(good), row(good), ...
            sampleIdx(~good), 'linear', 'extrap');
        values(ii, :) = row;
    end
end

function windowIdx = localForwardWindow(nSamples, startIdx, windowPts)
    windowPts = max(3, round(windowPts));
    if mod(windowPts, 2) == 0
        windowPts = windowPts - 1;
    end
    windowPts = min(windowPts, nSamples);
    startIdx = min(max(1, startIdx), nSamples);
    lastIdx = min(nSamples, startIdx + windowPts - 1);
    windowIdx = startIdx:lastIdx;
end

function [firstCollapseIdx, info] = findFirstCollapseIndex(R_data)
    R_data = R_data(:);
    nSamples = numel(R_data);
    info = struct('method', 'sustained-raw-rebound', ...
        'confirmationIdx', nSamples, 'reboundRise', 0, ...
        'collapseDrop', 0, 'usedFallback', false);
    if nSamples < 3
        firstCollapseIdx = numel(R_data);
        return
    end

    if any(~isfinite(R_data))
        error('Radius data contain nonfinite values; cannot locate collapse.');
    end

    initialRadius = R_data(1);
    totalDrop = max(0, initialRadius - min(R_data));
    minCollapseDrop = max(0.05 * max(abs(initialRadius), eps), ...
        0.10 * totalDrop);
    runningMin = cummin(R_data);
    noiseTol = max(1e-10, 0.002 * max(totalDrop, abs(initialRadius)));

    for candidateIdx = 2:nSamples-2
        collapseDrop = initialRadius - R_data(candidateIdx);
        if collapseDrop < minCollapseDrop || ...
                R_data(candidateIdx) > runningMin(candidateIdx) + noiseTol
            continue
        end

        nConfirm = min(5, nSamples - candidateIdx);
        reboundSegment = R_data(candidateIdx:candidateIdx + nConfirm);
        reboundRise = reboundSegment(end) - reboundSegment(1);
        minReboundRise = max(0.01 * max(abs(initialRadius), eps), ...
            0.05 * collapseDrop);
        nPositiveSteps = sum(diff(reboundSegment) > 0);
        sustainedRebound = min(reboundSegment(2:end)) >= ...
            reboundSegment(1) - noiseTol && ...
            nPositiveSteps >= ceil(0.6 * nConfirm) && ...
            reboundRise >= minReboundRise;

        if sustainedRebound
            [~, localOffset] = min(reboundSegment);
            firstCollapseIdx = candidateIdx + localOffset - 1;
            info.confirmationIdx = candidateIdx + nConfirm;
            info.reboundRise = R_data(info.confirmationIdx) - ...
                R_data(firstCollapseIdx);
            info.collapseDrop = initialRadius - R_data(firstCollapseIdx);
            return
        end
    end

    [~, firstCollapseIdx] = min(R_data);
    info.confirmationIdx = firstCollapseIdx;
    info.collapseDrop = initialRadius - R_data(firstCollapseIdx);
    info.usedFallback = true;
end

function paramSpec = buildParameterSpec(opt)
    paramSpec = repmat(struct('name', '', 'variableName', '', ...
        'bounds', [], 'optimizerBounds', [], 'logScale', false, ...
        'optimize', true, 'fixedValue', NaN), 5, 1);

    fixed = fixedParameterValues(opt);
    paramSpec(1) = makeParamSpec('G', opt.bounds.G, opt.logScale.G, ...
        fixed.G);
    paramSpec(2) = makeParamSpec('alph', opt.bounds.alph, ...
        opt.logScale.alph, fixed.alph);
    paramSpec(3) = makeParamSpec('mu', opt.bounds.mu, opt.logScale.mu, ...
        fixed.mu);
    paramSpec(4) = makeParamSpec('ani1', opt.bounds.ani(1, :), ...
        opt.logScale.ani(1), fixed.ani(1));
    paramSpec(5) = makeParamSpec('ani2', opt.bounds.ani(2, :), ...
        opt.logScale.ani(2), fixed.ani(2));
end

function fixed = fixedParameterValues(opt)
    fixed = struct('G', NaN, 'alph', NaN, 'mu', NaN, 'ani', [NaN, NaN]);
    if isfield(opt, 'fixed')
        fixed = copyFixedField(opt.fixed, fixed, 'G');
        fixed = copyFixedField(opt.fixed, fixed, 'alph');
        fixed = copyFixedField(opt.fixed, fixed, 'mu');
        if isfield(opt.fixed, 'ani') && ~isempty(opt.fixed.ani)
            fixed.ani = opt.fixed.ani(:).';
        end
    end
    if numel(fixed.ani) ~= 2
        error('opt.fixed.ani must contain exactly two values.');
    end
end

function fixed = copyFixedField(source, fixed, fieldName)
    if isfield(source, fieldName) && ~isempty(source.(fieldName))
        fixed.(fieldName) = source.(fieldName);
    end
end

function spec = makeParamSpec(name, bounds, logScale, fixedValue)
    if numel(bounds) ~= 2 || bounds(1) > bounds(2)
        error('Bounds for %s must be [lower upper] with lower <= upper.', ...
            name);
    end
    isFixed = isfinite(fixedValue) || bounds(1) == bounds(2);
    if bounds(1) == bounds(2)
        fixedValue = bounds(1);
    end
    if logScale && ~isFixed && any(bounds <= 0)
        error('Log-scaled parameter %s must have strictly positive bounds.', ...
            name);
    end
    spec.name = name;
    spec.bounds = bounds;
    spec.logScale = logScale;
    spec.optimize = ~isFixed;
    spec.fixedValue = fixedValue;
    if logScale
        spec.variableName = ['log10_' name];
        if isFixed
            spec.optimizerBounds = [];
        else
            spec.optimizerBounds = log10(bounds);
        end
    else
        spec.variableName = name;
        if isFixed
            spec.optimizerBounds = [];
        else
            spec.optimizerBounds = bounds;
        end
    end
end

function vars = buildBayesoptVariables(paramSpec)
    vars = [];
    for ii = 1:numel(paramSpec)
        if ~paramSpec(ii).optimize
            continue
        end
        thisVar = optimizableVariable(paramSpec(ii).variableName, ...
            paramSpec(ii).optimizerBounds);
        if isempty(vars)
            vars = thisVar;
        else
            vars(end+1, 1) = thisVar; %#ok<AGROW>
        end
    end
end

function [lbOpt, ubOpt] = optimizerBounds(paramSpec)
    activeSpec = activeParameterSpec(paramSpec);
    lbOpt = zeros(1, numel(activeSpec));
    ubOpt = zeros(1, numel(activeSpec));
    for ii = 1:numel(activeSpec)
        lbOpt(ii) = activeSpec(ii).optimizerBounds(1);
        ubOpt(ii) = activeSpec(ii).optimizerBounds(2);
    end
end

function initialX = makeLowerRegionLatinHypercubeInitialX(paramSpec, ...
    nPoints, regionFraction)
    paramSpec = activeParameterSpec(paramSpec);
    if nPoints < 1 || isempty(paramSpec)
        initialX = table();
        return
    end
    if regionFraction <= 0 || regionFraction > 1
        error('InitialRegionFraction must be in the interval (0, 1].');
    end

    nVars = numel(paramSpec);
    x = zeros(nPoints, nVars);
    for jj = 1:nVars
        lb = paramSpec(jj).optimizerBounds(1);
        ub = paramSpec(jj).optimizerBounds(2);
        upper = lb + regionFraction * (ub - lb);
        strata = ((0:nPoints-1)' + rand(nPoints, 1)) ./ nPoints;
        strata = strata(randperm(nPoints));
        x(:, jj) = lb + strata .* (upper - lb);
    end

    initialX = array2table(x, 'VariableNames', ...
        {paramSpec.variableName});
end

function bayesOpts = prepareBayesoptParallelPool(bayesOpts, scriptDir)
    if ~bayesOpts.UseParallel
        return
    end

    p = gcp('nocreate');
    if optionValue(bayesOpts, 'RestartParallelPool', false) && ~isempty(p)
        delete(p);
        p = [];
    end

    if isempty(p)
        try
            p = startRequestedPool(bayesOpts);
        catch ME
            warning('ModelToData:ParallelPoolStartFailed', ...
                ['Could not start a parallel pool for bayesopt: %s\n', ...
                'Continuing with UseParallel=false.'], ME.message);
            bayesOpts.UseParallel = false;
            return
        end
    end
    syncParallelWorkerPaths(p, scriptDir);

    nWorkers = p.NumWorkers;
    fprintf('Bayes-opt parallel pool: %d workers (%s).\n', ...
        nWorkers, class(p));
    if nWorkers < 2
        warning('Bayesopt has UseParallel=true but the pool has only one worker.');
    end

    if isempty(optionValue(bayesOpts, 'MinWorkerUtilization', []))
        bayesOpts.MinWorkerUtilization = nWorkers;
    elseif bayesOpts.MinWorkerUtilization > nWorkers
        warning(['MinWorkerUtilization (%d) exceeds pool size (%d); ', ...
            'using %d.'], bayesOpts.MinWorkerUtilization, nWorkers, ...
            nWorkers);
        bayesOpts.MinWorkerUtilization = nWorkers;
    end
    fprintf('Bayes-opt MinWorkerUtilization: %d.\n', ...
        bayesOpts.MinWorkerUtilization);

    if bayesOpts.NumSeedPoints < nWorkers
        fprintf(['Increasing NumSeedPoints from %d to %d so the initial ', ...
            'parallel batch can fill the pool.\n'], ...
            bayesOpts.NumSeedPoints, nWorkers);
        bayesOpts.NumSeedPoints = nWorkers;
    end
end

function syncParallelWorkerPaths(poolObj, scriptDir)
    workerPaths = {scriptDir, ...
        fullfile(scriptDir, 'src', 'common'), ...
        fullfile(scriptDir, 'src', 'forward_solver'), ...
        fullfile(scriptDir, 'src', 'characterization')};
    for ii = 1:numel(workerPaths)
        if isfolder(workerPaths{ii})
            addpath(workerPaths{ii});
        end
    end

    if isempty(poolObj)
        return
    end

    futures = cell(0, 1);
    for ii = 1:numel(workerPaths)
        if isfolder(workerPaths{ii})
            futures{end+1, 1} = parfevalOnAll(poolObj, @addpath, 0, ...
                workerPaths{ii}); %#ok<AGROW>
        end
    end
    for ii = 1:numel(futures)
        wait(futures{ii});
    end
    fprintf('Synchronized IMRv2 paths on parallel workers.\n');
end

function p = startRequestedPool(bayesOpts)
    profile = optionValue(bayesOpts, 'ParallelPoolProfile', '');
    nWorkers = optionValue(bayesOpts, 'NumWorkers', []);
    if isempty(profile) && isempty(nWorkers)
        p = parpool;
    elseif isempty(profile)
        p = parpool(nWorkers);
    elseif isempty(nWorkers)
        p = parpool(profile);
    else
        p = parpool(profile, nWorkers);
    end
end

function initializeDiagnosticLog(simOpts)
    logFile = optionValue(simOpts, 'DiagnosticLogFile', '');
    if isempty(logFile)
        return
    end
    appendLog = optionValue(simOpts, 'DiagnosticLogAppend', false);
    if appendLog
        mode = 'a';
    else
        mode = 'w';
    end
    try
        logDir = fileparts(logFile);
        if ~isempty(logDir) && ~isfolder(logDir)
            mkdir(logDir);
        end
        fid = fopen(logFile, mode);
        if fid < 0
            warning('ModelToData:DiagnosticLogOpenFailed', ...
                'Could not open diagnostic log: %s', logFile);
            return
        end
        cleanup = onCleanup(@() fclose(fid));
        if ~appendLog
            fprintf(fid, ['timestamp\tworker\tsuccess\tloss\trawLoss', ...
                '\telapsed\tcompleted_tstar\tdt_requested\tdt_radial', ...
                '\tdt_ep\ttimedOut\tidentifier\tmessage\tG\talph', ...
                '\tmu\tani1\tani2\n']);
        end
        clear cleanup
        fprintf('Objective diagnostic log: %s\n', logFile);
    catch ME
        warning('ModelToData:DiagnosticLogInitFailed', ...
            'Could not initialize diagnostic log %s: %s', ...
            logFile, ME.message);
    end
end

function bayesOpts = runOptimizerPreflight(bayesOpts, xData, paramSpec, simOpts)
    preflightPoint = makePreflightPoint(paramSpec, bayesOpts);
    clientLoss = NaN;
    workerLoss = NaN;

    if optionValue(bayesOpts, 'RunClientPreflight', false)
        fprintf('Running optimizer client preflight evaluation...\n');
        clientLoss = f_optimize_model_to_data_loss(preflightPoint, xData, ...
            paramSpec, simOpts);
        printPreflightLoss('client', clientLoss, simOpts);
    end

    if bayesOpts.UseParallel && optionValue(bayesOpts, ...
            'RunWorkerPreflight', false)
        poolObj = gcp('nocreate');
        if isempty(poolObj)
            warning('Worker preflight requested, but no parallel pool exists.');
        else
            fprintf('Running optimizer worker preflight evaluation...\n');
            try
                future = parfeval(poolObj, @f_optimize_model_to_data_loss, ...
                    1, preflightPoint, xData, paramSpec, simOpts);
                workerLoss = fetchOutputs(future);
                printPreflightLoss('worker', workerLoss, simOpts);
            catch ME
                warning('ModelToData:WorkerPreflightFailed', ...
                    'Worker preflight errored: %s', ME.message);
                if optionValue(bayesOpts, ...
                        'FallbackToSerialOnWorkerPreflightFailure', false)
                    bayesOpts.UseParallel = false;
                    fprintf(['Using serial bayesopt because the worker ', ...
                        'preflight errored.\n']);
                end
                return
            end
        end
    end

    if isfinite(clientLoss) && isfinite(workerLoss) && ...
            ~isPenaltyObjective(clientLoss, simOpts) && ...
            isPenaltyObjective(workerLoss, simOpts) && ...
            optionValue(bayesOpts, ...
            'FallbackToSerialOnWorkerPreflightFailure', false)
        bayesOpts.UseParallel = false;
        warning('ModelToData:WorkerPreflightPenalty', ...
            ['Client preflight succeeded but worker preflight returned ', ...
            'the failure penalty. Using serial bayesopt to avoid an ', ...
            'all-penalty parallel run. Check the diagnostic log for ', ...
            'the worker failure reason.']);
    end
end

function preflightPoint = makePreflightPoint(paramSpec, bayesOpts)
    activeSpec = activeParameterSpec(paramSpec);
    if isempty(activeSpec)
        preflightPoint = [];
        return
    end
    if isfield(bayesOpts, 'InitialX') && ~isempty(bayesOpts.InitialX)
        preflightPoint = bayesOpts.InitialX(1, :);
        return
    end

    x = zeros(1, numel(activeSpec));
    for ii = 1:numel(activeSpec)
        x(ii) = mean(activeSpec(ii).optimizerBounds);
    end
    preflightPoint = array2table(x, 'VariableNames', ...
        {activeSpec.variableName});
end

function printPreflightLoss(label, loss, simOpts)
    penaltyObjective = penaltyObjectiveValue(simOpts);
    if isPenaltyObjective(loss, simOpts)
        fprintf(['Optimizer %s preflight objective = %.6g ', ...
            '(failure penalty; diagnostic log has the reason)\n'], ...
            label, loss);
    else
        fprintf(['Optimizer %s preflight objective = %.6g ', ...
            '(penalty would be %.6g)\n'], label, loss, penaltyObjective);
    end
end

function tf = isPenaltyObjective(loss, simOpts)
    penaltyObjective = penaltyObjectiveValue(simOpts);
    tf = isfinite(loss) && abs(loss - penaltyObjective) <= ...
        1e-10 * max(1, abs(penaltyObjective));
end

function value = penaltyObjectiveValue(simOpts)
    rawPenalty = optionValue(simOpts, 'FailurePenalty', 1e6);
    if optionValue(simOpts, 'UseLogLoss', false)
        floorValue = optionValue(simOpts, 'LogLossFloor', 1e-12);
        value = log10(max(rawPenalty, floorValue));
    else
        value = rawPenalty;
    end
end

function value = optionValue(options, fieldName, defaultValue)
    if isfield(options, fieldName) && ~isempty(options.(fieldName))
        value = options.(fieldName);
    else
        value = defaultValue;
    end
end

function lossWeights = makePriorityLossWeights(n, lossOpts)
    radialWeight = optionValue(lossOpts, 'RadialWeight', 1);
    modeBaseWeight = optionValue(lossOpts, 'ModeBaseWeight', 1);
    modeReference = optionValue(lossOpts, 'ModeReference', 2);
    modePower = optionValue(lossOpts, 'ModeWeightPower', 2);
    modeOffset = optionValue(lossOpts, 'ModeWeightOffset', 0);

    modeDenom = abs(n(:).') + modeOffset;
    if modeReference <= 0 || any(modeDenom <= 0)
        error(['Mode weights require ModeReference and n + ' ...
            'ModeWeightOffset to be positive.']);
    end

    lossWeights = struct();
    lossWeights.radial = radialWeight;
    lossWeights.modes = modeBaseWeight .* ...
        (modeReference ./ modeDenom) .^ modePower;
    lossWeights.modeBaseWeight = modeBaseWeight;
    lossWeights.modeReference = modeReference;
    lossWeights.modePower = modePower;
    lossWeights.modeOffset = modeOffset;
end

function printLossWeights(n, fitModeIdx, lossWeights)
    fprintf(['Loss priority weights: radial=%.6g, mode formula=' ...
        '%.6g*(%.6g/(n+%.6g))^%.6g\n'], lossWeights.radial, ...
        lossWeights.modeBaseWeight, lossWeights.modeReference, ...
        lossWeights.modeOffset, lossWeights.modePower);
    fprintf('Training mode priority weights:\n');
    for ii = 1:numel(fitModeIdx)
        idx = fitModeIdx(ii);
        fprintf('  n=%g: %.6g\n', n(idx), lossWeights.modes(idx));
    end
end

function args = makeBayesoptArgs(bayesOpts)
    args = { ...
        'MaxObjectiveEvaluations', bayesOpts.MaxObjectiveEvaluations, ...
        'UseParallel', bayesOpts.UseParallel, ...
        'IsObjectiveDeterministic', bayesOpts.IsObjectiveDeterministic, ...
        'AcquisitionFunctionName', bayesOpts.AcquisitionFunctionName, ...
        'NumSeedPoints', bayesOpts.NumSeedPoints, ...
        'ExplorationRatio', bayesOpts.ExplorationRatio, ...
        'Verbose', bayesOpts.Verbose};
    if bayesOpts.UseParallel && isfield(bayesOpts, 'MinWorkerUtilization') ...
            && ~isempty(bayesOpts.MinWorkerUtilization)
        args = [args, {'MinWorkerUtilization', ...
            bayesOpts.MinWorkerUtilization}];
    end
    if bayesOpts.UseParallel && isfield(bayesOpts, 'ParallelMethod') && ...
            ~isempty(bayesOpts.ParallelMethod)
        args = [args, {'ParallelMethod', bayesOpts.ParallelMethod}];
    end
    if isfield(bayesOpts, 'InitialX') && ~isempty(bayesOpts.InitialX)
        args = [args, {'InitialX', bayesOpts.InitialX}];
    end
    if ~isempty(bayesOpts.PlotFcn)
        args = [args, {'PlotFcn', bayesOpts.PlotFcn}];
    end
end

function lsqopts = makeLsqOptions(refineOpts)
    lsqopts = optimoptions('lsqcurvefit', ...
        'Display', refineOpts.Display, ...
        'OptimalityTolerance', refineOpts.OptimalityTolerance, ...
        'FunctionTolerance', refineOpts.FunctionTolerance, ...
        'StepTolerance', refineOpts.StepTolerance, ...
        'MaxFunctionEvaluations', refineOpts.MaxFunctionEvaluations, ...
        'MaxIterations', refineOpts.MaxIterations, ...
        'UseParallel', false, ...
        'FiniteDifferenceType', refineOpts.FiniteDifferenceType, ...
        'Algorithm', refineOpts.Algorithm);
end

function starts = bestBayesStarts(results, paramSpec, nStarts)
    objectiveTrace = results.ObjectiveTrace;
    xTrace = tableToOptimizerMatrix(results.XTrace, paramSpec);
    [~, order] = sort(objectiveTrace, 'ascend');
    nKeep = min(nStarts, numel(order));
    starts = xTrace(order(1:nKeep), :);
end

function mat = tableToOptimizerMatrix(T, paramSpec)
    activeSpec = activeParameterSpec(paramSpec);
    mat = zeros(height(T), numel(activeSpec));
    for ii = 1:numel(activeSpec)
        mat(:, ii) = T.(activeSpec(ii).variableName);
    end
end

function activeSpec = activeParameterSpec(paramSpec)
    activeSpec = paramSpec([paramSpec.optimize]);
end

function printParameterSpec(paramSpec)
    fprintf('Optimized parameters: ');
    activeSpec = activeParameterSpec(paramSpec);
    if isempty(activeSpec)
        fprintf('(none; evaluating fixed parameter set only)\n');
    else
        fprintf('%s\n', strjoin({activeSpec.name}, ', '));
    end

    fixedSpec = paramSpec(~[paramSpec.optimize]);
    if isempty(fixedSpec)
        fprintf('Fixed parameters: (none)\n');
        return
    end
    fixedText = cell(1, numel(fixedSpec));
    for ii = 1:numel(fixedSpec)
        fixedText{ii} = sprintf('%s=%.6g', fixedSpec(ii).name, ...
            fixedSpec(ii).fixedValue);
    end
    fprintf('Fixed parameters: %s\n', strjoin(fixedText, ', '));
end

function loss = computePerturbationSubsetLoss(sim, xData, modeIdx)
    if isempty(modeIdx)
        loss = NaN;
        return
    end
    if ~isfield(sim, 'success') || ~sim.success
        loss = NaN;
        return
    end

    epRows = perturbationFitRows(xData);
    if numel(sim.t) ~= numel(xData.tfit_nd)
        loss = NaN;
        return
    end

    timeTol = 1e-8;
    if max(abs(sim.t(:) - xData.tfit_nd(:))) > timeTol
        loss = NaN;
        return
    end

    epModel = sim.epnm(epRows, modeIdx);
    epData = xData.ep_data(epRows, modeIdx);
    weights = xData.aEP(modeIdx);

    if any(~isfinite(epModel(:)))
        loss = NaN;
        return
    end

    weightedData = epData .* weights;
    weightedModel = epModel .* weights;
    denom = norm(weightedData(:));
    if denom < eps
        loss = NaN;
    else
        loss = sqrt(sum((weightedData(:) - weightedModel(:)).^2)) / denom;
    end
end

function epRows = perturbationFitRows(xData)
    if isfield(xData, 'epFitIdx') && ~isempty(xData.epFitIdx)
        epRows = xData.epFitIdx(:);
    else
        epRows = (1:numel(xData.tfit_nd)).';
    end
end

function verification = makeTimeExtractionVerification(sim, xData, runInfo)
    epRows = perturbationFitRows(xData);
    radialRequested = xData.tfit_nd(:);
    perturbationRequested = xData.tfit_nd(epRows);

    if isfield(sim, 't')
        returnedTimes = sim.t(:);
    else
        returnedTimes = [];
    end

    [radialReturned, radialAbsDt] = nearestReturnedTimes(returnedTimes, ...
        radialRequested);
    [perturbationReturned, perturbationAbsDt] = nearestReturnedTimes( ...
        returnedTimes, perturbationRequested);

    verification = struct();
    verification.radialTable = table(radialRequested, radialReturned, ...
        radialAbsDt, 'VariableNames', {'requestedTimeNd', ...
        'nearestReturnedTimeNd', 'absDtNd'});
    verification.perturbationTable = table(epRows, perturbationRequested, ...
        perturbationReturned, perturbationAbsDt, 'VariableNames', ...
        {'rowIndex', 'requestedTimeNd', 'nearestReturnedTimeNd', ...
        'absDtNd'});
    verification.maxRadialAbsDt = maxFinite(radialAbsDt);
    verification.maxPerturbationAbsDt = maxFinite(perturbationAbsDt);
    verification.firstPerturbationRow = epRows(1);
    verification.lastPerturbationRow = epRows(end);
    verification.firstCollapseIdx = xData.firstCollapseIdx;
    verification.firstCollapseTimeNd = xData.firstCollapseTimeNd;
    verification.optimizerPerturbationTimeNd = perturbationRequested;
    verification.runInfo = runInfo;
end

function [nearestTimes, absDt] = nearestReturnedTimes(returnedTimes, ...
    requestedTimes)
    returnedTimes = returnedTimes(:);
    requestedTimes = requestedTimes(:);
    nearestTimes = NaN(size(requestedTimes));
    absDt = NaN(size(requestedTimes));
    if isempty(returnedTimes) || isempty(requestedTimes)
        return
    end

    for ii = 1:numel(requestedTimes)
        [absDt(ii), nearestIdx] = min(abs(returnedTimes - ...
            requestedTimes(ii)));
        nearestTimes(ii) = returnedTimes(nearestIdx);
    end
end

function value = maxFinite(values)
    finiteValues = values(isfinite(values));
    if isempty(finiteValues)
        value = NaN;
    else
        value = max(finiteValues);
    end
end

function plotOptimizedFit(bestSim, xData, bestParams)
    if ~isfield(bestSim, 'success') || ~bestSim.success
        warning('Best simulation did not complete, so no fit plot was made.');
        return
    end

    figure('Name', 'Optimized model fit')
    plotl = ceil(sqrt(numel(xData.n) + 1));

    subplot(plotl, plotl, 1)
    hold on
    plot(bestSim.t, bestSim.R, '-', 'LineWidth', 1.5)
    plot(xData.tfit_nd, xData.R_data, 'o')
    plot(xData.firstCollapseTimeNd, ...
        xData.R_data(xData.firstCollapseIdx), 'kp', ...
        'MarkerFaceColor', 'y', 'MarkerSize', 10)
    xlabel("t^*")
    ylabel("R")
    title(sprintf('G=%.3g, \\mu=%.3g, \\alpha=%.3g, ani=[%.3g %.3g]', ...
        bestParams.G, bestParams.mu, bestParams.alph, bestParams.ani(1), ...
        bestParams.ani(2)))

    for ii = 1:numel(xData.n)
        subplot(plotl, plotl, ii + 1)
        hold on
        plot(bestSim.t, bestSim.epnm(:, ii), '-', 'LineWidth', 1.5)
        plot(xData.tfit_nd, xData.ep_data(:, ii), 'r^')
        xlabel("t^*")
        ylabel(sprintf('$\\epsilon_{%.0f}$', xData.n(ii)), ...
            'Interpreter', 'latex')
    end
end

function plotCollapseDetection(tfit_nd, R_data, collapseInfo)
    figure('Name', 'First-collapse detection')
    plot(tfit_nd, R_data, 'o-', 'DisplayName', 'experimental radius')
    hold on
    collapseIdx = find(R_data == min(R_data(1: ...
        collapseInfo.confirmationIdx)), 1, 'first');
    plot(tfit_nd(collapseIdx), R_data(collapseIdx), 'kp', ...
        'MarkerFaceColor', 'y', 'MarkerSize', 11, ...
        'DisplayName', 'first collapse')
    plot(tfit_nd(collapseInfo.confirmationIdx), ...
        R_data(collapseInfo.confirmationIdx), 'ks', ...
        'MarkerFaceColor', 'c', 'MarkerSize', 8, ...
        'DisplayName', 'rebound confirmation')
    xlabel("t^*")
    ylabel("R/R_{max}")
    title('First collapse from raw radius and sustained rebound')
    legend('Location', 'best')
    grid on
end

function plotInitialConditionCheck(texp, maxidx, amps_og, amps, epnm0, epnmd0, tc)
    figure('Name', 'Initial-condition check')
    nmodes = size(amps, 1);
    plotl = ceil(sqrt(nmodes));
    t0 = texp(maxidx);
    for ii = 1:nmodes
        subplot(plotl, plotl, ii)
        plot((texp - t0)./tc, amps_og(ii, :), 'o')
        hold on
        plot((texp - t0)./tc, amps(ii, :), '-')
        plot((texp - t0)./tc, 0 .* (texp - t0) + epnm0(ii), ':')
        plot((texp - t0)./tc, ((texp - t0) ./ tc) .* epnmd0(ii) + epnm0(ii), '-')
        % xlim([min(texp - t0), -min(texp - t0)])
    end
end
