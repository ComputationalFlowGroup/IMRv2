%% Evaluate one anisotropic IMR model against experimental data
clear
clc
close all

[scriptDir, projectDir] = locateProjectPaths();
addpath(fullfile(scriptDir, 'src', 'common'));
addpath(fullfile(scriptDir, 'src', 'forward_solver'));
addpath(fullfile(scriptDir, 'src', 'characterization'));

%% User settings
dataFile = fullfile(projectDir, 'data', 'SicongJinChicken', ...
    'PVA_Rt_data/', 'processed_14.mat');

material = "PVA";
maxmode = 22;
polyOrder = 3;
windowPts = 9;   % odd integer: 3, 5, 7, ...
icVelocityWindowPts = 9;  % forward polynomial derivative window from Rmax
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

% Physical model parameters for this single simulation.
modelParams.G = 10^(5.0);       % Pa
modelParams.alph = 10^(-1);
modelParams.mu = 10^(-0.8);        % Pa*s
modelParams.ani = [1, 2];

% Forward-solver controls. The simulation is evaluated at the experimental
% post-Rmax times, so tsteps is only a fallback for non-optimization calls.
opt.sim.tsteps = 3000;
opt.sim.RelTol = 1e-4;
opt.sim.AbsTol = 1e-5;
opt.sim.Nt = 75;
opt.sim.Method = 23;
opt.sim.MaxWallTime = 45;
opt.sim.UseHardTimeout = false;
opt.sim.TimeoutPollInterval = 1;
opt.sim.FailurePenalty = 1e5;
opt.sim.UseLogLoss = true;
opt.sim.LogLossFloor = 1e-12;
opt.sim.TimeMatchTolerance = 1e-8;
opt.sim.VerifyTimeExtraction = true;
opt.sim.PrintFailures = true;
opt.sim.PrintSuccess = true;

opt.outputFile = fullfile(scriptDir, 'single_model_to_data_eval.mat');
opt.makePlot = true;

% Load and process data
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

fitIdx = maxidx:numel(texp);
tfit_nd = (texp(fitIdx) - texp(maxidx)) ./ tc;
R_data = expR(fitIdx).' ./ Rmax;
ep_data = amps(:, fitIdx).';
tf_nd = max(tfit_nd);
[firstCollapseIdx, collapseInfo] = findFirstCollapseIndex(R_data);
epFitIdx = 1:firstCollapseIdx;
firstCollapseTimeNd = tfit_nd(firstCollapseIdx);
firstCollapseTimeSeconds = firstCollapseTimeNd * tc;

modeRows = 3:maxmode+1;
n = mode_extract_fft(modeRows, 10);
n = n(:).';
m = zeros(size(n));

nmodes = size(ep_data, 2);
modeEnergy = sum(ep_data(epFitIdx, :).^2, 1);
[~, modeEnergyOrder] = sort(modeEnergy, 'descend');
nFitModes = min(opt.fit.NumPerturbationModes, nmodes);
fitModeIdx = sort(modeEnergyOrder(1:nFitModes));
testModeIdx = setdiff(1:nmodes, fitModeIdx, 'stable');

lossWeights = makePriorityLossWeights(n, opt.loss);
aR = lossWeights.radial / norm(R_data);
sEP = vecnorm(ep_data(epFitIdx, :), 2, 1);
sEP(sEP < eps) = 1;
aEP = lossWeights.modes ./ sEP;
yData = [aR .* R_data; reshape(ep_data(epFitIdx, fitModeIdx) .* ...
    aEP(fitModeIdx), [], 1)];

xData = struct( ...
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
    'y_data', yData, ...
    'ultra', false);

fprintf('Single-simulation model-to-data evaluation\n');
fprintf('  G     = %.6g Pa\n', modelParams.G);
fprintf('  alph  = %.6g\n', modelParams.alph);
fprintf('  mu    = %.6g Pa*s\n', modelParams.mu);
fprintf('  ani   = [%.6g %.6g]\n', modelParams.ani(1), ...
    modelParams.ani(2));
fprintf('Training perturbation modes: %s\n', num2str(n(fitModeIdx)));
fprintf('Held-out perturbation modes: %s\n', num2str(n(testModeIdx)));
printLossWeights(n, fitModeIdx, lossWeights);
fprintf(['First collapse (raw radius minimum) at post-Rmax sample ', ...
    '%d/%d: t* = %.6g, dt = %.6g us, R/Rmax = %.6g\n'], ...
    firstCollapseIdx, numel(tfit_nd), firstCollapseTimeNd, ...
    firstCollapseTimeSeconds * 1e6, R_data(firstCollapseIdx));
fprintf(['  Sustained rebound confirmed through sample %d: t* = %.6g, ', ...
    'rise in R/Rmax = %.6g\n'], collapseInfo.confirmationIdx, ...
    tfit_nd(collapseInfo.confirmationIdx), collapseInfo.reboundRise);
fprintf('Perturbation loss uses %d samples from t* = %.6g to %.6g.\n', ...
    numel(epFitIdx), tfit_nd(epFitIdx(1)), tfit_nd(epFitIdx(end)));

% Run one simulation and compute losses
paramSpec = fixedParamSpecFromParams(modelParams);
[yModel, runInfo, sim] = f_optimize_model_to_data_predict([], xData, ...
    paramSpec, opt.sim);
isotropicParams = modelParams;
isotropicParams.ani = [0, 0];
[~, isotropicRunInfo, isotropicSim] = f_optimize_model_to_data_predict([], ...
    xData, fixedParamSpecFromParams(isotropicParams), opt.sim);

if runInfo.success
    rawObjectiveLoss = sqrt(sum((yData - yModel).^2)) / norm(yData);
    optimizerLoss = transformLossForOptimizer(rawObjectiveLoss, opt.sim);
    fitR2 = 1 - sum((yModel - yData).^2) / ...
        sum((yData - mean(yData)).^2);
else
    rawObjectiveLoss = opt.sim.FailurePenalty;
    optimizerLoss = transformLossForOptimizer(rawObjectiveLoss, opt.sim);
    fitR2 = NaN;
end

radialLoss = computeRadialLoss(sim, xData);
trainModeLoss = computePerturbationSubsetLoss(sim, xData, ...
    xData.fitModeIdx);
heldoutModeLoss = computePerturbationSubsetLoss(sim, xData, ...
    xData.testModeIdx);
allModeLoss = computePerturbationSubsetLoss(sim, xData, 1:nmodes);
timeVerification = makeTimeExtractionVerification(sim, xData, runInfo);

fprintf('\nLoss summary\n');
fprintf('  run success             = %d\n', runInfo.success);
fprintf('  run message             = %s\n', runInfo.message);
fprintf('  elapsed                 = %.3f s\n', runInfo.elapsed);
fprintf('  raw objective loss      = %.6g\n', rawObjectiveLoss);
fprintf('  optimizer loss          = %.6g\n', optimizerLoss);
fprintf('  R2                      = %.6g\n', fitR2);
fprintf('  radial loss             = %.6g\n', radialLoss);
fprintf('  train mode loss         = %.6g\n', trainModeLoss);
fprintf('  held-out mode loss      = %.6g\n', heldoutModeLoss);
fprintf('  all mode loss           = %.6g\n', allModeLoss);
fprintf('  max returned time error = %.3g\n', ...
    timeVerification.maxRadialAbsDt);
fprintf('  isotropic comparison success = %d (%s)\n', ...
    isotropicRunInfo.success, isotropicRunInfo.message);

if opt.makePlot
    plotSingleEvaluation(sim, isotropicSim, xData, modelParams);
end

save(opt.outputFile, 'opt', 'modelParams', 'paramSpec', 'xData', ...
    'yData', 'yModel', 'sim', 'runInfo', 'isotropicParams', ...
    'isotropicRunInfo', 'isotropicSim', 'timeVerification', ...
    'rawObjectiveLoss', 'optimizerLoss', 'fitR2', 'radialLoss', ...
    'trainModeLoss', 'heldoutModeLoss', 'allModeLoss');

%% Local helper functions
function [scriptDir, projectDir] = locateProjectPaths()
    candidates = {};
    candidates = addCandidate(candidates, fileparts(mfilename('fullpath')));
    candidates = addCandidate(candidates, fileparts(which( ...
        's_evaluate_model_to_data_loss')));
    candidates = addCandidate(candidates, fileparts(which( ...
        'f_call_IMRv2_exp')));
    candidates = addCandidate(candidates, pwd);
    candidates = addCandidate(candidates, fullfile(pwd, 'IMRv2'));
    candidates = addCandidate(candidates, fileparts(pwd));
    candidates = addCandidate(candidates, fullfile(fileparts(pwd), ...
        'IMRv2'));

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

function paramSpec = fixedParamSpecFromParams(modelParams)
    paramSpec = repmat(struct('name', '', 'variableName', '', ...
        'bounds', [], 'optimizerBounds', [], 'logScale', false, ...
        'optimize', false, 'fixedValue', NaN), 5, 1);

    paramSpec(1) = fixedParamSpec('G', 'log10_G', modelParams.G, true);
    paramSpec(2) = fixedParamSpec('alph', 'alph', modelParams.alph, false);
    paramSpec(3) = fixedParamSpec('mu', 'log10_mu', modelParams.mu, true);
    paramSpec(4) = fixedParamSpec('ani1', 'ani1', modelParams.ani(1), false);
    paramSpec(5) = fixedParamSpec('ani2', 'ani2', modelParams.ani(2), false);
end

function spec = fixedParamSpec(name, variableName, fixedValue, logScale)
    spec = struct('name', name, 'variableName', variableName, ...
        'bounds', [fixedValue, fixedValue], 'optimizerBounds', [], ...
        'logScale', logScale, 'optimize', false, ...
        'fixedValue', fixedValue);
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

function value = optionValue(options, fieldName, defaultValue)
    if isfield(options, fieldName) && ~isempty(options.(fieldName))
        value = options.(fieldName);
    else
        value = defaultValue;
    end
end

function loss = transformLossForOptimizer(rawLoss, simOpts)
    if simOpt(simOpts, 'UseLogLoss', false)
        floorValue = simOpt(simOpts, 'LogLossFloor', 1e-12);
        loss = log10(max(rawLoss, floorValue));
    else
        loss = rawLoss;
    end
end

function loss = computeRadialLoss(sim, xData)
    if ~isfield(sim, 'success') || ~sim.success
        loss = NaN;
        return
    end
    if ~timesMatch(sim, xData)
        loss = NaN;
        return
    end
    denom = norm(xData.R_data(:));
    if denom < eps
        loss = NaN;
    else
        loss = sqrt(sum((xData.R_data(:) - sim.R(:)).^2)) / denom;
    end
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
    if ~timesMatch(sim, xData)
        loss = NaN;
        return
    end

    epRows = perturbationFitRows(xData);
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

function tf = timesMatch(sim, xData)
    timeTol = 1e-8;
    tf = numel(sim.t) == numel(xData.tfit_nd) && ...
        max(abs(sim.t(:) - xData.tfit_nd(:))) <= timeTol;
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

function value = simOpt(simOpts, fieldName, defaultValue)
    if isfield(simOpts, fieldName) && ~isempty(simOpts.(fieldName))
        value = simOpts.(fieldName);
    else
        value = defaultValue;
    end
end

function plotSingleEvaluation(sim, isotropicSim, xData, modelParams)
    if ~isfield(sim, 'success') || ~sim.success
        warning('Simulation did not complete, so no fit plot was made.');
        return
    end
    hasIsotropic = isfield(isotropicSim, 'success') && isotropicSim.success;

    figure('Name', 'Single model-to-data evaluation')
    plotl = ceil(sqrt(numel(xData.n) + 1));

    subplot(plotl, plotl, 1)
    hold on
    hFit = plot(sim.t, sim.R, '-', 'LineWidth', 1.5);
    if hasIsotropic
        hIso = plot(isotropicSim.t, isotropicSim.R, 'k--', ...
            'LineWidth', 1.2);
    else
        hIso = gobjects(0);
    end
    hData = plot(xData.tfit_nd, xData.R_data, 'o');
    hCollapse = plot(xData.firstCollapseTimeNd, ...
        xData.R_data(xData.firstCollapseIdx), 'kp', ...
        'MarkerFaceColor', 'y', 'MarkerSize', 10);
    xlabel("t^*")
    ylabel("R")
    title(sprintf('G=%.3g, \\mu=%.3g, \\alpha=%.3g, ani=[%.3g %.3g]', ...
        modelParams.G, modelParams.mu, modelParams.alph, ...
        modelParams.ani(1), modelParams.ani(2)))
    if hasIsotropic
        legend([hFit, hIso, hData, hCollapse], ...
            {'model', 'ani=[0 0]', 'data', 'first collapse'}, ...
            'Location', 'best')
    else
        legend([hFit, hData, hCollapse], ...
            {'model', 'data', 'first collapse'}, 'Location', 'best')
    end

    for ii = 1:numel(xData.n)
        subplot(plotl, plotl, ii + 1)
        hold on
        plot(sim.t, sim.epnm(:, ii), '-', 'LineWidth', 1.5)
        if hasIsotropic
            plot(isotropicSim.t, isotropicSim.epnm(:, ii), 'k--', ...
                'LineWidth', 1.2)
        end
        plot(xData.tfit_nd, xData.ep_data(:, ii), 'r^')
        xlabel("t^*")
        ylabel(sprintf('$\\epsilon_{%.0f}$', xData.n(ii)), ...
            'Interpreter', 'latex')
    end
end
