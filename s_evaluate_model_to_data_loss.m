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
    'chicken_Rt_data/', 'processed_11.mat');

material = "PVA";
maxmode = 22;
polyOrder = 3;
windowPts = 8;   % adjusted to the nearest valid odd window below this value
icVelocityWindowPts = 5;  % forward polynomial derivative window from Rmax
icVelocityPolyOrder = 3;
opt.fit.NumPerturbationModes = 5;
opt.fit.NumRadialCollapses = 4;
opt.fit.CollapseMinProminenceFraction = 0.08;
opt.fit.CollapseMinSeparationSteps = 5;
opt.fit.StopAtRadiusEquilibrium = true;
opt.fit.RadiusEquilibriumRelTol = 0.02;
opt.fit.RadiusEquilibriumAbsTol = 0;
opt.fit.RadiusEquilibriumConsecutiveSteps = 10;

% Loss priority weights after per-trace normalization. Radial receives
% RadialWeight, and mode n receives
% ModeBaseWeight*(ModeReference/(n + ModeWeightOffset))^ModeWeightPower.
opt.loss.RadialWeight = 1;
opt.loss.ModeBaseWeight = 1;
opt.loss.ModeReference = 2;
opt.loss.ModeWeightPower = 1;
opt.loss.ModeWeightOffset = 0;

% Physical model parameters for this single simulation.
modelParams.G = 10^(5.25);       % Pa
modelParams.alph = 0;
modelParams.mu = 10^(-0.1);        % Pa*s
modelParams.ani = [0, 0];

% Forward-solver controls. The simulation is evaluated at the experimental
% post-Rmax times, so tsteps is only a fallback for non-optimization calls.
opt.sim.tsteps = 3000;
opt.sim.RelTol = 1e-4;
opt.sim.AbsTol = 1e-5;
opt.sim.Nt = 75;
opt.sim.Method = 23;
opt.sim.MaxWallTime = 120;
opt.sim.UseHardTimeout = false;
opt.sim.TimeoutPollInterval = 1;
opt.sim.FailurePenalty = 1e5;
opt.sim.UseLogLoss = true;
opt.sim.LogLossFloor = 1e-12;
opt.sim.TimeMatchTolerance = 1e-8;
opt.sim.AllowNearestTimeExtraction = true;
opt.sim.NearestTimeTolerance = 1e-6;
opt.sim.VerifyTimeExtraction = true;
opt.sim.PrintFailures = true;
opt.sim.PrintSuccess = true;

opt.outputFile = fullfile(scriptDir, 'single_model_to_data_eval.mat');
opt.makePlot = true;
opt.plotInitialConditionCheck = true;

% Load and process data
if ~isfile(dataFile)
    error('Could not find the data file: %s', dataFile);
end
loadedData = load(dataFile, 'amp_extract_fft', 'mode_extract_fft');
requiredDataVariables = {'amp_extract_fft', 'mode_extract_fft'};
for variableIdx = 1:numel(requiredDataVariables)
    variableName = requiredDataVariables{variableIdx};
    if ~isfield(loadedData, variableName)
        error('Data file does not contain %s: %s', variableName, dataFile);
    end
end
amp_extract_fft = loadedData.amp_extract_fft;
mode_extract_fft = loadedData.mode_extract_fft;
if size(mode_extract_fft, 2) ~= size(amp_extract_fft, 2)
    error(['mode_extract_fft and amp_extract_fft must contain the same ', ...
        'number of experimental time samples in %s.'], dataFile);
end

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

[collapseFitEndIdx, radialCollapseCutoffInfo] = ...
    findRadialCollapseCutoff(expR, maxidx, opt.fit);
if radialCollapseCutoffInfo.foundRequestedCount
    fitEndIdx = collapseFitEndIdx;
    [~, radiusEquilibriumInfo] = findRadiusEquilibriumCutoff( ...
        expR, Req, maxidx, opt.fit);
    radiusEquilibriumInfo.usedAsFallback = false;
else
    [fitEndIdx, radiusEquilibriumInfo] = findRadiusEquilibriumCutoff( ...
        expR, Req, maxidx, opt.fit);
    radiusEquilibriumInfo.usedAsFallback = true;
end
fitIdx = maxidx:fitEndIdx;
tfit_nd = (texp(fitIdx) - texp(maxidx)) ./ tc;
R_data = expR(fitIdx).' ./ Rmax;
ep_data = amps(:, fitIdx).';
tf_nd = max(tfit_nd);
[firstCollapseIdx, collapseInfo] = findFirstCollapseIndex(R_data);
epFitIdx = (1:numel(tfit_nd)).';
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
    'radialCollapseCutoffInfo', radialCollapseCutoffInfo, ...
    'radiusEquilibriumInfo', radiusEquilibriumInfo, ...
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
printFitWindowCutoff(radialCollapseCutoffInfo, radiusEquilibriumInfo, ...
    texp, maxidx);
fprintf('Perturbation loss uses %d samples from t* = %.6g to %.6g.\n', ...
    numel(epFitIdx), tfit_nd(epFitIdx(1)), tfit_nd(epFitIdx(end)));

if opt.plotInitialConditionCheck
    fitEndPostRmaxIdx = fitEndIdx - maxidx + 1;
    tpost_nd = (texp(maxidx:end) - texp(maxidx)) ./ tc;
    Rpost_data = expR(maxidx:end).' ./ Rmax;
    plotInitialConditionCheck(texp, maxidx, fitEndIdx, amps_og, amps, ...
        epnm0, epnmd0, tc);
    plotCollapseDetection(tpost_nd, Rpost_data, collapseInfo, ...
        fitEndPostRmaxIdx);
end

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

function [fitEndIdx, info] = findRadialCollapseCutoff( ...
        radius, startIdx, fitOpts)
    radius = radius(:);
    nSamples = numel(radius);
    fitEndIdx = nSamples;

    requestedCount = optionValue(fitOpts, 'NumRadialCollapses', 4);
    prominenceFraction = optionValue(fitOpts, ...
        'CollapseMinProminenceFraction', 0.08);
    minSeparation = optionValue(fitOpts, ...
        'CollapseMinSeparationSteps', 5);

    if ~isscalar(requestedCount) || ~isfinite(requestedCount) || ...
            requestedCount < 1 || requestedCount ~= round(requestedCount)
        error('opt.fit.NumRadialCollapses must be a positive integer.');
    end
    if ~isscalar(prominenceFraction) || ~isfinite(prominenceFraction) || ...
            prominenceFraction < 0
        error(['opt.fit.CollapseMinProminenceFraction must be finite ', ...
            'and nonnegative.']);
    end
    if ~isscalar(minSeparation) || ~isfinite(minSeparation) || ...
            minSeparation < 1 || minSeparation ~= round(minSeparation)
        error(['opt.fit.CollapseMinSeparationSteps must be a positive ', ...
            'integer.']);
    end
    if any(~isfinite(radius))
        error('Finite experimental radius data are required.');
    end

    startIdx = min(max(1, round(startIdx)), nSamples);
    postMaxRadius = radius(startIdx:end);
    radialRange = max(postMaxRadius) - min(postMaxRadius);
    minProminence = prominenceFraction * radialRange;
    [~, localCollapseIdx, ~, prominences] = findpeaks( ...
        -postMaxRadius, 'MinPeakProminence', minProminence, ...
        'MinPeakDistance', minSeparation);
    collapseIdx = startIdx + localCollapseIdx - 1;

    info = struct('requestedCount', requestedCount, ...
        'detectedCount', numel(collapseIdx), ...
        'foundRequestedCount', numel(collapseIdx) >= requestedCount, ...
        'collapseIdx', collapseIdx(:), 'prominences', prominences(:), ...
        'prominenceFraction', prominenceFraction, ...
        'minProminence', minProminence, ...
        'minSeparationSteps', minSeparation, 'fitEndIdx', nSamples);
    if info.foundRequestedCount
        fitEndIdx = collapseIdx(requestedCount);
        info.fitEndIdx = fitEndIdx;
    end
end

function [fitEndIdx, info] = findRadiusEquilibriumCutoff( ...
        radius, Req, startIdx, fitOpts)
    radius = radius(:);
    nSamples = numel(radius);
    fitEndIdx = nSamples;

    enabled = optionValue(fitOpts, 'StopAtRadiusEquilibrium', false);
    relTol = optionValue(fitOpts, 'RadiusEquilibriumRelTol', 0.02);
    absTol = optionValue(fitOpts, 'RadiusEquilibriumAbsTol', 0);
    nConsecutive = optionValue(fitOpts, ...
        'RadiusEquilibriumConsecutiveSteps', 10);

    if ~isscalar(relTol) || ~isfinite(relTol) || relTol < 0
        error('opt.fit.RadiusEquilibriumRelTol must be finite and nonnegative.');
    end
    if ~isscalar(absTol) || ~isfinite(absTol) || absTol < 0
        error('opt.fit.RadiusEquilibriumAbsTol must be finite and nonnegative.');
    end
    if ~isscalar(nConsecutive) || ~isfinite(nConsecutive) || ...
            nConsecutive < 1 || nConsecutive ~= round(nConsecutive)
        error(['opt.fit.RadiusEquilibriumConsecutiveSteps must be a ', ...
            'positive integer.']);
    end
    if any(~isfinite(radius)) || ~isscalar(Req) || ~isfinite(Req)
        error('Finite experimental radius data and Req are required.');
    end

    startIdx = min(max(1, round(startIdx)), nSamples);
    tolerance = max(absTol, relTol * abs(Req));
    info = struct('enabled', logical(enabled), 'found', false, ...
        'startIdx', NaN, 'confirmationIdx', NaN, ...
        'fitEndIdx', fitEndIdx, 'nConsecutive', nConsecutive, ...
        'relativeTolerance', relTol, 'absoluteTolerance', absTol, ...
        'radiusTolerance', tolerance, 'Req', Req, ...
        'nExcludedSamples', 0);
    if ~enabled
        return
    end

    nearEquilibrium = abs(radius(startIdx:end) - Req) <= tolerance;
    runLength = 0;
    for ii = 1:numel(nearEquilibrium)
        if nearEquilibrium(ii)
            runLength = runLength + 1;
        else
            runLength = 0;
        end
        if runLength >= nConsecutive
            info.startIdx = startIdx + ii - nConsecutive;
            info.confirmationIdx = startIdx + ii - 1;
            fitEndIdx = info.confirmationIdx;
            info.fitEndIdx = fitEndIdx;
            info.found = true;
            info.nExcludedSamples = nSamples - fitEndIdx;
            return
        end
    end
end

function printFitWindowCutoff(collapseInfo, equilibriumInfo, texp, maxidx)
    if collapseInfo.foundRequestedCount
        selectedIdx = collapseInfo.collapseIdx( ...
            collapseInfo.requestedCount);
        postRmaxIdx = selectedIdx - maxidx + 1;
        collapseTimesUs = (texp(collapseInfo.collapseIdx) - ...
            texp(maxidx)) .* 1e6;
        fprintf(['Fit window ends at radial collapse %d, post-Rmax ', ...
            'sample %d (t = %.6g us).\n'], collapseInfo.requestedCount, ...
            postRmaxIdx, collapseTimesUs(collapseInfo.requestedCount));
        fprintf('  Detected collapse times (us after Rmax): %s\n', ...
            num2str(collapseTimesUs(1:collapseInfo.requestedCount).', ...
            ' %.6g'));
        return
    end

    fprintf(['Detected only %d of %d requested radial collapses; using ', ...
        'the radius-equilibrium fallback.\n'], collapseInfo.detectedCount, ...
        collapseInfo.requestedCount);
    printRadiusEquilibriumCutoff(equilibriumInfo, texp, maxidx);
end

function printRadiusEquilibriumCutoff(info, texp, maxidx)
    if ~info.enabled
        fprintf(['Radius-equilibrium cutoff is disabled, so the full ', ...
            'post-Rmax record is used.\n']);
        return
    end
    if ~info.found
        fprintf(['Radius did not remain within %.3g m of Req for %d ', ...
            'consecutive samples; the full post-Rmax record is used.\n'], ...
            info.radiusTolerance, info.nConsecutive);
        return
    end

    startPostRmax = info.startIdx - maxidx + 1;
    confirmationPostRmax = info.confirmationIdx - maxidx + 1;
    fprintf(['Radius equilibrium begins at post-Rmax sample %d and is ', ...
        'confirmed at sample %d after %d consecutive points ', ...
        '(t = %.6g us, |R-Req| <= %.3g m).\n'], startPostRmax, ...
        confirmationPostRmax, info.nConsecutive, ...
        (texp(info.confirmationIdx) - texp(maxidx)) * 1e6, ...
        info.radiusTolerance);
    fprintf('  Excluding %d later samples from simulation and loss.\n', ...
        info.nExcludedSamples);
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

function value = simOpt(simOpts, fieldName, defaultValue)
    if isfield(simOpts, fieldName) && ~isempty(simOpts.(fieldName))
        value = simOpts.(fieldName);
    else
        value = defaultValue;
    end
end

function plotCollapseDetection(tpost_nd, Rpost_data, collapseInfo, ...
        fitEndPostRmaxIdx)
    figure('Name', 'Optimization-window and collapse check')
    plot(tpost_nd, Rpost_data, 'o-', 'DisplayName', 'experimental radius')
    hold on
    collapseIdx = find(Rpost_data == min(Rpost_data(1: ...
        collapseInfo.confirmationIdx)), 1, 'first');
    plot(tpost_nd(collapseIdx), Rpost_data(collapseIdx), 'kp', ...
        'MarkerFaceColor', 'y', 'MarkerSize', 11, ...
        'DisplayName', 'first collapse')
    plot(tpost_nd(collapseInfo.confirmationIdx), ...
        Rpost_data(collapseInfo.confirmationIdx), 'ks', ...
        'MarkerFaceColor', 'c', 'MarkerSize', 8, ...
        'DisplayName', 'rebound confirmation')
    xline(tpost_nd(fitEndPostRmaxIdx), 'r--', 'optimization ends', ...
        'LineWidth', 1.5, 'LabelVerticalAlignment', 'middle', ...
        'DisplayName', 'optimization cutoff')
    xlabel("t^*")
    ylabel("R/R_{max}")
    title('Experimental radius and optimization window')
    legend('Location', 'best')
    grid on
end

function plotInitialConditionCheck(texp, maxidx, fitEndIdx, amps_og, ...
        amps, epnm0, epnmd0, tc)
    figure('Name', 'Initial-condition check')
    nmodes = size(amps, 1);
    plotl = ceil(sqrt(nmodes));
    t0 = texp(maxidx);
    timeNd = (texp - t0) ./ tc;
    fitEndTimeNd = timeNd(fitEndIdx);
    for ii = 1:nmodes
        subplot(plotl, plotl, ii)
        plot(timeNd, amps_og(ii, :), 'o')
        hold on
        plot(timeNd, amps(ii, :), '-')
        plot(timeNd, 0 .* timeNd + epnm0(ii), ':')
        plot(timeNd, timeNd .* epnmd0(ii) + epnm0(ii), '-')
        xline(0, 'k:', 'LineWidth', 1)
        if ii == 1
            xline(fitEndTimeNd, 'r--', 'optimization ends', ...
                'LineWidth', 1.5, 'LabelVerticalAlignment', 'middle')
        else
            xline(fitEndTimeNd, 'r--', 'LineWidth', 1.5)
        end
        ylim([min(amps_og(ii, :)) max(amps_og(ii, :))])
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
