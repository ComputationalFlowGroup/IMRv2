function [y_out, runInfo, sim] = f_optimize_model_to_data_predict( ...
    z, xData, paramSpec, simOpts)
% Run one IMRv2 simulation and return the weighted model vector.

wantSimOutput = nargout >= 3;
runInfo = struct('success', false, 'message', '', 'identifier', '', ...
    'elapsed', NaN, 'completedTime', NaN, 'timedOut', false, ...
    'maxRequestedReturnedTimeMismatch', NaN, ...
    'radialLossTimeNd', [], 'perturbationLossTimeNd', [], ...
    'radialLossReturnedTimeNd', [], 'perturbationLossReturnedTimeNd', [], ...
    'maxRadialNearestTimeMismatch', NaN, ...
    'maxPerturbationNearestTimeMismatch', NaN, ...
    'firstCollapseTimeNd', NaN, 'nRadialLossTimes', 0, ...
    'nPerturbationLossTimes', 0, 'params', struct(), ...
    'usedNearestTimeExtraction', false);
sim = struct('success', false);

try
    params = f_unpack_model_to_data_params(z, paramSpec);
catch ME
    runInfo.identifier = nonemptyIdentifier(ME, ...
        'ModelToData:ParameterUnpackFailed');
    runInfo.message = exceptionSummary(ME);
    y_out = failedPrediction(xData, simOpts);
    return
end
runInfo.params = params;

ticRun = tic;
try
    [t, R, epnm, hardTimedOut] = runForwardWithOptionalHardTimeout( ...
        params, xData, simOpts);
catch ME
    runInfo.elapsed = toc(ticRun);
    runInfo.identifier = nonemptyIdentifier(ME, ...
        'ModelToData:ForwardSolveFailed');
    runInfo.message = exceptionSummary(ME);
    runInfo.timedOut = strcmp(ME.identifier, 'IMR:MaxWallTimeExceeded');
    y_out = failedPrediction(xData, simOpts);
    return
end

runInfo.elapsed = toc(ticRun);
if isempty(t)
    runInfo.completedTime = NaN;
else
    runInfo.completedTime = max(t);
end
timeTol = simOpt(simOpts, 'TimeMatchTolerance', ...
    max(1e-10, 1e-8 * max(1, abs(xData.tfit_nd(end)))));
runInfo.timedOut = hardTimedOut || (simOpt(simOpts, 'MaxWallTime', 0) > 0 && ...
    runInfo.elapsed >= 0.95 * simOpt(simOpts, 'MaxWallTime', 0) && ...
    (isempty(t) || max(t) < xData.tfit_nd(end) - timeTol));

if wantSimOutput
    sim = struct('success', false, 't', t(:), 'R', R(:), 'epnm', epnm, ...
        'params', params);
end
runInfo.radialLossTimeNd = xData.tfit_nd(:);
runInfo.nRadialLossTimes = numel(xData.tfit_nd);
epRows = perturbationFitRows(xData);
runInfo.perturbationLossTimeNd = xData.tfit_nd(epRows);
runInfo.nPerturbationLossTimes = numel(epRows);
if simOpt(simOpts, 'VerifyTimeExtraction', false) || ...
        simOpt(simOpts, 'PrintSuccess', false)
    [radialReturnedTime, radialAbsDt] = nearestReturnedTimes(t(:), ...
        runInfo.radialLossTimeNd);
    [perturbationReturnedTime, perturbationAbsDt] = nearestReturnedTimes( ...
        t(:), runInfo.perturbationLossTimeNd);
    runInfo.radialLossReturnedTimeNd = radialReturnedTime;
    runInfo.perturbationLossReturnedTimeNd = perturbationReturnedTime;
    runInfo.maxRadialNearestTimeMismatch = maxFinite(radialAbsDt);
    runInfo.maxPerturbationNearestTimeMismatch = maxFinite( ...
        perturbationAbsDt);
end
if isfield(xData, 'firstCollapseTimeNd')
    runInfo.firstCollapseTimeNd = xData.firstCollapseTimeNd;
end
if numel(t) < 2 || any(~isfinite(t(:))) || ...
        max(t) < xData.tfit_nd(end) - timeTol
    if runInfo.timedOut
        runInfo.message = 'Simulation exceeded the per-evaluation wall-time limit.';
    else
        runInfo.message = ['Simulation failed, returned nonfinite times, ', ...
            'or ended before the fit window.'];
    end
    y_out = failedPrediction(xData, simOpts);
    return
end

[RFit, epFit, runInfo] = extractRequestedTimeSamples(t, R, epnm, xData, ...
    epRows, simOpts, runInfo, timeTol);
if isempty(RFit)
    runInfo.identifier = 'IMR:UnexpectedTimeVector';
    y_out = failedPrediction(xData, simOpts);
    return
end

if any(~isfinite(RFit(:)))
    runInfo.identifier = 'IMR:NonfiniteRadialLossData';
    runInfo.message = 'Solver returned nonfinite radial values on the loss time grid.';
    y_out = failedPrediction(xData, simOpts);
    return
end

epLoss = epFit(epRows, xData.fitModeIdx);
if any(~isfinite(epLoss(:)))
    runInfo.identifier = 'IMR:NonfinitePerturbationLossData';
    runInfo.message = ['Solver returned nonfinite values for perturbation ', ...
        'modes used in the loss window.'];
    y_out = failedPrediction(xData, simOpts);
    return
end

y_out = [xData.aR .* RFit(:); reshape(epFit(epRows, xData.fitModeIdx) .* ...
    xData.aEP(xData.fitModeIdx), [], 1)];
runInfo.success = true;
runInfo.message = 'Simulation completed.';
if wantSimOutput
    sim.success = true;
end
end

function [RFit, epFit, runInfo] = extractRequestedTimeSamples(t, R, epnm, ...
    xData, epRows, simOpts, runInfo, timeTol)
RFit = [];
epFit = [];
requested = xData.tfit_nd(:);
returned = t(:);

if numel(returned) == numel(requested)
    dt = abs(returned - requested);
    runInfo.maxRequestedReturnedTimeMismatch = max(dt);
    runInfo.maxRadialNearestTimeMismatch = runInfo.maxRequestedReturnedTimeMismatch;
    runInfo.maxPerturbationNearestTimeMismatch = max(dt(epRows));
    if runInfo.maxRequestedReturnedTimeMismatch <= timeTol
        RFit = R(:);
        epFit = epnm;
        return
    end
end

if ~simOpt(simOpts, 'AllowNearestTimeExtraction', true)
    runInfo.message = sprintf(['Solver did not return the requested ', ...
        'experimental time vector; max |dt*| = %.3g.'], ...
        runInfo.maxRequestedReturnedTimeMismatch);
    return
end

[nearestIdx, nearestDt] = nearestReturnedTimeIndices(returned, requested);
nearestTol = simOpt(simOpts, 'NearestTimeTolerance', ...
    max(timeTol, 1e-6 * max(1, abs(requested(end)))));
runInfo.maxRequestedReturnedTimeMismatch = maxFinite(nearestDt);
runInfo.maxRadialNearestTimeMismatch = runInfo.maxRequestedReturnedTimeMismatch;
runInfo.maxPerturbationNearestTimeMismatch = maxFinite(nearestDt(epRows));

if any(~isfinite(nearestDt)) || ...
        runInfo.maxRequestedReturnedTimeMismatch > nearestTol
    runInfo.message = sprintf(['Solver did not return usable samples near ', ...
        'the requested experimental time vector; max nearest |dt*| = %.3g ', ...
        '(tol %.3g).'], runInfo.maxRequestedReturnedTimeMismatch, ...
        nearestTol);
    return
end

RFit = R(nearestIdx);
epFit = epnm(nearestIdx, :);
runInfo.usedNearestTimeExtraction = true;
end

function epRows = perturbationFitRows(xData)
if isfield(xData, 'epFitIdx') && ~isempty(xData.epFitIdx)
    epRows = xData.epFitIdx(:);
else
    epRows = (1:numel(xData.tfit_nd)).';
end
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
    [absDt(ii), nearestIdx] = min(abs(returnedTimes - requestedTimes(ii)));
    nearestTimes(ii) = returnedTimes(nearestIdx);
end
end

function [nearestIdx, absDt] = nearestReturnedTimeIndices(returnedTimes, ...
    requestedTimes)
returnedTimes = returnedTimes(:);
requestedTimes = requestedTimes(:);
nearestIdx = nan(size(requestedTimes));
absDt = nan(size(requestedTimes));
if isempty(returnedTimes) || isempty(requestedTimes)
    return
end

for ii = 1:numel(requestedTimes)
    [absDt(ii), nearestIdx(ii)] = min(abs(returnedTimes - requestedTimes(ii)));
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

function [t, R, epnm, hardTimedOut] = runForwardWithOptionalHardTimeout( ...
    params, xData, simOpts)
hardTimedOut = false;
maxWallTime = simOpt(simOpts, 'MaxWallTime', 0);
useHardTimeout = simOpt(simOpts, 'UseHardTimeout', false);

if useHardTimeout && maxWallTime > 0
    future = parfeval(backgroundPool, @runForwardSimulation, 3, ...
        params, xData, simOpts);
    cleanup = onCleanup(@() cancelUnfinishedFuture(future));
    pollInterval = simOpt(simOpts, 'TimeoutPollInterval', 0.25);
    waitStart = tic;
    while ~strcmp(future.State, 'finished') && toc(waitStart) < maxWallTime
        pause(pollInterval);
    end

    if ~strcmp(future.State, 'finished')
        cancel(future);
        hardTimedOut = true;
        t = [];
        R = [];
        epnm = [];
        return
    end

    [t, R, epnm] = fetchOutputs(future);
else
    [t, R, epnm] = runForwardSimulation(params, xData, simOpts);
end
end

function [t, R, epnm] = runForwardSimulation(params, xData, simOpts)
[t, R, epnm] = f_call_IMRv2_exp(xData.Rmax, xData.Req, ...
    xData.epnm0, xData.epnmd0, xData.epnmeq, xData.n, xData.m, ...
    params.mu, params.G, params.alph, params.ani, xData.sig, ...
    xData.p_a, xData.f_a, xData.tf_nd, numel(xData.tfit_nd), ...
    xData.ultra, ...
    'RelTol', simOpt(simOpts, 'RelTol', 1e-6), ...
    'AbsTol', simOpt(simOpts, 'AbsTol', 1e-7), ...
    'Nt', simOpt(simOpts, 'Nt', 75), ...
    'Method', simOpt(simOpts, 'Method', 45), ...
    'MaxWallTime', simOpt(simOpts, 'MaxWallTime', 0), ...
    'OutputTimeNd', xData.tfit_nd);
end

function cancelUnfinishedFuture(future)
if ~isempty(future) && ~strcmp(future.State, 'finished')
    cancel(future);
end
end

function message = exceptionSummary(ME)
message = ME.message;
if ~isempty(ME.cause)
    causeMessages = cellfun(@(cause) cause.message, ME.cause, ...
        'UniformOutput', false);
    message = strjoin([{message}, causeMessages(:).'], ' | Cause: ');
end
end

function identifier = nonemptyIdentifier(ME, fallback)
identifier = ME.identifier;
if isempty(identifier)
    identifier = fallback;
end
end

function value = simOpt(simOpts, fieldName, defaultValue)
if isfield(simOpts, fieldName) && ~isempty(simOpts.(fieldName))
    value = simOpts.(fieldName);
else
    value = defaultValue;
end
end

function y_out = failedPrediction(xData, simOpts)
penalty = simOpt(simOpts, 'FailurePenalty', 1e6);
residual = sqrt(penalty);
y_out = xData.y_data + residual .* ones(size(xData.y_data));
end
