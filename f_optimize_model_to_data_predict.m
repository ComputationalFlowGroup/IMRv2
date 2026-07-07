function [y_out, runInfo, sim] = f_optimize_model_to_data_predict( ...
    z, xData, paramSpec, simOpts)
% Run one IMRv2 simulation and return the weighted model vector.

params = f_unpack_model_to_data_params(z, paramSpec);
wantSimOutput = nargout >= 3;
runInfo = struct('success', false, 'message', '', 'identifier', '', ...
    'elapsed', NaN, 'completedTime', NaN, 'timedOut', false, ...
    'maxRequestedReturnedTimeMismatch', NaN, ...
    'radialLossTimeNd', [], 'perturbationLossTimeNd', [], ...
    'radialLossReturnedTimeNd', [], 'perturbationLossReturnedTimeNd', [], ...
    'maxRadialNearestTimeMismatch', NaN, ...
    'maxPerturbationNearestTimeMismatch', NaN, ...
    'firstCollapseTimeNd', NaN, 'nRadialLossTimes', 0, ...
    'nPerturbationLossTimes', 0);
sim = struct('success', false);

ticRun = tic;
try
    [t, R, epnm, hardTimedOut] = runForwardWithOptionalHardTimeout( ...
        params, xData, simOpts);
catch ME
    runInfo.elapsed = toc(ticRun);
    runInfo.identifier = ME.identifier;
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
if numel(t) == numel(xData.tfit_nd)
    runInfo.maxRequestedReturnedTimeMismatch = max(abs(t(:) - xData.tfit_nd(:)));
    runInfo.maxRadialNearestTimeMismatch = runInfo.maxRequestedReturnedTimeMismatch;
    runInfo.maxPerturbationNearestTimeMismatch = max(abs(t(epRows) - ...
        xData.tfit_nd(epRows)));
else
    runInfo.maxRequestedReturnedTimeMismatch = Inf;
    runInfo.maxRadialNearestTimeMismatch = Inf;
    runInfo.maxPerturbationNearestTimeMismatch = Inf;
end

if numel(t) < 2 || max(t) < xData.tfit_nd(end) - timeTol || ...
        any(~isfinite(t(:))) || any(~isfinite(R(:))) || ...
        any(~isfinite(epnm(:)))
    if runInfo.timedOut
        runInfo.message = 'Simulation exceeded the per-evaluation wall-time limit.';
    else
        runInfo.message = 'Simulation failed or ended before the fit window.';
    end
    y_out = failedPrediction(xData, simOpts);
    return
end

timesAligned = numel(t) == numel(xData.tfit_nd) && ...
    runInfo.maxRequestedReturnedTimeMismatch <= timeTol;
if ~timesAligned
    runInfo.identifier = 'IMR:UnexpectedTimeVector';
    runInfo.message = sprintf(['Solver did not return the requested ', ...
        'experimental time vector; max |dt*| = %.3g.'], ...
        runInfo.maxRequestedReturnedTimeMismatch);
    y_out = failedPrediction(xData, simOpts);
    return
end

RFit = R(:);
epFit = epnm;

if any(~isfinite(RFit(:))) || any(~isfinite(epFit(:)))
    runInfo.message = 'Solver returned nonfinite values on the experimental time grid.';
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
