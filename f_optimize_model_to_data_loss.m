function loss = f_optimize_model_to_data_loss(z, xData, paramSpec, simOpts)
% Scalar normalized loss used by bayesopt.

[yModel, runInfo] = f_optimize_model_to_data_predict(z, xData, paramSpec, ...
    simOpts);

if ~runInfo.success
    rawLoss = simOpt(simOpts, 'FailurePenalty', 1e6);
    loss = transformLossForOptimizer(rawLoss, simOpts);
    if simOpt(simOpts, 'PrintFailures', false)
        fprintf(['[%s] Penalty loss %.6g (raw %.6g): %s %s ', ...
            '(elapsed %.3f s, completed t*=%.6g)\n'], ...
            currentWorkerLabel(), loss, rawLoss, runInfo.identifier, ...
            runInfo.message, runInfo.elapsed, runInfo.completedTime);
    end
    return
end

rawLoss = sqrt(sum((xData.y_data - yModel).^2)) / norm(xData.y_data);
loss = transformLossForOptimizer(rawLoss, simOpts);
if ~isfinite(loss)
    rawLoss = simOpt(simOpts, 'FailurePenalty', 1e6);
    loss = transformLossForOptimizer(rawLoss, simOpts);
    if simOpt(simOpts, 'PrintFailures', false)
        fprintf(['[%s] Penalty loss %.6g (raw %.6g): nonfinite ', ...
            'objective value.\n'], currentWorkerLabel(), loss, rawLoss);
    end
elseif simOpt(simOpts, 'PrintSuccess', false)
    fprintf(['[%s] Loss %.6g (raw %.6g) completed in %.3f s ', ...
        '(|dt_R^*|=%.3g, |dt_ep^*|=%.3g).\n'], ...
        currentWorkerLabel(), loss, rawLoss, runInfo.elapsed, ...
        runInfo.maxRadialNearestTimeMismatch, ...
        runInfo.maxPerturbationNearestTimeMismatch);
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

function value = simOpt(simOpts, fieldName, defaultValue)
if isfield(simOpts, fieldName) && ~isempty(simOpts.(fieldName))
    value = simOpts.(fieldName);
else
    value = defaultValue;
end
end

function label = currentWorkerLabel()
task = getCurrentTask();
if isempty(task)
    label = 'client';
else
    label = sprintf('worker-%d', task.ID);
end
end
