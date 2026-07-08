function loss = f_optimize_model_to_data_loss(z, xData, paramSpec, simOpts)
% Scalar normalized loss used by bayesopt.

[yModel, runInfo] = f_optimize_model_to_data_predict(z, xData, paramSpec, ...
    simOpts);

if ~runInfo.success
    rawLoss = simOpt(simOpts, 'FailurePenalty', 1e6);
    loss = transformLossForOptimizer(rawLoss, simOpts);
    writeEvaluationDiagnostic(simOpts, runInfo, rawLoss, loss);
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
    runInfo.success = false;
    runInfo.identifier = 'ModelToData:NonfiniteObjective';
    runInfo.message = 'Objective value was nonfinite after loss transform.';
    writeEvaluationDiagnostic(simOpts, runInfo, rawLoss, loss);
    if simOpt(simOpts, 'PrintFailures', false)
        fprintf(['[%s] Penalty loss %.6g (raw %.6g): nonfinite ', ...
            'objective value.\n'], currentWorkerLabel(), loss, rawLoss);
    end
else
    writeEvaluationDiagnostic(simOpts, runInfo, rawLoss, loss);
    if simOpt(simOpts, 'PrintSuccess', false)
        fprintf(['[%s] Loss %.6g (raw %.6g) completed in %.3f s ', ...
            '(|dt_R^*|=%.3g, |dt_ep^*|=%.3g).\n'], ...
            currentWorkerLabel(), loss, rawLoss, runInfo.elapsed, ...
            runInfo.maxRadialNearestTimeMismatch, ...
            runInfo.maxPerturbationNearestTimeMismatch);
    end
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

function writeEvaluationDiagnostic(simOpts, runInfo, rawLoss, loss)
logFile = simOpt(simOpts, 'DiagnosticLogFile', '');
if isempty(logFile)
    return
end

try
    logDir = fileparts(logFile);
    if ~isempty(logDir) && ~isfolder(logDir)
        mkdir(logDir);
    end
    fid = fopen(logFile, 'a');
    if fid < 0
        return
    end
    cleanup = onCleanup(@() fclose(fid));
    [G, alph, mu, ani1, ani2] = paramsToLogValues(runInfo);
    fprintf(fid, ['%s\t%s\t%d\t%.16g\t%.16g', ...
        '\t%.6g\t%.16g\t%.6g\t%.6g\t%.6g', ...
        '\t%d\t%s\t%s\t%.16g\t%.16g\t%.16g\t%.16g\t%.16g\n'], ...
        char(datetime("now", "Format", "yyyy-MM-dd'T'HH:mm:ss")), ...
        currentWorkerLabel(), logical(runInfo.success), ...
        loss, rawLoss, runInfo.elapsed, runInfo.completedTime, ...
        runInfo.maxRequestedReturnedTimeMismatch, ...
        runInfo.maxRadialNearestTimeMismatch, ...
        runInfo.maxPerturbationNearestTimeMismatch, ...
        logical(runInfo.timedOut), sanitizeLogText(runInfo.identifier), ...
        sanitizeLogText(runInfo.message), G, alph, mu, ani1, ani2);
    clear cleanup
catch
    % Diagnostics should never change the optimizer objective.
end
end

function [G, alph, mu, ani1, ani2] = paramsToLogValues(runInfo)
G = NaN;
alph = NaN;
mu = NaN;
ani1 = NaN;
ani2 = NaN;
if ~isfield(runInfo, 'params') || ~isstruct(runInfo.params) || ...
        isempty(fieldnames(runInfo.params))
    return
end
params = runInfo.params;
try
    G = params.G;
    alph = params.alph;
    mu = params.mu;
    ani1 = params.ani(1);
    ani2 = params.ani(2);
catch
end
end

function text = sanitizeLogText(text)
if isempty(text)
    text = '';
elseif isstring(text)
    text = char(text);
elseif ~ischar(text)
    text = char(string(text));
end
text = regexprep(text, '[\r\n\t]+', ' ');
end

function label = currentWorkerLabel()
task = getCurrentTask();
if isempty(task)
    label = 'client';
else
    label = sprintf('worker-%d', task.ID);
end
end
