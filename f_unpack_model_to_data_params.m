function params = f_unpack_model_to_data_params(z, paramSpec)
% Convert optimizer variables into physical model parameters.

params = struct('G', NaN, 'alph', NaN, 'mu', NaN, 'ani', [NaN, NaN]);
activeIdx = 0;
for ii = 1:numel(paramSpec)
    optimizeThis = specField(paramSpec(ii), 'optimize', true);
    if optimizeThis
        activeIdx = activeIdx + 1;
        value = optimizerValue(z, paramSpec(ii), activeIdx);
        if specField(paramSpec(ii), 'logScale', false)
            value = 10.^value;
        end
    else
        value = specField(paramSpec(ii), 'fixedValue', NaN);
    end

    switch paramSpec(ii).name
        case 'G'
            params.G = value;
        case 'alph'
            params.alph = value;
        case 'mu'
            params.mu = value;
        case 'ani1'
            params.ani(1) = value;
        case 'ani2'
            params.ani(2) = value;
        otherwise
            error('Unknown optimization parameter "%s".', paramSpec(ii).name);
    end
end

validateOptimizerVariableCount(z, activeIdx);
end

function value = optimizerValue(z, spec, activeIdx)
if istable(z)
    value = z.(spec.variableName)(1);
elseif isstruct(z)
    value = z.(spec.variableName);
else
    values = z(:).';
    value = values(activeIdx);
end
end

function validateOptimizerVariableCount(z, nActive)
if istable(z) || isstruct(z)
    return
end
if numel(z) ~= nActive
    error('Expected %d active optimizer variables, received %d.', ...
        nActive, numel(z));
end
end

function value = specField(spec, fieldName, defaultValue)
if isfield(spec, fieldName)
    value = spec.(fieldName);
else
    value = defaultValue;
end
end
