% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% s = imr.fluid_properties(profile, P, T): liquid and ideal-vapor SI slices.
% s = imr.fluid_properties(profile, [], T): saturation pressure and derivative.
% profile is a data struct or a zero-argument MATLAB function handle.
% P (Pa) and T (K) accept compatible arrays via implicit expansion. A zero
% vapor pressure gives rho = 0 and v = Inf. No property calls occur in ODEs.
function s = fluid_properties(profile, P, T)
    if isa(profile,'function_handle')
        profile = profile();
    end
    if ~(isstruct(profile) && isscalar(profile))
        error('IMR:fluidProfile', 'A fluid profile must be a scalar struct or function handle.');
    end
    validateattributes(T, {'numeric'}, {'real','finite','nonempty','positive'}, mfilename, 'T');
    T = double(T);
    if isempty(P)
        s = f_saturation(profile.saturation, T);
        return
    end
    validateattributes(P, {'numeric'}, {'real','finite','nonempty'}, mfilename, 'P');
    f_range(P, profile.pressure_range, 'IMR:fluidPressure', 'P (Pa)');
    f_range(T, profile.temperature_range, 'IMR:fluidTemperature', 'T (K)');
    if any(T(:) < profile.sound_temperature_range(1) | T(:) > profile.sound_temperature_range(2))
        warning('IMR:fluidExtrapolation', '%s sound speed is outside its measured temperature range.', profile.name);
    end
    broadcast = double(P) + T;
    shape = size(broadcast);
    P = double(P) + zeros(shape);
    T = T + zeros(shape);
    l = profile.liquid;
    v = profile.vapor;
    s.P = P;
    s.T = T;
    switch l.model
        case 'constant'
        s.liquid.rho = l.rho + zeros(shape);
        s.liquid.c = l.c + zeros(shape);
        case 'murnaghan'
        rho0 = l.rho0_a*T + l.rho0_b;
        c0 = l.c0_0 + l.c0_1*T;
        K0 = rho0.*c0.^2/l.gamma;
        n = 1 + 2*K0*l.s1./c0;
        base = 1 + n.*P./K0;
        if any(base(:) <= 0)
            error('IMR:fluidDomain', 'The Murnaghan base must be positive.');
        end
        s.liquid.rho0 = rho0;
        s.liquid.c0 = c0;
        s.liquid.K0 = K0;
        s.liquid.K0_prime = n;
        s.liquid.rho = rho0.*base.^(1./n);
        s.liquid.c = sqrt(l.gamma*(K0 + n.*P)./s.liquid.rho);
        otherwise
        error('IMR:fluidModel', 'Unknown liquid model: %s.', l.model);
    end
    s.liquid.h = l.Cp*T + l.q;
    switch v.model
        case 'ideal'
        s.vapor.rho = P./(v.Rs*T);
        s.vapor.v = v.Rs*T./P;
        s.vapor.c = sqrt(v.gamma*v.Rs*T);
        otherwise
        error('IMR:fluidModel', 'Unknown vapor model: %s.', v.model);
    end
    switch profile.caloric.model
        case 'linear_enthalpy'
        s.vapor.h = v.Cp*T + v.q;
        s.latent_heat = s.vapor.h - s.liquid.h;
        case 'constant_latent_heat'
        s.latent_heat = profile.caloric.latent_heat + zeros(shape);
        s.vapor.h = s.liquid.h + s.latent_heat;
        otherwise
        error('IMR:fluidModel', 'Unknown caloric model: %s.', profile.caloric.model);
    end
end

function s = f_saturation(law, T)
    f_range(T, law.temperature_range, 'IMR:fluidSaturation', 'Saturation T (K)');
    s.model = law.model;
    s.c2 = 0; % No exponential wall law for a table.
    switch law.model
        case 'exp_inverse_temperature'
        [s.pressure,s.dpressure_dT] = imr.physics.saturation_pressure(T,law.c1,law.c2);
        s.c2 = law.c2;
        case 'table_logp_inverseT'
        % Normalize table orientation so each query is compared with every
        % knot, including when a custom profile supplies column vectors.
        x = 1./flipud(law.T(:))';
        y = log(flipud(law.P(:)))';
        query = 1./T(:);
        logp = interp1(x, y, query, 'linear');
        % Right-hand T derivative at interior knots; left-hand at Tmax.
        slopes = diff(y(:))./diff(x(:));
        index = sum(query > x,2);
        index = min(max(index,1),numel(slopes));
        s.pressure = reshape(exp(logp),size(T));
        s.dpressure_dT = reshape(-slopes(index(:)).*exp(logp)./T(:).^2,size(T));
        otherwise
        error('IMR:fluidModel', 'Unknown saturation model: %s.', law.model);
    end
end

function f_range(value, bounds, id, label)
    if any(value(:) < bounds(1) | value(:) > bounds(2))
        error(id, '%s must be in [%g, %g].', label, bounds(1), bounds(2));
    end
end
