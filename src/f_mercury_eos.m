% file f_mercury_eos.m
% brief evaluates the fitted mercury liquid and ideal-vapor properties

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Usage: hg = f_mercury_eos(P, T), with absolute P in Pa and T in K.
% Scalars and compatible arrays expand using MATLAB implicit expansion.
% Output fields are SI. At P = 0, ideal vapor v is Inf and rho is zero.
% This provides property slices, not a complete thermodynamic p(rho,e).
function hg = f_mercury_eos(P, T)
    validateattributes(P, {'numeric'}, {'real','finite','nonempty'}, ...
        mfilename, 'P');
    validateattributes(T, {'numeric'}, {'real','finite','nonempty'}, ...
        mfilename, 'T');
    p = f_mercury_parameters();
    if any(P(:) < p.pressure_range(1) | P(:) > p.pressure_range(2))
        error('IMR:mercuryPressure', ...
            'Mercury P must be in [0, 6.23e9] Pa. Tension is unsupported.');
    end
    if any(T(:) < p.temperature_range(1) | T(:) > p.temperature_range(2))
        error('IMR:mercuryTemperature', 'Mercury T must be in [275, 550] K.');
    end
    if any(T(:) < p.sound_temperature_range(1) | T(:) > p.sound_temperature_range(2))
        warning('IMR:mercuryExtrapolation', ...
            'Mercury sound speed is extrapolated outside 296-513 K.');
    end
    % Expand all reported fields to the same shape, including caloric ones.
    shape = size(double(P) + double(T));
    P = double(P) + zeros(shape);
    T = double(T) + zeros(shape);
    l = p.liquid;
    v = p.vapor;
    rho0 = l.rho0_a*T + l.rho0_b;
    c0 = l.c0_0 + l.c0_1*T;
    K0 = rho0.*c0.^2/l.gamma;
    n = 1 + 2*K0*l.s1./c0;
    base = 1 + n.*P./K0;
    if any(base(:) <= 0)
        error('IMR:mercuryDomain', 'Mercury Murnaghan base must be positive.');
    end
    hg.P = P;
    hg.T = T;
    hg.liquid.rho0 = rho0;
    hg.liquid.c0 = c0;
    hg.liquid.K0 = K0;
    hg.liquid.K0_prime = n;
    hg.liquid.rho = rho0.*base.^(1./n);
    hg.liquid.c = sqrt(l.gamma*(K0 + n.*P)./hg.liquid.rho);
    hg.liquid.h = l.Cp*T + l.q;
    hg.vapor.rho = P./(v.Rs*T);
    hg.vapor.v = v.Rs*T./P;
    hg.vapor.c = sqrt(v.gamma*v.Rs*T);
    hg.vapor.h = v.Cp*T + v.q;
    hg.latent_heat = hg.vapor.h - hg.liquid.h;
end
