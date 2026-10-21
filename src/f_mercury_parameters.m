% file f_mercury_parameters.m
% brief returns the versioned liquid and vapor mercury fit coefficients

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% All coefficients are SI. No fit, file read or external runtime is needed.
% See docs/model/mercury.md for equations, units and original data citations.
% Gamma is provisional and s2 is a surface-fit diagnostic. Neither enters
% the implemented liquid density/acoustic evaluators or the KM equation.
function p = f_mercury_parameters()
    p.source = 'mrodrig6/mercury_eos';
    p.commit = '5d51892ee4eef7d898d33ac81a02bee52d554165';
    p.temperature_range = [275 550];
    p.sound_temperature_range = [296 513];
    p.pressure_range = [0 6.23e9];

    p.liquid.gamma = 1.129196028264567;
    p.liquid.Gamma = 0.129196028264567;
    p.liquid.Cv = 121.07661322485443;
    p.liquid.Cp = 136.71923076923076;
    p.liquid.q = 5763.949807692312;
    p.liquid.c0_0 = 1597.66970844695;
    p.liquid.c0_1 = -0.46575634646134273;
    p.liquid.s1 = 1.8089270168866193e-7;
    p.liquid.s2 = -9.353830846791366e-18;
    p.liquid.rho0_a = -2.4016642246642292;
    p.liquid.rho0_b = 14248.445722568176;
    p.liquid.K0_at_Tref = 23194789101.321957;
    p.liquid.K0_prime_at_Tref = 6.970306780987923;
    p.liquid.Tref = 412.5;
    p.liquid.cs_rms_surface_percent = 0.49595540334060295;
    p.liquid.cs_rms_canonical_percent = 0.7124894179683818;

    p.vapor.gamma = 1.6668832801308178;
    p.vapor.Gamma = 0.6668832801308178;
    p.vapor.Cv = 62.154204843565914;
    p.vapor.Cp = 103.60380484356591;
    p.vapor.q = 321831.5461763998;
    p.vapor.Rs = 41.4496;
    p.vapor.cs_rms_percent = 0.03343476679216099;
    % Fokin, Popov and Naurzakov (2011), as transcribed in source data.
    p.vapor.saturation_T = [298.15 300 350 400 450 500 550];
    p.vapor.saturation_P = [0.26117 0.30427 10.148 138.75 ...
        1049.1 5247.7 19457];
end
