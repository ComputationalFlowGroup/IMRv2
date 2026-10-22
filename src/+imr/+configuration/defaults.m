% file defaults.m
% brief contains function imr.configuration.defaults, the default (dimensionless) inputs

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief All quantities are dimensionless. Lengths are scaled by the initial
% radius R0, pressures by the far-field pressure P8, velocities by
% Uc = sqrt(P8/rho8), time by t0 = R0/Uc, and temperatures by the far-field
% temperature T8. The values below correspond to a 50 um bubble in water
% at 1 atm and 298.15 K; src/imr.nondimensionalize.m converts dimensional
% (SI) properties into these numbers.
% Numerical/model defaults live here; material data live in imr.fluids.water.
% Physical groups below are derived through the SI converter. imr.configuration.prepare_case
% starts from these values and overrides explicit name/value pairs.
function d = defaults()
    % One SI material source and one conversion path; no duplicated physical
    % constants. This call does not recurse: the converter does not use d.
    pairs = imr.nondimensionalize();
    physical = struct(pairs{:});
    d = struct();

    % equation options

    % 1 : RP, 2 : K-M P, 3 : K-M E, 4 : Gil, 5 : K-M MG, 6 : Gil MG
    d.radial = 1;
    % 0 : polytropic assumption, 1: thermal PDE in bubble
    d.bubtherm = 0;
    % 0 : cold fluid, 1: warm fluid assumption
    d.medtherm = 0;
    % 0 : none, 1 : NHKV, 2 : qNHKV, 3 : linear Maxwell, Jeffreys, Zener
    % (neo-Hookean), 4 : same with quadratic neo-Hookean, 5 : UCM or OldB
    d.stress = 0;
    % this value must be [0, 0.5]
    d.eps3 = 0;
    % 0 : ignore vapor pressure, 1 : vapor pressure
    d.vapor = 0;
    % mass transfer, default is no mass transfer
    d.masstrans = 0;

    % solver options

    % final time
    d.TFin = 3.903439778613622;
    % time span vector
    d.TVector = linspace(0,d.TFin,128);
    % stiff time integrator: 15 ode15s, 23 ode23tb, 3 CVODES (the ode
    % object, MATLAB R2024a or later)
    d.method = 23;
    % force spectral collocation solution
    d.spectral = 0;
    % minimum number of timesteps
    d.divisions = 0;
    % number of points in bubble, thermal PDE
    d.Nt = 25;
    % number of points outside of bubble, thermal PDE
    d.Mt = 25;
    % number of stress modes of imr.solve_spectral: the stress 5 field
    % (spectral = 1), and the stress 3, 4 polymer stress field with medium
    % heat transfer
    d.Nv = 150;
    % grid stretching length of the stress field of imr.solve_spectral
    d.Lv = 3;
    % grid stretching length of the medium temperature field
    d.Lt = 2;
    % relative tolerance of the time integrator
    d.rtol = 1e-8;
    % absolute tolerance of the time integrator
    d.atol = 1e-8;
    % Jacobian for the stiff integrators, both solvers
    % 0 : dense numerical, 1 : sparsity pattern (nodal solver; numerical in
    % the spectral solver), 2 : analytic, 3 : batched finite differences
    % (nodal solver)
    d.jacobian = 2;
    % differentiation matrices of the nodal solver imr.solve_nodal
    % 0 : finite difference, 1 : Chebyshev
    d.scheme = 0;

    % initial conditions

    % 0 : start from R0 = 1, 1 : Rayleigh collapse from the maximum radius
    d.collapse = 0;
    % initial bubble wall velocity, Rdot(0)/Uc
    d.Rdotzero = physical.rdotzero;
    % equilibrium radius for pre-stress bubble, Req/R0, see Estrada JMPS 2017
    d.Req = physical.req;

    % output options

    % display progress while code running; note this prints on every
    % right-hand-side evaluation and adds noticeable wall-clock time
    d.progdisplay = 0;

    % acoustic options

    % far-field sound speed, C8/Uc
    d.Cstar = physical.cstar;
    % Tait state equation parameter, GAM/P8
    d.GAMa = physical.gama;
    % state equation parameter
    d.nstate = 7.15;
    % hugoniot slope
    d.hugoniot_s = 1.65;

    % pressure wave options

    % pressure amplitude, pA/P8
    d.ee = physical.ee;
    % frequency, omega*t0
    d.om = physical.om;
    % Gaussian width, TW/t0
    d.tw = physical.tw;
    % delay, DT/t0
    d.dt = physical.dt;
    % power shift for waveform
    d.mn = 0;
    % wave type oscillating bubble, see imr.physics.far_field_pressure
    d.wave_type = 0;

    % stress options

    % Weber number, P8*R0/(2*S)
    d.We = physical.we;
    % Reynolds number, P8*R0/(mu8*Uc)
    d.Re8 = physical.re;
    % Cauchy number, P8/G
    d.Ca = physical.ca;
    % Deborah number, lambda1*Uc/R0
    d.De = physical.de;
    % relaxation time ratio, lambda2/lambda1
    d.LAM = physical.lam;
    % qKV term
    d.alphax = 0.25;

    % viscosity variables

    % non-Newtonian Reynolds number, P8*R0/((mu0-mu8)*Uc), 0 : no contribution
    d.DRe = physical.dre;
    % power factor 1
    d.v_a = 0;
    % power factor 2
    d.v_nc = 0;
    % viscosity time scale, v_lambda/t0
    d.v_lambda_star = physical.v_lambda;
    % viscosity model
    d.nu_model = 0;

    % thermal options

    % ratio of specific heats
    d.kappa = 1.47;
    % thermal Fourier number, Dm/(Uc*R0)
    d.Foh = physical.foh;
    % Brinkman number, Uc^2/(Cp*T8)
    d.Br = physical.br;
    % gas and vapor thermal conductivity coefficients, scaled by K8
    d.alpha_g = physical.alpha_g;
    d.beta_g = physical.beta_g;
    d.alpha_v = physical.alpha_v;
    d.beta_v = physical.beta_v;
    % thermal conduction number, T8*K8/(P8*R0*Uc)
    d.chi = physical.chi;
    % medium to gas thermal conductivity ratio, Km/K8
    d.Km_star = physical.kmstar;

    % mass transfer options

    % mass Fourier number, D0/(Uc*R0)
    d.Fom = physical.fom;
    % latent heat of evaporation, L_heat/Uc^2
    d.L_heat_star = physical.lheat;
    % gas constants of the vapor and non-condensible gas, R/(P8/(rho8*T8))
    d.Rv_star = physical.rv;
    d.Ra_star = physical.ra;

    % pressure variables

    % saturated vapor pressure at the far-field temperature, Pv/P8
    d.Pv_sat = physical.pv;
    % vapor pressure temperature exponent, Pv(T)/P8 = Pv_sat*exp(thetav*(1-1/T))
    d.thetav = physical.thetav;

end
