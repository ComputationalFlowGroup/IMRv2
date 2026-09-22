% file default_case.m
% brief contains the default (dimensionless) simulation case

% brief All quantities are dimensionless. Lengths are scaled by the initial
% radius R0, pressures by the far-field pressure P8, velocities by
% Uc = sqrt(P8/rho8), time by t0 = R0/Uc, and temperatures by the far-field
% temperature T8. The values below correspond to a 50 um bubble in water
% at 1 atm and 298.15 K; toolchain/f_dim2nondim.m converts dimensional
% (SI) properties into these numbers.

% equation options

% 1 : RP, 2 : K-M P, 3 : K-M E, 4 : Gil, 5 : K-M MG, 6 : Gil MG
radial          = 1;
% 0 : polytropic assumption, 1: thermal PDE in bubble
bubtherm        = 0;
% 0 : cold fluid, 1: warm fluid assumption
medtherm        = 0;
% 0 : none, 1 : NHKV, 2 : qNHKV, 3 : linear Maxwell, Jeffreys, Zener
% (neo-Hookean), 4 : same with quadratic neo-Hookean, 5 : UCM or OldB
stress          = 0;
% this value must be [0, 0.5]
eps3            = 0;
% 0 : ignore vapor pressure, 1 : vapor pressure
vapor           = 0;
% mass transfer, default is no mass transfer
masstrans       = 0;

% solver options

% final time
TFin            = 3.903439778613622;
% time span vector
TVector         = linspace(0,TFin,128);
% ode45 setting for the time stepper
method          = 23;
% force spectral collocation solution
spectral        = 0;
% minimum number of timesteps
divisions       = 0;
% number of points in bubble, thermal PDE
Nt              = 25;
% number of points outside of bubble, thermal PDE
Mt              = 25;
% number of points outside of bubble, viscoelastic stress PDE
Nv              = 150;
% characteristic length for grid stretching, leave at 3
Lv              = 3;
% characteristic length for grid stretching, leave at 3
Lt              = 2;
% relative tolerance of the time integrator
rtol            = 1e-8;
% absolute tolerance of the time integrator
atol            = 1e-8;
% Jacobian for the stiff integrators (finite difference solver only)
% 0 : dense numerical, 1 : sparsity pattern, 2 : analytic (falls back to
% the sparsity pattern for configurations it does not cover)
jacobian        = 2;

% initial conditions

collapse        = 0;
% initial bubble wall velocity, Rdot(0)/Uc
Rdotzero        = 0;
% equilibrium radius for pre-stress bubble, Req/R0, see Estrada JMPS 2017
Req             = 0.2;

% output options

% display progress while code running; note this prints on every
% right-hand-side evaluation and adds noticeable wall-clock time
progdisplay     = 0;

% acoustic options

% far-field sound speed, C8/Uc
Cstar           = 152.071002414908;
% Tait state equation parameter, GAM/P8
GAMa            = 3009.257340241796;
% state equation parameter
nstate          = 7.15;
% hugoniot slope
hugoniot_s      = 1.65;

% pressure wave options

% pressure amplitude, pA/P8
ee              = 0;
% frequency, omega*t0
om              = 0;
% Gaussian width, TW/t0
tw              = 0;
% delay, DT/t0
dt              = 0;
% power shift for waveform
mn              = 0;
% wave type oscillating bubble, see f_pinfinity
wave_type       = 0;

% stress options

% Weber number, P8*R0/(2*S)
We              = 36.1875;
% Reynolds number, P8*R0/(mu8*Uc)
Re8             = 623.3655014295973;
% Cauchy number, P8/G
Ca              = 101.325;
% Deborah number, lambda1*Uc/R0
De              = 0.01951719889306811;
% relaxation time ratio, lambda2/lambda1
LAM             = 0.1;
% qKV term
alphax          = 0.25;

% viscosity variables

% non-Newtonian Reynolds number, P8*R0/((mu0-mu8)*Uc), 0 : no contribution
DRe             = 0;
% power factor 1
v_a             = 0;
% power factor 2
v_nc            = 0;
% viscosity time scale, v_lambda/t0
v_lambda_star   = 0;
% viscosity model
nu_model        = 0;

% thermal options

% ratio of specific heats
kappa           = 1.47;
% thermal Fourier number, Dm/(Uc*R0)
Foh             = 0.0002533864560818251;
% Brinkman number, Uc^2/(Cp*T8)
Br              = 7.639413243495116e-05;
% gas and vapor thermal conductivity coefficients, scaled by K8
alpha_g         = 0.5761007932661034;
beta_g          = 0.4263395891806357;
alpha_v         = 0.3600629957913147;
beta_v          = 0.6374966217619463;
% thermal conduction number, T8*K8/(P8*R0*Uc)
chi             = 0.1647900432365406;
% medium to gas thermal conductivity ratio, Km/K8
Km_star         = 20.12762009007293;

% mass transfer options

% mass Fourier number, D0/(Uc*R0)
Fom             = 0.04959728111106163;
% latent heat of evaporation, L_heat/Uc^2
L_heat_star     = 23781.93575129533;
% gas constants of the vapor and non-condensible gas, R/(P8/(rho8*T8))
Rv_star         = 1444.950104049612;
Ra_star         = 898.6805465194676;

% pressure variables

% saturated vapor pressure at the far-field temperature, Pv/P8
Pv_sat          = 0.03076018438379325;
% vapor pressure temperature exponent, Pv(T)/P8 = Pv_sat*exp(thetav*(1-1/T))
thetav          = 17.44088546033876;
