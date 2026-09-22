% file f_call_params.m
% brief contains function f_call_params

% brief This function sets up the simulation environment. The function can
% receive an input file. The input file must be in the formatted similar to
% the default_case.m file. Otherwise, if no input file is provided, the
% code will read in the default case. The default case can be altered based
% on the input arguments when invoking either the finite difference or
% spectral module.
function [eqns_opts, solve_opts, init_opts, init_stress, tspan_opts, out_opts, ...
    acos_opts, wave_opts, sigma_opts, thermal_opts, mass_opts] = ...
    f_call_params(varargin)

disp('--- Inertial Microcavitation Rheometry forward solver ---');
% check that all inputs are matched
if mod(nargin,2) == 1
    error('input error: unmatched inputs');
end

% checking for the casefile, if any
defaultread = true;
for n = 1:2:nargin
    if strcmpi(varargin{n},'casefile') == 1
        cfname = varargin{n+1};
        fprintf('Case file: using %s\n', cfname);
        try
            run(cfname);
        catch
            error('failed to read the given file');
        end
        defaultread = false;
        break;
    end
end

% otherwise, reading the default casefile
if defaultread
    disp('Case file: Using default case file');
    run('default_case.m');
end

% overrides defaults with options and dimensional inputs %

% load inputs
misscount = 0;
p0set = 0;
dmset = 0;
tflag = 0;
visflag = 0;
for n = 1:2:nargin
    if strcmpi(varargin{n},'p0') == 1
        p0set = 1;
    end
    switch lower(varargin{n})
        % equation options
        case 'radial',      radial = varargin{n+1};
        case 'bubtherm',    bubtherm = varargin{n+1};
        case 'medtherm',    medtherm = varargin{n+1};
        case 'stress',      stress = varargin{n+1};
        case 'eps3',        eps3 = varargin{n+1};
        case 'vapor',       vapor = varargin{n+1};
        case 'masstrans',   masstrans = varargin{n+1};
        
        % solver options
        case 'method',      method = varargin{n+1};
        case 'spectral',    spectral = varargin{n+1};
        case 'divisions',   divisions = varargin{n+1};
        case 'nv',          Nv = varargin{n+1};
        case 'nt',          Nt = varargin{n+1};
        case 'mt',          Mt = varargin{n+1};
        case 'lv',          Lv = varargin{n+1};
        case 'lt',          Lt = varargin{n+1};
        case 'reltol',      rtol = varargin{n+1};
        case 'abstol',      atol = varargin{n+1};
        case 'jacobian',    jacobian = varargin{n+1};
        case 'tfin',        TFin = varargin{n+1};
        tflag = tflag + 1;
        case 'tvector',     TVector = varargin{n+1};
        tflag = tflag + 1;
        TFin = 0;
        
        % initial options
        case 'collapse',    collapse = varargin{n+1};
        case 'r0',          R0 = varargin{n+1};
        case 'u0',          U0 = varargin{n+1};
        case 'req',         Req = varargin{n+1};
        case 'stress0',     Szero = varargin{n+1};
        
        % output options
        case 'dimout',      dimensionalout = varargin{n+1};
        case 'progdisplay', progdisplay = varargin{n+1};
        
        % acoustic options
        case 'rho8',        rho8 = varargin{n+1};
        case 'gam',         GAM = varargin{n+1};
        case 'nstate',      nstate = varargin{n+1};
        case 'p8',          P8 = varargin{n+1};
        case 'c8',          C8 = varargin{n+1};
        case 'hugoniot_s',  hugoniot_s = varargin{n+1};
        
        % pressure waveform options
        case 'pa',          pA = varargin{n+1};
        case 'omega',       omega = varargin{n+1};
        case 'tw',          TW = varargin{n+1};
        case 'dt',          DT = varargin{n+1};
        case 'mn',          mn = varargin{n+1};
        case 'wave_type',   wave_type = varargin{n+1};
        
        % stress options
        case 'mu',          mu8 = varargin{n+1};
        visflag = visflag + 1;
        case 'g',           G = varargin{n+1};
        case 'lambda1',     lambda1 = varargin{n+1};
        case 'lambda2',     lambda2 = varargin{n+1};
        case 'alphax',      alphax = varargin{n+1};
        case 'surft',       S = varargin{n+1};
        
        % non-Newtonian viscosity options
        case 'du',          Dmu         = varargin{n+1};
        case 'mu0',         muo         = varargin{n+1};
        Dmu = muo - mu8;
        case 'v_a',         v_a         = varargin{n+1};
        case 'v_nc',        v_nc        = varargin{n+1};
        case 'v_lambda',    v_lambda    = varargin{n+1};
        case 'nu_model',    nu_model     = varargin{n+1};
        
        % thermal options
        case 't8',          T8 = varargin{n+1};
        case 'kappa',       kappa = varargin{n+1};
        case 'atg',         ATg = varargin{n+1};
        case 'btg',         BTg = varargin{n+1};
        case 'atv',         ATv = varargin{n+1};
        case 'btv',         BTv = varargin{n+1};
        case 'km',          Km = varargin{n+1};
        case 'dm',          Dm = varargin{n+1};
        dmset = 1;
        
        % mass transfer options
        case 'dmass',       D0 = varargin{n+1};
        case 'lheat',       L_heat = varargin{n+1};
        case 'rv',          Rv = varargin{n+1};
        case 'ra',          Ra = varargin{n+1};
        
        % pressure options
        case 'pv',          Pv = varargin{n+1};
        case 'p0',          P0 = varargin{n+1};
        
        % handled in the pre-scan; matched here only, so that it is
        % counted once by the second (nondimensional) switch
        case 'casefile'
        
        otherwise,          misscount = misscount + 1;
        
    end
end

% vapor pressure at the far-field temperature
Pv = vapor*f_pvsat(T8);
% initial gas pressure: this block is the single authoritative computation
% of P0. A user-supplied 'p0' is kept as given. Otherwise P0 follows from
% mechanical equilibrium at Req, compressed from Req to R0 polytropically
% (exponent 3*kappa) for bubtherm = 0 and isothermally (exponent 3) for
% bubtherm = 1, where the thermal PDE resolves the temperature.
if p0set == 0
    if bubtherm == 0
        P0 = (P8 + 2*S/Req - Pv)*(Req/R0)^(3*kappa);
    else
        P0 = (P8 + 2*S/Req - Pv)*(Req/R0)^3;
    end
end
if dmset == 0
    % (m^2/s) thermal diffusivity
    Dm = Km / (rho8*Cp);
end

% loading waveform data
if wave_type < 0
    % waveform files live next to this function in src/
    srcdir = fileparts(mfilename('fullpath'));
    if wave_type == -1
        waveform_dir = fullfile(srcdir,'d_hn.mat');
    elseif wave_type == -2
        waveform_dir = fullfile(srcdir,'d_mn.mat');
    elseif wave_type == -3
        waveform_dir = fullfile(srcdir,'d_ml.mat');
    else
        error('INPUT ERROR: wave_type must be an integer in -3..3');
    end
    tmp = load(waveform_dir,'poly','dpoly');
    wave_poly  = tmp.poly;
    wave_dpoly = tmp.dpoly;
else
    wave_poly = [];
    wave_dpoly = [];
end

check = 1-isnumeric(radial);
if check || radial > 6 || radial <= 0
    error('INPUT ERROR: radial must be an integer in 1..6');
end
check = 1-isnumeric(bubtherm);
if check || bubtherm ~= 0 && bubtherm ~= 1
    error('INPUT ERROR: bubtherm must be 0 or 1');
end
check = 1-isnumeric(medtherm);
if check || medtherm ~= 0 && medtherm ~= 1
    error('INPUT ERROR: medtherm must be 0 or 1');
end
if bubtherm == 0 && medtherm ~= 0
    error('INPUT ERROR: medtherm must be 0 if bubtherm = 0');
end
if medtherm == 1 && bubtherm ~= 1
    error('INPUT ERROR: bubtherm must be 1 if medtherm = 1');
end
check = 1-isnumeric(stress);
if check || stress > 5 || stress < 0
    error('INPUT ERROR: stress must be between 0 to 5');
end
check = 1-isnumeric(vapor);
if check || vapor ~= 0 && vapor ~= 1
    error('INPUT ERROR: vapor must be 0 or 1');
end
check = isnumeric(masstrans);
if check && masstrans ~= 0 && masstrans ~= 1
    error('INPUT ERROR: masstrans must be 0 or 1');
end
if check && masstrans == 1 && vapor == 0
    error('INPUT ERROR: if masstrans is 1 you must have vapor be 1')
end
if masstrans == 1 && bubtherm == 0
    error(['INPUT ERROR: masstrans = 1 requires bubtherm = 1 ' ...
        '(isothermal vapor diffusion is not implemented)']);
end
if nu_model > 2
    error(['INPUT ERROR: nu_model 3-7 are not yet implemented; ' ...
        'supported: 0 (Newtonian), 1 (Carreau), 2 (Carreau-Yasuda)']);
end
if stress >= 3 && medtherm == 1 && spectral == 0
    warning('IMR:dissipation', ['Medium heating for stress>=3 currently ' ...
        'uses the Kelvin-Voigt closed form (inconsistent for ' ...
        'history-dependent stress); see docs.']);
end
check = 1-isnumeric(wave_type);
if check || wave_type > 3 || wave_type < -3
    error('INPUT ERROR: wave_type must be an integer in -3..3');
end
if (tflag > 1)
    error('INPUT ERROR: Only tvector or tfin can be specified, not both');
end
if TVector == 0
    TVector = [0 TFin];
end
check = 1-isnumeric(collapse);
if check || collapse ~= 0 && collapse ~= 1
    error('INPUT ERROR: collapse must be 0 or 1');
end
if collapse == 1 && (bubtherm ~= 1 || medtherm ~= 1 || masstrans ~= 1 || vapor ~= 1)
    error('INPUT ERROR: collapse must have full model with vapor');
end

% intermediate calculated variables

% far-field thermal conductivity (W/(m K))
K8          = 0.5*(ATg*T8+BTg + ATv*T8+BTv);
% dimensional parameter for gas constants
Rnondim     = P8/(rho8*T8);
% characteristic velocity (m/s)
Uc          = sqrt(P8/rho8);

% final non-dimensional variables
Pref        = P8;
% dimensionless infinity pressure
P0_star     = P0/Pref;
% dimensionless vapor pressure
Pv_star     = vapor*f_pvsat(T8)/P8;
% characteristic time (s)
t0          = R0/Uc;

% dimensionless waveform parameters
tvector     = TVector./t0;
% non-dimensional frequency
om          = omega*t0;
ee          = pA/Pref;
tw          = TW/t0;
dt          = DT/t0;
% acoustic properties

% bulk liquid stiffness
GAMa        = GAM/P8;
% speed of sound
Cstar       = C8/Uc;
% thermal properties
chi         = T8*K8/(P8*R0*Uc);
iota        = Km/(K8*Lt);
Foh         = Dm/(Uc*R0);
alpha_g     = ATg*T8/K8;
beta_g      = BTg/K8;
alpha_v     = ATv*T8/K8;
beta_v      = BTv/K8;
Br          = Uc^2/(Cp*T8);
% mass diffusion
Fom         = D0/(Uc*R0);
% mass of vapor
mv0         = vapor*Pv*(4/3*pi*R0^3)/Rv/T8;
% mass of non-condensible gas
ma0         = P0*(4/3*pi*R0^3)/Ra/T8;
Mnondim     = rho8*(4/3*pi*R0^3);
mv0         = mv0/Mnondim;
ma0         = ma0/Mnondim;
Rv_star     = Rv/Rnondim;
Ra_star     = Ra/Rnondim;

% non-dimensional latent heat
L_heat_star = L_heat/(Uc)^2;

% viscoelastic properties

% Cauchy number
Ca      = Pref/G;
% Reynolds number
Re8     = Pref*R0/(mu8*Uc);
if Dmu ~= 0
    DRe = Pref*R0/(Dmu*Uc);
else
    DRe = 0;
end
if DRe==0
    iDRe = 0;
else
    iDRe = 1/DRe;
end

% Weber number
We = Pref*R0/(2*S);
% relaxation time
v_lambda_star = v_lambda/t0;
% Weissenberg number
if lambda1 > 0
    LAM = lambda2/lambda1;
else
    % no relaxation time; LAM unused for stress <= 2
    LAM = 0;
end
% Deborah number
De      = lambda1*Uc/R0;
% dimensionless initial conditions
Rzero   = 1;
Req_zero = Req/R0;
Rdotzero= U0/Uc;

% overwrite defaults with nondimensional inputs
if isempty(varargin) == 0
    for n = 1:2:nargin
        switch lower(varargin{n})
            % dimensionless state variables
            case 'cstar',   Cstar = varargin{n+1};
            case 'gama',    GAMa = varargin{n+1};
            % dimensionless waveform parameters
            case 'om',      om = varargin{n+1};
            case 'ee',      ee = varargin{n+1};
            case 'twx',     tw = varargin{n+1};
            case 'dtx',     dt = varargin{n+1};
            % dimensionless numbers
            case 'we',      We = varargin{n+1};
            case 're',      Re8 = varargin{n+1};
            case 'dre',     DRe = varargin{n+1};
            case 'ca',      Ca = varargin{n+1};
            case 'lam',     LAM = varargin{n+1};
            case 'de',      De = varargin{n+1};
            case 'foh',     Foh = varargin{n+1};
            case 'br',      Br = varargin{n+1};
            case 'fom',     Fom = varargin{n+1};
            % dimensionless thermal quantities
            case 'alpha_g', alpha_g = varargin{n+1};
            case 'beta_g',  beta_g = varargin{n+1};
            case 'alpha_v', alpha_v = varargin{n+1};
            case 'beta_v',  beta_v = varargin{n+1};
            case 'chi',     chi = varargin{n+1};
            case 'iota',    iota = varargin{n+1};
            % dimensionless initial conditions
            case 'rzero',   Rzero = varargin{n+1};
            case 'rdotzero',Rdotzero = varargin{n+1};
            case 'p0star',  P0_star = varargin{n+1};
            otherwise
            misscount = misscount + 1;
        end
    end
end
% check that all inputs were accounted for
if misscount ~= nargin/2
    error(['INPUT ERROR: ' num2str(misscount-nargin/2) ' unrecognized input(s)']);
end

% final setting adjustments

% 1 : N-H, 2: qN-H, 3: linear Maxwell, Jeffreys, Zener, 5: UCM or OldB, 6: PTT, 7: Giesekus
if stress == 0 || stress == 1 || stress == 2 || stress == 3 || stress == 4
    spectral = 0;
    if stress == 3 || stress == 4
        Nv = 1;
    else
        Nv = 0;
    end
elseif stress == 5
    Nv = (Nv - 1)*(spectral == 1) + 1;
    Ca = -1;
elseif stress == 6
    Nv = (Nv - 1)*(spectral == 1) + 1;
    Ca = -1;
    spectral = 1;
elseif stress == 7
    Nv = (Nv - 1)*(spectral == 1) + 1;
    Ca = -1;
    spectral = 1;
end

if Ca == -1
    Ca = Inf;
end

if stress == 1 || stress == 2
    JdotA = 4/Re8;
elseif stress == 3 || stress == 4
    % the stress integral carries -4*LAM/Re8*Rdot/R, so the Rddot
    % coefficient moved to the Keller-Miksis denominator carries LAM too
    JdotA = 4*LAM/Re8;
elseif stress == 5 || stress == 6
    JdotA = 4*LAM/Re8;
else
    JdotA = 0;
end
if spectral == 1
    JdotA = 0;
end

% Keller-Miksis equation
% elseif nhkv_pld==1
%     %JdotA = 4/Re8*(2^alpha_g+1)/3*(abs(Rdot)/R)^(alpha_g-1);
%     JdotA = 4/Re8/3*(2^alpha_g+1)*sign(Rdot)*(abs(Rdot)/R)^(alpha_g)*R^2/Rdot^2;
%     if isnan(SdotA)
%         SdotA=4/Re8;
%     end
% end

if stress == 0 || stress == 1 || stress == 2 || stress == 3 || stress == 4
    zeNO = 0;
else
    zeNO = 1;
end

if (stress == 3 || stress == 4 || stress == 5) && De == 0
    error('INPUT ERROR: De cannot equal zero for stress = 3, 4, or 5');
end

% inertial Rayleigh collapse out of equilibrium initial conditions
opts = optimset('display','off');
if collapse
    % calculate the equilibrium radii ratio for initial stress state
    fun = @(x) Pv_star*(1+(ma0/x)*(Ra_star/Rv_star))-1-(1/We)*(Pv_star/(Rv_star*x))^(1/3);
    MTotal0 = ma0 + mv0;
    MVE = fzero(fun,MTotal0,opts);
    while (isnan(MVE))
        MTotal0 = MTotal0/1.11;
        MVE = fzero(fun,MTotal0,opts);
    end
end
Pb_star = P0_star + Pv_star;
% initial concentration, mass air / mass vapor
theta = Rv_star/Ra_star*(Pb_star/Pv_star - 1);
kv0 = 1/(1+theta);
if collapse
    Req_zero = (Rv_star*MVE/Pv_star)^(1/3);
    fprintf('New equilibrium radius, Req = %.8f\n',Req_zero);
    Rdotzero = 0;
    %-(1-Pb_star)/(Cstar);
end

% initial stress field for bubble collapse
if stress < 3
    Szero = [];
elseif stress == 3 || stress == 4
    if collapse
        if ~exist('Szero','var')
            [Szero] = f_init_stress(Req_zero, Re8, Ca, De, We, Cstar, Pv_star);
        end
    else
        Szero = 0;
    end
elseif stress == 5
    % TODO initial max stress for UCM and Oldroyd-B
    Szero = zeros((Nv - 1)*(spectral == 1) + 2,1);
end

% out parameters

% equation settings
eqns_opts = [radial bubtherm medtherm stress eps3 masstrans];
% solver options
solve_opts = [method spectral divisions Nv Nt Mt Lv Lt rtol atol jacobian];
% dimensionless initial conditions
init_opts = [Rzero Rdotzero Pb_star P8 T8 Pv_star Req_zero alphax];
% dimensionaless initial stress
init_stress = Szero;
% time span options
tspan_opts = tvector;
% output options
out_opts = [dimensionalout progdisplay t0 R0 Uc];

% physical parameters

% acoustic parameters
acos_opts = [Cstar GAMa kappa nstate hugoniot_s];
% dimensionless waveform parameters
wave_opts = {[om ee tw dt mn wave_type], wave_poly, wave_dpoly};
% dimensionless viscoelastic
sigma_opts = [We Re8 DRe v_a v_nc Ca LAM De JdotA nu_model v_lambda_star zeNO iDRe];
% dimensionless thermal
thermal_opts = [Foh Br alpha_g beta_g alpha_v beta_v chi iota];
% dimensionaless mass transfer
mass_opts = [Fom kv0 Rv_star Ra_star L_heat_star mv0 ma0];

end
