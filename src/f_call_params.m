% file f_call_params.m
% brief contains function f_call_params

% brief This function sets up the simulation environment. All inputs are
% dimensionless (see default_case.m for the scalings); dimensional values
% must be converted beforehand, e.g. with toolchain/f_dim2nondim.m. The
% function can receive an input file formatted similar to default_case.m,
% otherwise the default case is read. The case can be altered with
% name/value pairs when invoking either the finite difference or spectral
% module. The inputs are checked by f_checkers. The first argument names
% the calling solver, 'nodal' or 'spectral'.
function [eqns_opts, solve_opts, init_opts, init_stress, tspan_opts, out_opts, ...
    acos_opts, wave_opts, sigma_opts, thermal_opts, mass_opts] = ...
    f_call_params(solver, varargin)

disp('--- Inertial Microcavitation Rheometry forward solver ---');
% check that all inputs are matched
if mod(numel(varargin),2) == 1
    error('input error: unmatched inputs');
end

% checking for the casefile, if any
defaultread = true;
for n = 1:2:numel(varargin)
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
% case files written before the nodal solver have no scheme setting
if ~exist('scheme','var')
    scheme = 0;
end

% overrides defaults with the dimensionless inputs

p0set = 0;
stress0set = 0;
tvecset = 0;
tfinset = 0;
wave_poly = [];
wave_dpoly = [];
unknown = {};
for n = 1:2:numel(varargin)
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
        case 'scheme',      scheme = varargin{n+1};
        case 'tfin',        TFin = varargin{n+1};
        tfinset = 1;
        case 'tvector',     TVector = varargin{n+1};
        tvecset = 1;
        
        % initial options
        case 'collapse',    collapse = varargin{n+1};
        case 'req',         Req = varargin{n+1};
        case 'rdotzero',    Rdotzero = varargin{n+1};
        case 'p0star',      P0_star = varargin{n+1};
        p0set = 1;
        case 'stress0',     Szero = varargin{n+1};
        stress0set = 1;
        
        % output options
        case 'progdisplay', progdisplay = varargin{n+1};
        
        % acoustic options
        case 'cstar',       Cstar = varargin{n+1};
        case 'gama',        GAMa = varargin{n+1};
        case 'nstate',      nstate = varargin{n+1};
        case 'hugoniot_s',  hugoniot_s = varargin{n+1};
        
        % pressure waveform options
        case 'ee',          ee = varargin{n+1};
        case 'om',          om = varargin{n+1};
        case 'tw',          tw = varargin{n+1};
        case 'dt',          dt = varargin{n+1};
        case 'mn',          mn = varargin{n+1};
        case 'wave_type',   wave_type = varargin{n+1};
        
        % stress options
        case 'we',          We = varargin{n+1};
        case 're',          Re8 = varargin{n+1};
        case 'ca',          Ca = varargin{n+1};
        case 'de',          De = varargin{n+1};
        case 'lam',         LAM = varargin{n+1};
        case 'alphax',      alphax = varargin{n+1};
        
        % non-Newtonian viscosity options
        case 'dre',         DRe = varargin{n+1};
        case 'v_a',         v_a = varargin{n+1};
        case 'v_nc',        v_nc = varargin{n+1};
        case 'v_lambda',    v_lambda_star = varargin{n+1};
        case 'nu_model',    nu_model = varargin{n+1};
        
        % thermal options
        case 'kappa',       kappa = varargin{n+1};
        case 'foh',         Foh = varargin{n+1};
        case 'br',          Br = varargin{n+1};
        case 'alpha_g',     alpha_g = varargin{n+1};
        case 'beta_g',      beta_g = varargin{n+1};
        case 'alpha_v',     alpha_v = varargin{n+1};
        case 'beta_v',      beta_v = varargin{n+1};
        case 'chi',         chi = varargin{n+1};
        case 'kmstar',      Km_star = varargin{n+1};
        
        % mass transfer options
        case 'fom',         Fom = varargin{n+1};
        case 'lheat',       L_heat_star = varargin{n+1};
        case 'rv',          Rv_star = varargin{n+1};
        case 'ra',          Ra_star = varargin{n+1};
        
        % vapor pressure options
        case 'pv',          Pv_sat = varargin{n+1};
        case 'thetav',      thetav = varargin{n+1};
        
        % handled in the pre-scan
        case 'casefile'
        
        otherwise,          unknown{end+1} = varargin{n}; %#ok<AGROW>
        
    end
end
if ~isempty(unknown)
    error(['INPUT ERROR: unrecognized input(s): ' strjoin(unknown,', ') ...
        '. IMR takes dimensionless inputs only; toolchain/f_dim2nondim.m ' ...
        'converts dimensional values.']);
end
if tfinset && ~tvecset
    TVector = [0 TFin];
end

% loading the (dimensionless) experimental waveform data
if wave_type < 0
    % waveform files live next to this function in src/
    srcdir = fileparts(mfilename('fullpath'));
    if wave_type == -1
        waveform_dir = fullfile(srcdir,'d_hn.mat');
    elseif wave_type == -2
        waveform_dir = fullfile(srcdir,'d_mn.mat');
    elseif wave_type == -3
        waveform_dir = fullfile(srcdir,'d_ml.mat');
    end
    if exist('waveform_dir','var')
        tmp = load(waveform_dir,'poly','dpoly');
        wave_poly  = tmp.poly;
        wave_dpoly = tmp.dpoly;
    end
end

% derived dimensionless quantities

% vapor pressure at the far-field temperature
Pv_star = vapor*Pv_sat;
% initial gas pressure: a user-supplied 'p0star' is kept as given.
% Otherwise P0 follows from mechanical equilibrium at Req, compressed from
% Req to R0 polytropically (exponent 3*kappa) for bubtherm = 0 and
% isothermally (exponent 3) for bubtherm = 1, where the thermal PDE
% resolves the temperature
if p0set == 0
    if bubtherm == 0
        P0_star = (1 + 1/(We*Req) - Pv_star)*Req^(3*kappa);
    else
        P0_star = (1 + 1/(We*Req) - Pv_star)*Req^3;
    end
end
Pb_star = P0_star + Pv_star;
% medium to gas conductivity ratio over the grid stretching length
iota = Km_star/Lt;
% non-Newtonian viscosity contribution
if DRe == 0
    iDRe = 0;
else
    iDRe = 1/DRe;
end
% initial vapor concentration, mass of vapor and non-condensible gas
theta = Rv_star/Ra_star*(Pb_star/Pv_star - 1);
kv0 = 1/(1+theta);
mv0 = Pv_star/Rv_star;
ma0 = P0_star/Ra_star;
% dimensionless initial conditions
Rzero = 1;
Req_zero = Req;
% dimensionless time
tvector = TVector;

% compatibility of the conditions and dimensionless numbers; settings
% tied to the stress model (spectral, Nv, Ca) come back adjusted
chk = f_checkers(struct('solver',solver, ...
    'radial',radial, 'bubtherm',bubtherm, 'medtherm',medtherm, ...
    'stress',stress, 'eps3',eps3, 'vapor',vapor, 'masstrans',masstrans, ...
    'collapse',collapse, 'method',method, 'spectral',spectral, ...
    'divisions',divisions, 'Nv',Nv, 'Nt',Nt, 'Mt',Mt, 'Lv',Lv, 'Lt',Lt, ...
    'rtol',rtol, 'atol',atol, 'jacobian',jacobian, 'scheme',scheme, ...
    'tflag',tfinset + tvecset, 'tvector',tvector, 'progdisplay',progdisplay, ...
    'Req',Req_zero, 'Rdotzero',Rdotzero, 'P0_star',P0_star, ...
    'Pv_star',Pv_star, 'Pv_sat',Pv_sat, 'thetav',thetav, ...
    'Cstar',Cstar, 'GAMa',GAMa, 'nstate',nstate, 'hugoniot_s',hugoniot_s, ...
    'ee',ee, 'om',om, 'tw',tw, 'dt',dt, 'mn',mn, 'wave_type',wave_type, ...
    'wave_poly',{wave_poly}, 'wave_dpoly',{wave_dpoly}, ...
    'We',We, 'Re8',Re8, 'Ca',Ca, 'De',De, 'LAM',LAM, 'alphax',alphax, ...
    'DRe',DRe, 'v_a',v_a, 'v_nc',v_nc, 'v_lambda_star',v_lambda_star, ...
    'nu_model',nu_model, 'kappa',kappa, 'Foh',Foh, 'Br',Br, ...
    'alpha_g',alpha_g, 'beta_g',beta_g, 'alpha_v',alpha_v, ...
    'beta_v',beta_v, 'chi',chi, 'Km_star',Km_star, 'Fom',Fom, ...
    'L_heat_star',L_heat_star, 'Rv_star',Rv_star, 'Ra_star',Ra_star));
spectral = chk.spectral;
Nv = chk.Nv;
Ca = chk.Ca;

% derived stress-model quantities

if stress == 1 || stress == 2
    JdotA = 4/Re8;
elseif stress == 3 || stress == 4 || stress == 5
    % the stress integral carries -4*LAM/Re8*Rdot/R, so the Rddot
    % coefficient moved to the Keller-Miksis denominator carries LAM too
    JdotA = 4*LAM/Re8;
else
    JdotA = 0;
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
    Req_zero = (Rv_star*MVE/Pv_star)^(1/3);
    fprintf('New equilibrium radius, Req = %.8f\n',Req_zero);
    Rdotzero = 0;
end

% initial stress state: none for stress 0-2, one state for stress 3-4,
% and two states (Z1, Z2) or Nv Chebyshev coefficients of each of tau_rr
% and tau_thetatheta (spectral stress representation) for stress 5
if stress < 3
    nstress = 0;
elseif stress == 3 || stress == 4
    nstress = 1;
else
    nstress = 2*Nv;
end
if stress0set
    % a user-supplied 'stress0' is kept as given
    Szero = Szero(:);
    if numel(Szero) ~= nstress
        error('INPUT ERROR: stress0 must have %d element(s) for stress = %d', ...
            nstress, stress);
    end
elseif (stress == 3 || stress == 4) && collapse
    % stress at the maximum radius of the collapse
    Szero = f_init_stress(Req_zero, Re8, Ca, De, We, Cstar, Pv_star);
else
    % TODO initial max stress for UCM and Oldroyd-B
    Szero = zeros(nstress,1);
end

% out parameters

% equation settings
eqns_opts = [radial bubtherm medtherm stress eps3 masstrans];
% solver options
solve_opts = [method spectral divisions Nv Nt Mt Lv Lt rtol atol jacobian ...
    scheme];
% dimensionless initial conditions
init_opts = [Rzero Rdotzero Pb_star Pv_star Req_zero alphax];
% dimensionless initial stress
init_stress = Szero;
% time span options
tspan_opts = tvector;
% output options
out_opts = progdisplay;

% physical parameters

% acoustic parameters
acos_opts = [Cstar GAMa kappa nstate hugoniot_s];
% dimensionless waveform parameters
wave_opts = {[om ee tw dt mn wave_type], wave_poly, wave_dpoly};
% dimensionless viscoelastic
sigma_opts = [We Re8 DRe v_a v_nc Ca LAM De JdotA nu_model v_lambda_star iDRe];
% dimensionless thermal
thermal_opts = [Foh Br alpha_g beta_g alpha_v beta_v chi iota];
% dimensionless mass transfer
mass_opts = [Fom kv0 Rv_star Ra_star L_heat_star mv0 ma0 thetav];

end
