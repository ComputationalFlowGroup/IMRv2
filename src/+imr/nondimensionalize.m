% file nondimensionalize.m
% brief contains function imr.nondimensionalize

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This side function converts dimensional (SI) physical properties
% into the dimensionless inputs taken by the IMR solvers (imr.solve_nodal,
% imr.solve_spectral). It is not called by the solvers: the solvers accept only
% dimensionless inputs, and this helper is provided for users who start
% from dimensional values.
%
% Usage:
%   [varin, ref] = imr.nondimensionalize('r0',50e-6,'req',5e-6,'mu',1e-2,'g',1e3);
%   [t, R] = imr.solve_nodal(varin{:}, 'radial', 2, 'Nt', 50);
%   % t*ref.t0 is the time in s, R*ref.R0 the radius in m
%
% Dimensional inputs (SI units) and their defaults:
%   r0 (m) 50e-6, req (m) r0/5, u0 (m/s) 0, tfin (s), tvector (s)
%   rho8 (kg/m^3) 1064, p8 (Pa) 101325, c8 (m/s) 1484, gam (Pa) 3049.13e5
%   surft (N/m) 0.07, g (Pa) 1e3, mu (Pa s) 8.3283e-4, mu0 (Pa s) mu,
%   du (Pa s) mu0 - mu, lambda1 (s) 1e-7, lambda2 (s) 1e-8, v_lambda (s) 0
%   t8 (K) 298.15, atg (W/m-K^2) 5.28e-5, btg (W/m-K) 1.165e-2,
%   atv (W/m-K^2) 3.30e-5, btv (W/m-K) 1.742e-2, km (W/m-K) 0.55,
%   cp (J/kg-K) 4.181e3, dm (m^2/s) km/(rho8*cp), dmass (m^2/s) 24.2e-6,
%   lheat (J/kg) 2264.76e3, rv (J/kg-K) Ru/mwv, ra (J/kg-K) Ru/mwg,
%   pv (Pa) water saturation at t8, p0 (Pa) computed by the solver if not given
%   pa (Pa) 0, omega (rad/s) 0, tw (s) 0, dt (s) 0
% Any other name/value pair (model flags, grid sizes, kappa, alphax, mn,
% nstate, hugoniot_s, v_a, v_nc, ...) is passed through unchanged.
%
% Outputs:
%   varin   cell array of dimensionless name/value pairs for the solvers
%   ref     reference scales: R0 (m), Uc (m/s), t0 (s), P8 (Pa), T8 (K),
%           rho8 (kg/m^3), to restore dimensions of the solver outputs
function [varin, ref] = nondimensionalize(varargin)

    if mod(nargin,2) == 1
        error('imr.nondimensionalize: unmatched inputs');
    end
    [varargin, fluid, saturation] = f_prepare_fluid(varargin);
    water = imr.fluids.water();

    % dimensional defaults (SI)
    R0 = 50e-6;
    Req = [];
    U0 = 0;
    TFin = [];
    TVector = [];
    rho8 = water.liquid.rho;
    P8 = 101325;
    C8 = water.liquid.c;
    GAM = water.liquid.tait_pressure;
    S = water.transport.surft;
    G = 1e3;
    mu8 = water.transport.mu;
    muo = [];
    Dmu = [];
    lambda1 = 1e-7;
    lambda2 = 1e-8;
    v_lambda = 0;
    T8 = 298.15;
    [ATg, BTg] = f_gas_conductivity();
    ATv = water.transport.atv;
    BTv = water.transport.btv;
    Km = water.transport.km;
    Cp = water.liquid.Cp;
    Dm = [];
    D0 = water.transport.dmass;
    L_heat = water.caloric.latent_heat;
    % (J/mol-K) universal gas constant, (kg/mol) molecular weights
    Ru = 8.3144598;
    Rv = water.vapor.Rs;
    Ra = Ru/28.966E-3;
    Pv = [];
    P0 = [];
    pA = 0;
    omega = 0;
    TW = 0;
    DT = 0;

    % read the dimensional inputs, pass everything else through
    passthrough = {};
    for n = 1:2:numel(varargin)
        val = varargin{n+1};
        switch lower(varargin{n})
            case 'r0',          R0 = val;
            case 'req',         Req = val;
            case 'u0',          U0 = val;
            case 'tfin',        TFin = val;
            case 'tvector',     TVector = val;
            case 'rho8',        rho8 = val;
            case 'p8',          P8 = val;
            case 'c8',          C8 = val;
            case 'gam',         GAM = val;
            case 'surft',       S = val;
            case 'g',           G = val;
            case 'mu',          mu8 = val;
            case 'mu0',         muo = val;
            case 'du',          Dmu = val;
            case 'lambda1',     lambda1 = val;
            case 'lambda2',     lambda2 = val;
            case 'v_lambda',    v_lambda = val;
            case 't8',          T8 = val;
            case 'atg',         ATg = val;
            case 'btg',         BTg = val;
            case 'atv',         ATv = val;
            case 'btv',         BTv = val;
            case 'km',          Km = val;
            case 'cp',          Cp = val;
            case 'dm',          Dm = val;
            case 'dmass',       D0 = val;
            case 'lheat',       L_heat = val;
            case 'rv',          Rv = val;
            case 'ra',          Ra = val;
            case 'pv',          Pv = val;
            case 'p0',          P0 = val;
            case 'pa',          pA = val;
            case 'omega',       omega = val;
            case 'tw',          TW = val;
            case 'dt',          DT = val;
            otherwise
            passthrough = [passthrough, varargin(n:n+1)]; %#ok<AGROW>
        end
    end

    % dependent defaults
    if isempty(Req)
        Req = R0/5;
    end
    if isempty(muo)
        muo = mu8;
    end
    if isempty(Dmu)
        Dmu = muo - mu8;
    end
    if isempty(Dm)
        Dm = Km/(rho8*Cp);
    end
    if isempty(fluid)
        saturation = imr.fluid_properties(water,[],T8);
    end
    if isempty(Pv)
        Pv = saturation.pressure;
    end

    % reference scales
    Uc = sqrt(P8/rho8);
    t0 = R0/Uc;
    Rnondim = P8/(rho8*T8);
    K8 = 0.5*(ATg*T8 + BTg + ATv*T8 + BTv);

    % dimensionless numbers
    if Dmu ~= 0
        DRe = P8*R0/(Dmu*Uc);
    else
        DRe = 0;
    end
    if lambda1 > 0
        LAM = lambda2/lambda1;
    else
        LAM = 0;
    end

    varin = { ...
        'req',      Req/R0, ...
        'rdotzero', U0/Uc, ...
        'cstar',    C8/Uc, ...
        'gama',     GAM/P8, ...
        'we',       P8*R0/(2*S), ...
        're',       P8*R0/(mu8*Uc), ...
        'dre',      DRe, ...
        'ca',       P8/G, ...
        'de',       lambda1*Uc/R0, ...
        'lam',      LAM, ...
        'v_lambda', v_lambda/t0, ...
        'foh',      Dm/(Uc*R0), ...
        'br',       Uc^2/(Cp*T8), ...
        'alpha_g',  ATg*T8/K8, ...
        'beta_g',   BTg/K8, ...
        'alpha_v',  ATv*T8/K8, ...
        'beta_v',   BTv/K8, ...
        'chi',      T8*K8/(P8*R0*Uc), ...
        'kmstar',   Km/K8, ...
        'fom',      D0/(Uc*R0), ...
        'lheat',    L_heat/Uc^2, ...
        'rv',       Rv/Rnondim, ...
        'ra',       Ra/Rnondim, ...
        'pv',       Pv/P8, ...
        'thetav',   -saturation.c2/T8, ...
        'om',       omega*t0, ...
        'ee',       pA/P8, ...
        'tw',       TW/t0, ...
        'dt',       DT/t0};
    if ~isempty(P0)
        varin = [varin, {'p0star', P0/P8}];
    end
    if ~isempty(TVector)
        varin = [varin, {'tvector', TVector/t0}];
    end
    if ~isempty(TFin)
        varin = [varin, {'tfin', TFin/t0}];
    end
    varin = [varin, passthrough];

    ref = struct('R0',R0,'Uc',Uc,'t0',t0,'P8',P8,'T8',T8,'rho8',rho8);
    if ~isempty(fluid)
        ref.fluid = fluid;
    end

end

% Resolve a data profile once, before nondimensionalization. All dispatch is
% by model/capability, never by material name. No profile enters the ODE.
function [args, fluid, saturation] = f_prepare_fluid(args)
    fluid = [];
    saturation = [];
    names = cell(size(args(1:2:end)));
    for k = 1:numel(names)
        name = args{2*k-1};
        if ~(ischar(name) && isrow(name)) && ~(isstring(name) && isscalar(name))
            error('IMR:fluidInputs', 'Input names must be text scalars.');
        end
        names{k} = lower(char(name));
    end
    index = find(strcmp(names,'fluid'),1,'last');
    if isempty(index)
        return
    end
    p = args{2*index};
    if isa(p,'function_handle')
        p = p();
    end
    if ~(isstruct(p) && isscalar(p))
        error('IMR:fluidProfile', 'fluid must be a profile struct or function handle.');
    end
    allowed = {'t8','p8','r0','req','u0','tfin','tvector','surft','g', ...
        'mu','mu0','du','lambda1','lambda2','v_lambda','atg','btg', ...
        'atv','btv','km','cp','dm','dmass','lheat','rv','ra','pv','p0', ...
        'rho8','c8','gam','pa','omega','tw','dt','kappa','alphax', ...
        'mn','nu_model','v_a','v_nc','nstate','hugoniot_s','radial','stress', ...
        'bubtherm','medtherm','masstrans','collapse','vapor','method', ...
        'spectral','scheme','divisions','nt','mt','nv','lv','lt', ...
        'reltol','abstol','jacobian','progdisplay','stress0','wave_type','eps3'};
    % Only measured/reference transport belongs to a material profile.
    % Gas conductivity, stress choices and forcing remain simulation inputs.
    opts = p.transport;
    [opts.atg, opts.btg] = f_gas_conductivity();
    if isfield(p.liquid,'tait_pressure')
        opts.gam = p.liquid.tait_pressure;
    end
    for k = 1:numel(names)
        if strcmp(names{k},'fluid')
            continue
        end
        if ~any(strcmp(names{k},allowed)) || any(strcmp(names{k},p.derived_inputs))
            error('IMR:fluidOverride', ...
                'Unsupported fluid input %s. Use SI inputs consistent with the profile.', names{k});
        end
        opts.(names{k}) = args{2*k};
    end
    f_required(p.required,names);
    stress = f_option(opts,'stress',0);
    validateattributes(stress, {'numeric'}, {'scalar','real','integer','>=',0,'<=',5});
    if stress > 0
        f_required(p.required_viscous,names);
    end
    % Keep stress, elasticity, relaxation, forcing and bubble gas choices
    % outside material data. If the ordinary radial default is unsupported,
    % choose the first supported branch unless the caller selected one.
    radial_default = 1;
    if ~any(radial_default == p.capabilities.radial)
        radial_default = p.capabilities.radial(1);
        if ~isfield(opts,'radial')
            opts.radial = radial_default;
        end
    end
    radial = f_option(opts,'radial',radial_default);
    validateattributes(radial, {'numeric'}, {'scalar','real','integer'});
    if ~any(radial == p.capabilities.radial)
        error('IMR:fluidRadial', 'The profile supports radial = %s.', mat2str(p.capabilities.radial));
    end
    flags = {'bubtherm','medtherm','vapor','masstrans','collapse'};
    selected = zeros(size(flags));
    for k = 1:numel(flags)
        val = f_option(opts,flags{k},0);
        validateattributes(val, {'numeric'}, {'scalar','real','integer','>=',0,'<=',1});
        selected(k) = val;
    end
    unsupported_mass = (selected(4) && ~p.capabilities.mass_transfer) || ...
        (selected(5) && ~p.capabilities.collapse);
    if unsupported_mass
        error('IMR:fluidMassTransfer', 'This profile does not support the selected mass transfer/collapse model.');
    end
    if selected(2) && ~selected(1)
        error('IMR:fluidThermal', 'medtherm = 1 requires bubtherm = 1.');
    end
    if selected(1)
        f_required(p.required_thermal,names);
    end
    if selected(2)
        f_required(p.required_medium,names);
    end
    % Numerical scales for inactive terms with no material data. These are
    % not estimates of physical transport: active models were required to
    % provide the missing values above. Never put these in material files.
    inactive = struct('mu',0,'km',0,'dmass',0,'atv',0,'btv',1);
    fields = fieldnames(inactive);
    for k = 1:numel(fields)
        if ~isfield(opts,fields{k})
            opts.(fields{k}) = inactive.(fields{k});
        end
    end
    T8 = f_option(opts,'t8',298.15);
    P8 = f_option(opts,'p8',101325);
    validateattributes(T8, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(P8, {'numeric'}, {'scalar','real','finite','positive'});
    state = imr.fluid_properties(p,P8,T8);
    derived = struct('rho8',state.liquid.rho,'c8',state.liquid.c, ...
        'cp',p.liquid.Cp,'rv',p.vapor.Rs,'lheat',state.latent_heat);
    fields = fieldnames(derived);
    for k = 1:numel(fields)
        if ~any(strcmp(names,fields{k}))
            opts.(fields{k}) = derived.(fields{k});
        end
    end
    % The existing wall closure analytically inverts an exponential law.
    % A table may set ambient vapor pressure but cannot enable that closure.
    if p.capabilities.mass_transfer && ~strcmp(p.saturation.model,'exp_inverse_temperature')
        error('IMR:fluidProfile', 'Mass transfer needs an invertible exponential saturation law.');
    end
    saturation = struct('pressure',0,'c2',0);
    if selected(3) || p.capabilities.mass_transfer
        saturation = imr.fluid_properties(p,[],T8);
    end
    if ~isfield(opts,'pv')
        opts.pv = saturation.pressure;
    end
    if selected(3) && opts.pv >= P8
        error('IMR:fluidBoiling', 'Ambient vapor pressure must be below p8.');
    end
    % Validate dimensional scales here instead of letting NaN/complex
    % dimensionless groups propagate into solver initialization.
    positive = {'rho8','c8','cp','rv'};
    for k = 1:numel(positive)
        validateattributes(opts.(positive{k}), {'numeric'}, {'scalar','real','finite','positive'});
    end
    nonnegative = {'surft','mu','pv','lheat'};
    for k = 1:numel(nonnegative)
        if isfield(opts,nonnegative{k})
            validateattributes(opts.(nonnegative{k}), {'numeric'}, {'scalar','real','finite','nonnegative'});
        end
    end
    if selected(1)
        transport = {'atg','btg','atv','btv'};
        for k = 1:numel(transport)
            validateattributes(opts.(transport{k}), {'numeric'}, {'scalar','real','finite','nonnegative'});
        end
        gas_K = opts.atg*T8 + opts.btg;
        vapor_K = opts.atv*T8 + opts.btv;
        if opts.atg <= 0 || gas_K <= 0 || vapor_K <= 0
            error('IMR:fluidThermal', 'Thermal runs need positive conductivities, positive gas slope and kappa > 1.');
        end
        if isfield(opts,'kappa') && opts.kappa <= 1
            error('IMR:fluidThermal', 'Profile thermal preparation needs kappa > 1.');
        end
    end
    if selected(2)
        validateattributes(opts.km, {'numeric'}, {'scalar','real','finite','positive'});
    end
    fluid = struct();
    fluid.profile = p;
    fluid.properties = state;
    fluid.SI = opts;
    fluid.Pv = opts.pv;
    keep = ~repelem(strcmp(names,'fluid'),2);
    fluid.input_pairs = [{'fluid',p}, args(keep)];
    fluid.assumption = 'Fixed ambient liquid density and sound speed for pressure-form KM';
    fields = fieldnames(opts);
    args = reshape([fields, struct2cell(opts)].',1,[]);
end

% Reference noncondensable gas, separate from the liquid/vapor material.
function [a, b] = f_gas_conductivity()
    a = 5.28e-5; % W/(m K^2).
    b = 1.165e-2; % W/(m K).
end

function value = f_option(opts, name, fallback)
    value = fallback;
    if isfield(opts,name)
        value = opts.(name);
    end
end

function f_required(required, names)
    for k = 1:numel(required)
        if ~any(strcmp(names,required{k}))
            error('IMR:fluidMissingInput', 'This profile requires explicit SI input %s.', required{k});
        end
    end
end
