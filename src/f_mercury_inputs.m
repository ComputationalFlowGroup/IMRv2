% file f_mercury_inputs.m
% brief prepares SI mercury properties for pressure-form Keller-Miksis

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% [varin, ref] = f_mercury_inputs('mu',0,'surft',0,'kappa',1.4, ...)
% [t, R] = f_imr_nodal(varin{:}) or f_imr_spectral(varin{:}).
% All physical inputs are SI, as in f_dim2nondim. mu and surft must be
% supplied explicitly (zero selects the corresponding idealized limit).
% kappa describes the noncondensable gas, not liquid or Hg-vapor gamma.
% rho8, C8, Cp, Rv and ambient Pv are derived, never water defaults.
% The standard KM equation holds ambient density and sound speed fixed.
function [varin, ref] = f_mercury_inputs(varargin)
    if mod(nargin,2) ~= 0
        error('IMR:mercuryInputs', 'Every mercury input needs a value.');
    end
    names = cell(1,nargin/2);
    values = varargin(2:2:end);
    allowed = {'t8','p8','r0','req','u0','tfin','tvector','surft','g', ...
        'mu','mu0','du','lambda1','lambda2','v_lambda','atg','btg', ...
        'atv','btv','km','dm','ra','p0','pa','omega','tw','dt', ...
        'kappa','alphax','mn','nu_model','v_a','v_nc','radial','stress', ...
        'bubtherm','medtherm','masstrans','collapse','vapor','method', ...
        'spectral','scheme','divisions','nt','mt','nv','lv','lt', ...
        'reltol','abstol','jacobian','progdisplay','stress0'};
    for k = 1:numel(names)
        name = varargin{2*k-1};
        if ~(ischar(name) && isrow(name)) && ~(isstring(name) && isscalar(name))
            error('IMR:mercuryInputs', 'Mercury input names must be text scalars.');
        end
        names{k} = lower(char(name));
        if ~any(strcmp(names{k}, allowed))
            error('IMR:mercuryOverride', ...
                ['Unsupported mercury input %s. Use SI inputs, not derived ' ...
                'density, acoustics, caloric or dimensionless overrides.'], names{k});
        end
    end
    f_required('surft', names);
    f_required('kappa', names);
    stress = f_value('stress',1,names,values);
    validateattributes(stress, {'numeric'}, {'scalar','integer','>=',0,'<=',5});
    if stress > 0
        f_required('mu', names);
    end
    radial = f_value('radial',2,names,values);
    if ~isequal(radial,2)
        error('IMR:mercuryRadial', ...
            'The mercury helper supports radial = 2 (KM pressure) only.');
    end
    flags = {'bubtherm','medtherm','vapor','masstrans','collapse'};
    selected = zeros(size(flags));
    for k = 1:numel(flags)
        selected(k) = f_value(flags{k},0,names,values);
        validateattributes(selected(k), {'numeric'}, {'integer','>=',0,'<=',1});
    end
    if selected(4) || selected(5)
        error('IMR:mercuryMassTransfer', ...
            ['Mercury dynamic mass transfer and collapse initialization ' ...
            'are not supported. Use masstrans = collapse = 0.']);
    end
    if selected(2) && ~selected(1)
        error('IMR:mercuryThermal', 'medtherm = 1 requires bubtherm = 1.');
    end
    if selected(1)
        transport = {'atg','btg','atv','btv'};
        for k = 1:numel(transport)
            f_required(transport{k}, names);
            value = f_value(transport{k},[],names,values);
            validateattributes(value, {'numeric'}, {'scalar','real','finite','nonnegative'});
        end
    end
    if selected(2)
        f_required('km', names);
        validateattributes(f_value('km',[],names,values), ...
            {'numeric'}, {'scalar','real','finite','positive'});
    end
    T8 = f_value('t8',298.15,names,values);
    P8 = f_value('p8',101325,names,values);
    validateattributes(T8, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(P8, {'numeric'}, {'scalar','real','finite','positive'});
    mu = f_value('mu',0,names,values);
    surft = f_value('surft',[],names,values);
    kappa = f_value('kappa',[],names,values);
    validateattributes(mu, {'numeric'}, {'scalar','real','finite','nonnegative'});
    validateattributes(surft, {'numeric'}, {'scalar','real','finite','nonnegative'});
    validateattributes(kappa, {'numeric'}, {'scalar','real','finite','>=',1});
    hg = f_mercury_eos(P8,T8);
    if selected(1)
        gas_K = f_value('atg',[],names,values)*T8 + f_value('btg',[],names,values);
        vapor_K = f_value('atv',[],names,values)*T8 + f_value('btv',[],names,values);
        gas_slope = f_value('atg',[],names,values);
        if gas_K <= 0 || vapor_K <= 0 || gas_slope <= 0 || kappa <= 1
            error('IMR:mercuryThermal', ...
                ['Thermal runs need positive conductivities, positive gas ' ...
                'conductivity slope and kappa > 1.']);
        end
    end
    p = f_mercury_parameters();
    Pv = 0;
    if selected(3)
        Pv = f_mercury_pvsat(T8);
        if Pv >= P8
            error('IMR:mercuryBoiling', ...
                'Ambient mercury vapor pressure must be below p8.');
        end
    end
    % Defaults here select a Newtonian liquid, not the converter's elastic
    % water case. Explicit user pairs still override these generic choices.
    args = {'radial',2,'stress',1,'g',0,'lambda1',0,'lambda2',0,'mu',0};
    args = [args, varargin, {'rho8',hg.liquid.rho,'c8',hg.liquid.c, ...
        'cp',p.liquid.Cp,'rv',p.vapor.Rs,'pv',Pv, ...
        'lheat',hg.latent_heat,'thetav',0}];
    [varin, ref] = f_dim2nondim(args{:});
    % No wall saturation law is enabled. thetav = 0 must not be used with
    % mass transfer by modifying the returned pairs to bypass this guard.
    ref.mercury.source = p.source;
    ref.mercury.commit = p.commit;
    ref.mercury.coefficients = p;
    ref.mercury.properties = hg;
    ref.mercury.Pv = Pv;
    ref.mercury.input_pairs = varargin;
    ref.mercury.assumption = 'Fixed ambient liquid density and sound speed';
end

function value = f_value(name, fallback, names, values)
    k = find(strcmp(names,name),1,'last');
    value = fallback;
    if ~isempty(k)
        value = values{k};
    end
end

function f_required(name, names)
    if ~any(strcmp(names,name))
        error('IMR:mercuryMissingInput', ...
            'Mercury requires explicit SI input %s. No water value is assumed.', name);
    end
end
