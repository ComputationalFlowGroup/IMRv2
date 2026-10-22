% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Prepare one dimensionless case; no spatial operators or integration.
% c groups checked model, numerics, initial state, forcing and physics.
% c.inputs records requested public inputs; c.layout describes effective
% states after compatibility adjustments (for example stress mode counts).
function c = prepare_case(solver,varargin)
    if ~any(strcmp(solver,{'nodal','spectral'}))
        error('IMR:solver','Solver must be nodal or spectral.');
    end
    disp('--- Inertial Microcavitation Rheometry forward solver ---');
    if ~isempty(varargin) && isstruct(varargin{1})
        if ~isscalar(varargin{1})
            error('INPUT ERROR: an input struct must be scalar.');
        end
        varargin = [f_struct_pairs(varargin{1}),varargin(2:end)];
    end
    if mod(numel(varargin),2) ~= 0
        error('INPUT ERROR: unmatched inputs');
    end
    p = imr.configuration.defaults();
    schema = imr.configuration.schema();
    names = {schema.name};
    requested = imr.default_inputs();
    specified = struct('p0star',false,'stress0',false,'tfin',false,'tvector',false);
    unknown = {};
    for k = 1:2:numel(varargin)
        key = varargin{k};
        valid_name = (ischar(key) && isrow(key)) || (isstring(key) && isscalar(key));
        if ~valid_name
            error('INPUT ERROR: parameter names must be text scalars.');
        end
        name = lower(char(key));
        value = varargin{k+1};
        index = find(strcmp(name,names),1);
        if ~isempty(index)
            p.(schema(index).field) = value;
        elseif strcmp(name,'p0star')
            p.P0_star = value;
        elseif strcmp(name,'stress0')
            p.Szero = value;
        elseif strcmp(name,'casefile')
            error('INPUT ERROR: casefile is no longer supported. Pass dimensionless name/value pairs.');
        else
            unknown{end+1} = name; %#ok<AGROW>
        end
        requested.(name) = value;
        if isfield(specified,name)
            specified.(name) = true;
        end
    end
    if ~isempty(unknown)
        error('INPUT ERROR: unrecognized input(s): %s. Use imr.nondimensionalize for SI inputs.',strjoin(unknown,', '));
    end
    p.method = f_method_code(p.method);
    if specified.tfin && ~specified.tvector
        valid_end = isnumeric(p.TFin) && isscalar(p.TFin) && isreal(p.TFin) && isfinite(p.TFin);
        if ~valid_end || p.TFin <= 0
            error('INPUT ERROR: tfin must be a positive, finite, real scalar.');
        end
        p.TVector = [0 double(p.TFin)];
    end
    p.wave_poly = [];
    p.wave_dpoly = [];
    if isnumeric(p.wave_type) && isscalar(p.wave_type) && any(p.wave_type == -3:-1)
        files = {'d_ml.mat','d_mn.mat','d_hn.mat'};
        data = load(fullfile(fileparts(mfilename('fullpath')),'..','..','data',files{p.wave_type+4}),'poly','dpoly');
        p.wave_poly = data.poly;
        p.wave_dpoly = data.dpoly;
    end
    p.solver = solver;
    p.stress0set = specified.stress0;
    p.pressure0set = specified.p0star;
    p.tflag = specified.tfin+specified.tvector;
    p.tvector = p.TVector;
    % Check primitive inputs before deriving pressures, mixture fractions
    % or transport ratios. The checker returns checked initial pressures.
    p = imr.configuration.validate_case(p);
    p.iota = p.Km_star/p.Lt;
    if p.DRe == 0
        p.iDRe = 0;
    elseif p.stress >= 3
        p.iDRe = p.LAM/p.DRe;
    else
        p.iDRe = 1/p.DRe;
    end
    p.kv0 = 0;
    p.mv0 = 0;
    p.ma0 = 0;
    if p.bubtherm && p.vapor && p.Pv_star > 0
        p.mv0 = p.Pv_star/p.Rv_star;
        p.ma0 = p.P0_star/p.Ra_star;
        p.kv0 = p.mv0/(p.mv0+p.ma0);
        alpha_m = p.kv0*p.alpha_v+(1-p.kv0)*p.alpha_g;
        beta_m = p.kv0*p.beta_v+(1-p.kv0)*p.beta_g;
        if ~(isfinite(p.kv0) && isfinite(alpha_m+beta_m) && alpha_m+beta_m > 0)
            error('INPUT ERROR: the initial mixture needs a finite vapor fraction and positive conductivity.');
        end
        if p.masstrans && p.kv0 >= 1
            error('INPUT ERROR: mass transfer requires a representable nonzero gas mass fraction.');
        end
    end
    if p.stress == 1 || p.stress == 2
        p.JdotA = 4/p.Re8;
    elseif p.stress >= 3
        p.JdotA = 4*p.LAM/p.Re8;
    else
        p.JdotA = 0;
    end
    Rdot = p.Rdotzero;
    Req = p.Req;
    if p.collapse
        % At equilibrium Pg0/Req^3 + Pv = 1 + 1/(We*Req).
        % Solve in positive radius, not vapor mass: scalar-start fzero
        % could probe negative mass, and its NaN retry loop was unbounded.
        pressure_gap = 1-p.Pv_star;
        balance = @(r) pressure_gap*r^3+r^2/p.We-p.P0_star;
        if pressure_gap > 0
            gas_radius = exp((log(p.P0_star)-log(pressure_gap))/3);
            if isinf(p.We)
                Req = gas_radius;
            else
                cap_radius = exp((log(p.P0_star)+log(p.We))/2);
                upper = min(gas_radius,cap_radius);
                % Enlarge by roundoff so the upper endpoint brackets.
                upper = upper*(1+16*eps);
                scaled_balance = @(x) balance(upper*x)/p.P0_star;
                Req = upper*fzero(scaled_balance,[0 1],optimset('display','off'));
            end
        elseif pressure_gap == 0 && isfinite(p.We)
            Req = exp((log(p.P0_star)+log(p.We))/2);
        elseif pressure_gap < 0 && isfinite(p.We)
            % The smaller root is stable; the cubic maximum brackets it.
            upper = 2/(3*(-pressure_gap)*p.We);
            if balance(upper) <= 0
                error('INPUT ERROR: collapse has no stable positive equilibrium radius.');
            end
            scaled_balance = @(x) balance(upper*x)/p.P0_star;
            Req = upper*fzero(scaled_balance,[0 1],optimset('display','off'));
        else
            error('INPUT ERROR: collapse has no stable positive equilibrium radius.');
        end
        if ~(isfinite(Req) && Req > 0)
            error('INPUT ERROR: collapse equilibrium radius must be positive and finite.');
        end
        Rdot = 0;
        fprintf('New equilibrium radius, Req = %.8f\n',Req);
    end
    nstress = 0;
    if p.stress == 3 || p.stress == 4
        nstress = 1;
    elseif p.stress == 5
        nstress = 2*p.Nv;
    end
    if specified.stress0
        valid_stress = isnumeric(p.Szero) && isreal(p.Szero) && all(isfinite(p.Szero(:)));
        if ~valid_stress || (~isempty(p.Szero) && ~isvector(p.Szero))
            error('INPUT ERROR: stress0 must be a finite, real numeric vector.');
        end
        Szero = double(p.Szero(:));
        if numel(Szero) ~= nstress
            error('INPUT ERROR: stress0 must have %d element(s) for stress = %d',nstress,p.stress);
        end
    elseif (p.stress == 3 || p.stress == 4) && p.collapse
        Szero = imr.physics.initial_stress(Req,p.Re8,p.Ca,p.De,p.We,p.Cstar,p.Pv_star,p.alphax*(p.stress == 4));
    else
        Szero = zeros(nstress,1);
    end
    c.solver = solver;
    c.model.radial = p.radial;
    c.model.bubtherm = p.bubtherm;
    c.model.medtherm = p.medtherm;
    c.model.stress = p.stress;
    c.model.eps3 = p.eps3;
    c.model.masstrans = p.masstrans;
    c.model.vapor = p.vapor;
    c.model.collapse = p.collapse;
    c.numerics.method = p.method;
    c.numerics.spectral = p.spectral;
    c.numerics.divisions = p.divisions;
    c.numerics.Nv = p.Nv;
    c.numerics.Nt = p.Nt;
    c.numerics.Mt = p.Mt;
    c.numerics.Lv = p.Lv;
    c.numerics.Lt = p.Lt;
    c.numerics.rtol = p.rtol;
    c.numerics.atol = p.atol;
    c.numerics.jacobian = p.jacobian;
    c.numerics.scheme = p.scheme;
    c.acoustic.Cstar = p.Cstar;
    c.acoustic.GAMa = p.GAMa;
    c.acoustic.kappa = p.kappa;
    c.acoustic.nstate = p.nstate;
    c.acoustic.hugoniot_s = p.hugoniot_s;
    c.forcing.om = p.om;
    c.forcing.ee = p.ee;
    c.forcing.tw = p.tw;
    c.forcing.dt = p.dt;
    c.forcing.mn = p.mn;
    c.forcing.wave_type = p.wave_type;
    c.forcing.wave_poly = p.wave_poly;
    c.forcing.wave_dpoly = p.wave_dpoly;
    c.stress.We = p.We;
    c.stress.Re8 = p.Re8;
    c.stress.DRe = p.DRe;
    c.stress.v_a = p.v_a;
    c.stress.v_nc = p.v_nc;
    c.stress.Ca = p.Ca;
    c.stress.LAM = p.LAM;
    c.stress.De = p.De;
    c.stress.JdotA = p.JdotA;
    c.stress.nu_model = p.nu_model;
    c.stress.v_lambda_star = p.v_lambda_star;
    c.stress.iDRe = p.iDRe;
    c.stress.alphax = p.alphax;
    c.thermal.Foh = p.Foh;
    c.thermal.Br = p.Br;
    c.thermal.alpha_g = p.alpha_g;
    c.thermal.beta_g = p.beta_g;
    c.thermal.alpha_v = p.alpha_v;
    c.thermal.beta_v = p.beta_v;
    c.thermal.chi = p.chi;
    c.thermal.iota = p.iota;
    c.mass.Fom = p.Fom;
    c.mass.kv0 = p.kv0;
    c.mass.Rv_star = p.Rv_star;
    c.mass.Ra_star = p.Ra_star;
    c.mass.L_heat_star = p.L_heat_star;
    c.mass.mv0 = p.mv0;
    c.mass.ma0 = p.ma0;
    c.mass.thetav = p.thetav;
    c.initial = struct('R',1,'Rdot',Rdot,'P',p.Pb_star, ...
        'Pg',p.P0_star,'Pv',p.Pv_star,'Req',Req,'stress',Szero);
    c.time = p.tvector;
    c.output.progress = p.progdisplay;
    requested.tvector = c.time;
    requested.tfin = c.time(end);
    c.inputs = requested;
    c.layout = imr.solvers.state_layout(c);
end

% integrator code of the method input, which may also be given by name:
% 'cvodes' (3), 'ode15s' (15) or 'ode23tb' (23)
function code = f_method_code(method)
    if ischar(method) || isstring(method)
        valid_text = (ischar(method) && isrow(method)) || (isstring(method) && isscalar(method));
        if ~valid_text
            error('INPUT ERROR: method must be a text scalar or numeric solver code.');
        end
        names = {'cvodes', 'ode15s', 'ode23tb'};
        codes = [3 15 23];
        k = find(strcmpi(char(method), names), 1);
        if isempty(k)
            error(['INPUT ERROR: method must be ''cvodes'', ''ode15s'' or ' ...
                '''ode23tb'' (or 3, 15, 23)']);
        end
        code = codes(k);
    else
        code = method;
    end
end

% name/value pairs of the fields of p that differ from imr.default_inputs; a time
% vector that ends at tfin makes tfin redundant, so only tvector is passed
function pairs = f_struct_pairs(p)
    d = imr.default_inputs();
    f = fieldnames(p);
    changed = false(size(f));
    for k = 1:numel(f)
        changed(k) = ~isfield(d, f{k}) || ~isequaln(p.(f{k}), d.(f{k}));
    end
    if isfield(p, 'tvector') && isfield(p, 'tfin')
        it = find(strcmp(f, 'tvector'));
        jf = find(strcmp(f, 'tfin'));
        valid_time = isnumeric(p.tvector) && isvector(p.tvector) && ~isempty(p.tvector);
        if changed(it) && changed(jf) && valid_time && isequal(p.tvector(end),p.tfin)
            changed(jf) = false;
        end
    end
    f = f(changed);
    pairs = cell(1, 2*numel(f));
    pairs(1:2:end) = f;
    for k = 1:numel(f)
        pairs{2*k} = p.(f{k});
    end
end
