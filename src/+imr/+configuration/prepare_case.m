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
        p.TVector = [0 p.TFin];
    end
    p.wave_poly = [];
    p.wave_dpoly = [];
    if isnumeric(p.wave_type) && isscalar(p.wave_type) && any(p.wave_type == -3:-1)
        files = {'d_ml.mat','d_mn.mat','d_hn.mat'};
        data = load(fullfile(fileparts(mfilename('fullpath')),'..','..','data',files{p.wave_type+4}),'poly','dpoly');
        p.wave_poly = data.poly;
        p.wave_dpoly = data.dpoly;
    end
    p.Pv_star = p.vapor*p.Pv_sat;
    if ~specified.p0star
        exponent = 3;
        if p.bubtherm == 0
            exponent = 3*p.kappa;
        end
        p.P0_star = (1+1/(p.We*p.Req)-p.Pv_star)*p.Req^exponent;
    end
    p.Pb_star = p.P0_star+p.Pv_star;
    p.iota = p.Km_star/p.Lt;
    if p.DRe == 0
        p.iDRe = 0;
    elseif p.stress >= 3
        p.iDRe = p.LAM/p.DRe;
    else
        p.iDRe = 1/p.DRe;
    end
    theta = p.Rv_star/p.Ra_star*(p.Pb_star/p.Pv_star-1);
    p.kv0 = 1/(1+theta);
    p.mv0 = p.Pv_star/p.Rv_star;
    p.ma0 = p.P0_star/p.Ra_star;
    p.solver = solver;
    p.stress0set = specified.stress0;
    p.tflag = specified.tfin+specified.tvector;
    p.tvector = p.TVector;
    p = imr.configuration.validate_case(p);
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
        opts = optimset('display','off');
        fun = @(x) p.Pv_star*(1+(p.ma0/x)*(p.Ra_star/p.Rv_star))-1- ...
            (1/p.We)*(p.Pv_star/(p.Rv_star*x))^(1/3);
        total = p.ma0+p.mv0;
        vapor_mass = fzero(fun,total,opts);
        while isnan(vapor_mass)
            total = total/1.11;
            vapor_mass = fzero(fun,total,opts);
        end
        Req = (p.Rv_star*vapor_mass/p.Pv_star)^(1/3);
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
        Szero = p.Szero(:);
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
