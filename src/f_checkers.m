% file f_checkers.m
% brief contains function f_checkers

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function checks the compatibility of the model conditions and
% the ranges of the dimensionless numbers before a run. It is called by
% f_call_params with a struct p holding every setting and dimensionless
% number. Incompatible or out-of-range inputs raise an error; inputs that
% are valid but suspicious (e.g. values that look dimensional) raise a
% warning. Finally, the settings tied to the stress model (spectral, Nv,
% Ca) are made consistent with it, and the adjusted struct is returned.
function p = f_checkers(p)
    
    % model flags and options
    
    f_integer_in(p.radial, 1:6, 'radial', ['1 (RP), 2 (K-M pressure), ' ...
        '3 (K-M enthalpy, Tait), 4 (Gilmore, Tait), 5 (K-M, Mie-Gruneisen), ' ...
        '6 (Gilmore, Mie-Gruneisen)']);
    f_integer_in(p.bubtherm, 0:1, 'bubtherm', '0 or 1');
    f_integer_in(p.medtherm, 0:1, 'medtherm', '0 or 1');
    f_integer_in(p.stress, 0:5, 'stress', 'an integer in 0..5');
    f_integer_in(p.vapor, 0:1, 'vapor', '0 or 1');
    f_integer_in(p.masstrans, 0:1, 'masstrans', '0 or 1');
    f_integer_in(p.collapse, 0:1, 'collapse', '0 or 1');
    f_integer_in(p.nu_model, 0:7, 'nu_model', ['0 (Newtonian), 1 (Carreau), ' ...
        '2 (Carreau-Yasuda), 3 (Powell-Eyring), 4 (modified Powell-Eyring), ' ...
        '5 (Cross), 6 (simplified Cross), or 7 (modified Cross)']);
    f_integer_in(p.wave_type, -3:3, 'wave_type', 'an integer in -3..3');
    f_integer_in(p.method, [15 23], 'method', '15 (ode15s) or 23 (ode23tb)');
    f_integer_in(p.spectral, 0:1, 'spectral', '0 or 1');
    f_integer_in(p.jacobian, 0:2, 'jacobian', '0, 1, or 2');
    f_integer_in(p.scheme, 0:1, 'scheme', '0 (finite difference) or 1 (Chebyshev)');
    f_integer_in(p.progdisplay, 0:1, 'progdisplay', '0 or 1');
    if ~(isnumeric(p.eps3) && isscalar(p.eps3) && p.eps3 >= 0 && p.eps3 <= 0.5)
        error('INPUT ERROR: eps3 must be in [0, 0.5]');
    end
    
    % compatibility of the model conditions
    
    if p.medtherm == 1 && p.bubtherm == 0
        error('INPUT ERROR: medtherm = 1 requires bubtherm = 1');
    end
    if p.masstrans == 1 && p.vapor == 0
        error('INPUT ERROR: masstrans = 1 requires vapor = 1');
    end
    if p.masstrans == 1 && p.bubtherm == 0
        error(['INPUT ERROR: masstrans = 1 requires bubtherm = 1 ' ...
            '(isothermal vapor diffusion is not implemented)']);
    end
    if p.radial >= 5 && (p.bubtherm == 1 || p.medtherm == 1 || p.masstrans == 1)
        error(['INPUT ERROR: the Mie-Gruneisen models (radial = 5, 6) support ' ...
            'only the polytropic bubble, set bubtherm = medtherm = masstrans = 0']);
    end
    if p.radial >= 3 && p.medtherm == 1
        error(['INPUT ERROR: the enthalpy forms (radial = 3, 4) take ' ...
            'the liquid at the far-field temperature, so medium heat ' ...
            'transfer (medtherm = 1) is not available. Use radial = 1 or 2']);
    end
    if strcmp(p.solver,'spectral') && p.masstrans == 1
        error(['INPUT ERROR: masstrans = 1 is not available in the spectral ' ...
            'solver, use f_imr_nodal']);
    end
    if strcmp(p.solver,'spectral') && p.scheme ~= 0
        error(['INPUT ERROR: scheme selects the differentiation matrices of ' ...
            'f_imr_nodal. f_imr_spectral takes no scheme']);
    end
    if strcmp(p.solver,'nodal') && p.spectral == 1 && p.stress == 5
        error(['INPUT ERROR: the spectral stress representation (spectral = 1) ' ...
            'is available in f_imr_spectral only']);
    end
    if p.collapse == 1 && (p.bubtherm ~= 1 || p.medtherm ~= 1 || ...
        p.masstrans ~= 1 || p.vapor ~= 1)
    error('INPUT ERROR: collapse must have full model with vapor');
end
if p.stress >= 3 && p.De <= 0
    error('INPUT ERROR: De must be positive for stress = 3, 4, or 5');
end
if p.stress >= 3 && p.medtherm == 1 && p.spectral == 0
    warning('IMR:dissipation', ['Medium heating for stress>=3 currently ' ...
        'uses the Kelvin-Voigt closed form (inconsistent for ' ...
        'history-dependent stress), see docs.']);
end

% time span

if p.tflag > 1
    error('INPUT ERROR: Only tvector or tfin can be specified, not both');
end
t = p.tvector;
if ~isnumeric(t) || numel(t) < 2 || any(~isfinite(t)) || any(diff(t) <= 0)
    error('INPUT ERROR: tvector must be finite and strictly increasing');
end
if t(1) < 0
    error('INPUT ERROR: tvector must start at a nonnegative time');
end
f_dimensional_hint(t(end) < 1e-3, 'the final time', t(end));

% solver and grid options

f_nonnegative_integer(p.Nt, 'Nt');
f_nonnegative_integer(p.Mt, 'Mt');
f_nonnegative_integer(p.Nv, 'Nv');
if p.bubtherm == 1 && p.Nt < 3
    error('INPUT ERROR: Nt must be at least 3 when bubtherm = 1');
end
if p.medtherm == 1 && p.Mt < 3
    error('INPUT ERROR: Mt must be at least 3 when medtherm = 1');
end
if ~(isnumeric(p.divisions) && p.divisions >= 0 && p.divisions == round(p.divisions))
    error('INPUT ERROR: divisions must be a nonnegative integer');
end
f_positive({'Lv','Lt','rtol','atol'}, {p.Lv, p.Lt, p.rtol, p.atol});

% initial conditions

f_positive({'Req'}, {p.Req});
f_dimensional_hint(p.Req < 1e-3, 'Req (= Req/R0)', p.Req);
f_finite({'Rdotzero'}, {p.Rdotzero});
if ~(isfinite(p.P0_star) && p.P0_star >= 0)
    error(['INPUT ERROR: the initial gas pressure P0_star = %g is negative; ' ...
        'check We, Req and the vapor pressure (1 + 1/(We*Req) must exceed ' ...
        'Pv_star)'], p.P0_star);
end

% dimensionless numbers

% Reynolds, Weber, and Cauchy numbers may be Inf (no viscosity, surface
% tension, or elasticity, respectively)
f_positive({'Re','We','Ca'}, {p.Re8, p.We, p.Ca}, true);
f_nonnegative({'De','LAM','alphax','kappa - 1'}, ...
    {p.De, p.LAM, p.alphax, p.kappa - 1});
if p.stress >= 3 && p.LAM > 1
    warning('IMR:parameters', ['LAM = lambda2/lambda1 = %g exceeds 1; the ' ...
        'retardation time is normally shorter than the relaxation time'], p.LAM);
end
f_positive({'Cstar','GAMa','nstate','hugoniot_s'}, ...
    {p.Cstar, p.GAMa, p.nstate, p.hugoniot_s});
f_dimensional_hint(p.Cstar < 1, 'Cstar (= C8/Uc)', p.Cstar);

% non-Newtonian viscosity
if p.nu_model > 0 && p.stress == 0
    error(['INPUT ERROR: nu_model = %d needs a viscous stress model; ' ...
        'stress = 0 has no viscous stress, use stress = 1 to 5'], p.nu_model);
end
if p.nu_model > 0
    f_nonnegative({'DRe','v_lambda','v_nc'}, {p.DRe, p.v_lambda_star, p.v_nc});
    if p.DRe == 0
        warning('IMR:parameters', ['nu_model = %d with DRe = 0: the ' ...
            'non-Newtonian contribution is zero'], p.nu_model);
    end
    if (p.nu_model == 2 || p.nu_model == 7) && ~(p.v_a > 0)
        error('INPUT ERROR: nu_model = %d requires v_a > 0', p.nu_model);
    end
    if (p.nu_model == 5 || p.nu_model == 7) && ~(p.v_nc > 0)
        error('INPUT ERROR: nu_model = %d requires v_nc > 0', p.nu_model);
    end
    if p.nu_model == 4 && p.v_nc > 1
        error(['INPUT ERROR: nu_model = 4 (modified Powell-Eyring) requires ' ...
            'v_nc <= 1, the zero-shear viscosity is unbounded otherwise']);
    end
    if p.nu_model == 4 && ~(p.v_lambda_star > 0)
        error('INPUT ERROR: nu_model = 4 requires v_lambda > 0');
    end
end

% vapor pressure
if p.vapor == 1
    f_nonnegative({'Pv'}, {p.Pv_sat});
    if p.Pv_star >= 1
        warning('IMR:parameters', ['Pv = %g is not below the far-field ' ...
            'pressure, the liquid is boiling'], p.Pv_star);
    end
end

% thermal transport
if p.bubtherm == 1
    f_positive({'chi'}, {p.chi});
    f_nonnegative({'alpha_g','beta_g','alpha_v','beta_v'}, ...
        {p.alpha_g, p.beta_g, p.alpha_v, p.beta_v});
end
if p.medtherm == 1
    f_positive({'Foh','kmstar'}, {p.Foh, p.Km_star});
    f_nonnegative({'Br'}, {p.Br});
end

% mass transfer
if p.masstrans == 1
    f_positive({'Fom','lheat','rv','ra','thetav','Pv'}, ...
        {p.Fom, p.L_heat_star, p.Rv_star, p.Ra_star, p.thetav, p.Pv_sat});
end
if p.collapse == 1
    f_positive({'rv','ra'}, {p.Rv_star, p.Ra_star});
end

% pressure waveform
f_finite({'ee','dt'}, {p.ee, p.dt});
if p.wave_type == 1
    f_positive({'tw'}, {p.tw});
    f_dimensional_hint(p.tw < 1e-3, 'tw (= TW/t0)', p.tw);
elseif p.wave_type == 2
    f_positive({'om'}, {p.om});
    f_nonnegative({'mn'}, {p.mn});
elseif p.wave_type == 3
    f_nonnegative({'tw'}, {p.tw});
elseif p.wave_type < 0
    if ~(isstruct(p.wave_poly) && isstruct(p.wave_dpoly))
        error('INPUT ERROR: waveform data for wave_type = %d not loaded', ...
            p.wave_type);
    end
    tb = p.wave_poly.breaks;
    if t(end) > tb(end) || t(1) < tb(1)
        warning('IMR:parameters', ['tvector [%g, %g] extends outside the ' ...
            'experimental waveform data [%g, %g]; the waveform is ' ...
            'extrapolated'], t(1), t(end), tb(1), tb(end));
    end
end

% settings consistent with the stress model

if p.stress <= 4
    % closed-form or single-mode stress models, no spectral stress
    % representation
    p.spectral = 0;
    if p.stress == 3 || p.stress == 4
        p.Nv = 1;
    else
        p.Nv = 0;
    end
else
    % UCM and Oldroyd-B: Nv stress modes for the spectral representation,
    % one otherwise, and no elastic contribution
    p.Nv = (p.Nv - 1)*(p.spectral == 1) + 1;
    p.Ca = Inf;
end

end

% helper checks

function f_integer_in(val, allowed, name, desc)
    if ~(isnumeric(val) && isscalar(val) && any(val == allowed))
        error('INPUT ERROR: %s must be %s', name, desc);
    end
end

function f_nonnegative_integer(val, name)
    if ~(isnumeric(val) && isscalar(val) && val >= 0 && val == round(val))
        error('INPUT ERROR: %s must be a nonnegative integer', name);
    end
end

function f_positive(names, vals, allowinf)
    if nargin < 3
        allowinf = false;
    end
    for k = 1:numel(names)
        v = vals{k};
        ok = isnumeric(v) && isscalar(v) && isreal(v) && v > 0 && ...
            (allowinf || isfinite(v));
        if ~ok
            error('INPUT ERROR: %s must be positive%s (got %g)', names{k}, ...
                repmat(' and finite', 1, ~allowinf), v);
        end
    end
end

function f_nonnegative(names, vals)
    for k = 1:numel(names)
        v = vals{k};
        if ~(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && v >= 0)
            error('INPUT ERROR: %s must be nonnegative and finite (got %g)', ...
                names{k}, v);
        end
    end
end

function f_finite(names, vals)
    for k = 1:numel(names)
        v = vals{k};
        if ~(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v))
            error('INPUT ERROR: %s must be finite (got %g)', names{k}, v);
        end
    end
end

% dimensionless inputs are O(1)-scaled; values this small usually mean a
% dimensional (SI) value was passed
function f_dimensional_hint(cond, name, val)
    if cond
        warning('IMR:dimensional', ['%s = %g looks dimensional; IMR takes ' ...
            'dimensionless inputs (see f_dim2nondim.m)'], name, val);
    end
end
