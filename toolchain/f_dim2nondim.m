% file f_dim2nondim.m
% brief contains function f_dim2nondim

% Copyright (C) Mauro Rodriguez
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This side function converts dimensional (SI) physical properties
% into the dimensionless inputs taken by the IMR solvers (f_imr_nodal,
% f_imr_spectral). It is not called by the solvers: the solvers accept only
% dimensionless inputs, and this helper is provided for users who start
% from dimensional values.
%
% Usage:
%   [varin, ref] = f_dim2nondim('r0',50e-6,'req',5e-6,'mu',1e-2,'g',1e3);
%   [t, R] = f_imr_nodal(varin{:}, 'radial', 2, 'Nt', 50);
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
%   pv (Pa) f_pvsat(t8), p0 (Pa) computed by the solver if not given
%   pa (Pa) 0, omega (rad/s) 0, tw (s) 0, dt (s) 0
% Any other name/value pair (model flags, grid sizes, kappa, alphax, mn,
% nstate, hugoniot_s, v_a, v_nc, ...) is passed through unchanged.
%
% Outputs:
%   varin   cell array of dimensionless name/value pairs for the solvers
%   ref     reference scales: R0 (m), Uc (m/s), t0 (s), P8 (Pa), T8 (K),
%           rho8 (kg/m^3), to restore dimensions of the solver outputs
function [varin, ref] = f_dim2nondim(varargin)
    
    if mod(nargin,2) == 1
        error('f_dim2nondim: unmatched inputs');
    end
    
    % dimensional defaults (SI)
    R0 = 50e-6;
    Req = [];
    U0 = 0;
    TFin = [];
    TVector = [];
    rho8 = 1064;
    P8 = 101325;
    C8 = 1484;
    GAM = 3049.13*1e5;
    S = 0.07;
    G = 1e3;
    mu8 = 8.3283e-4;
    muo = [];
    Dmu = [];
    lambda1 = 1e-7;
    lambda2 = 1e-8;
    v_lambda = 0;
    T8 = 298.15;
    ATg = 5.28e-5;
    BTg = 1.165e-2;
    ATv = 3.30e-5;
    BTv = 1.742e-2;
    Km = 0.55;
    Cp = 4.181e3;
    Dm = [];
    D0 = 24.2e-6;
    L_heat = 2264.76e3;
    % (J/mol-K) universal gas constant, (kg/mol) molecular weights
    Ru = 8.3144598;
    Rv = Ru/18.01528E-3;
    Ra = Ru/28.966E-3;
    Pv = [];
    P0 = [];
    pA = 0;
    omega = 0;
    TW = 0;
    DT = 0;
    
    % read the dimensional inputs, pass everything else through
    passthrough = {};
    for n = 1:2:nargin
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
    [Pvsat, ~, pv_c2] = f_pvsat(T8);
    if isempty(Pv)
        Pv = Pvsat;
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
        'thetav',   -pv_c2/T8, ...
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
    
end
