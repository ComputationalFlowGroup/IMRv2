% file f_imr_nodal.m
% brief contains function f_imr_nodal

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function features the nodal solver of the PDEs involving
% thermal transport, vapor mass transfer and viscoelasticity coupled to the
% radial (Rayleigh-Plesset family) dynamics. The fields are carried as node
% values, and the spatial discretization is selected by the differentiation
% matrices of f_diff_mat: scheme = 0 second-order finite differences,
% scheme = 1 Chebyshev collocation. Everything else is shared.
function varargout = f_imr_nodal(varargin)
    
    % problem initialization
    [eqns_opts, solve_opts, init_opts, init_stress, tspan_opts, out_opts, ...
        acos_opts, wave_opts, sigma_opts, thermal_opts, mass_opts] ...
        = f_call_params('nodal', varargin{:});
    
    % defaults
    
    % non-Newtonian viscosity parameters
    fnu_m           = 0;
    intfnu          = 0;
    dintfnu         = 0;
    ddintfnu        = 0;
    
    % equations settings
    radial          = eqns_opts(1);
    bubtherm        = eqns_opts(2);
    medtherm        = eqns_opts(3);
    stress          = eqns_opts(4);
    eps3            = eqns_opts(5);
    masstrans       = eqns_opts(6);
    
    % solver options
    method          = solve_opts(1);
    spectral        = solve_opts(2);
    divisions       = solve_opts(3);
    Nv              = solve_opts(4);
    Nt              = solve_opts(5);
    Mt              = solve_opts(6);
    % Lv            = solve_opts(7), only used by the spectral solver
    Lt              = solve_opts(8);
    rtol            = solve_opts(9);
    atol            = solve_opts(10);
    jacobian        = solve_opts(11);
    scheme          = solve_opts(12);
    
    % dimensionless initial conditions
    Rzero           = init_opts(1);
    Rdotzero        = init_opts(2);
    Pb_star         = init_opts(3);
    Pv_star         = init_opts(4);
    Req             = init_opts(5);
    alphax          = init_opts(6);
    
    % dimensionless initial stress
    Szero           = init_stress;
    
    % time span options
    tspan           = tspan_opts;
    
    % output options
    progdisplay     = out_opts(1);
    
    % physical parameters
    
    % acoustic parameters
    Cstar           = acos_opts(1);
    GAMa            = acos_opts(2);
    kappa           = acos_opts(3);
    nstate          = acos_opts(4);
    hugoniot_s      = acos_opts(5);
    
    % dimensionless waveform parameters
    om              = wave_opts{1}(1);
    ee              = wave_opts{1}(2);
    tw              = wave_opts{1}(3);
    dt              = wave_opts{1}(4);
    mn              = wave_opts{1}(5);
    wave_type       = wave_opts{1}(6);
    wave_poly       = wave_opts{2};
    wave_dpoly      = wave_opts{3};
    
    pvarargin = {[om,ee,tw,dt,mn,wave_type], wave_poly, wave_dpoly};
    % the impulse waveform is the constant pressure ee (f_pinfinity returns
    % exactly ee and 0), so the right-hand side skips the call
    wave_const = (wave_type == 0);
    
    % dimensionless viscoelastic
    We              = sigma_opts(1);
    Re8             = sigma_opts(2);
    DRe             = sigma_opts(3);
    v_a             = sigma_opts(4);
    v_nc            = sigma_opts(5);
    Ca              = sigma_opts(6);
    LAM             = sigma_opts(7);
    De              = sigma_opts(8);
    JdotA           = sigma_opts(9);
    nu_model        = sigma_opts(10);
    v_lambda_star   = sigma_opts(11);
    iDRe            = sigma_opts(12);
    iWe             = 1/We;
    
    % dimensionless thermal
    Foh             = thermal_opts(1);
    Br              = thermal_opts(2);
    alpha_g         = thermal_opts(3);
    beta_g          = thermal_opts(4);
    alpha_v         = thermal_opts(5);
    beta_v          = thermal_opts(6);
    chi             = thermal_opts(7);
    iota            = thermal_opts(8);
    
    % dimensionless mass transfer
    Fom             = mass_opts(1);
    kv0             = mass_opts(2);
    Rv_star         = mass_opts(3);
    Rg_star         = mass_opts(4);
    L_heat_star     = mass_opts(5);
    thetav          = mass_opts(8);
    
    % pre_process
    
    % node positions and differentiation matrices: y = r/R in [0, 1] from
    % the bubble center to the wall, xi in [1, -1] from the bubble wall to
    % infinity in the medium
    [y, dmatrix_Tkv] = f_diff_mat(Nt,1,0,scheme);
    [~, ddmatrix_Tkv] = f_diff_mat(Nt,2,0,scheme);
    [xi, dmatrix_Tm] = f_diff_mat(Mt,1,1,scheme);
    [~, ddmatrix_Tm] = f_diff_mat(Mt,2,1,scheme);
    
    % outside the bubble
    yT   = ((2./(xi+1)-1)*Lt+1);
    yT2  = yT.^2;
    yT3  = yT.^3;
    iyT3 = yT.^-3;
    iyT4 = yT.^-4;
    iyT6 = yT.^-6;
    % constant factors of the medium temperature equation, formed with the
    % same operations as in the right-hand side, so results are unchanged
    xi1_sq = (1+xi).^2;
    xi1_p4 = (xi+1).^4;
    om_yT3 = 1-yT3;
    xi_adv = (xi+1)/(2*Lt)-1./yT;
    
    % precomputations move these to f_call_params
    %LDR = LAM*De/Re8;
    sam = 1 + GAMa;
    no = (nstate-1)/nstate;
    nog = (nstate-1)/2;
    kapover = (kappa-1)/kappa;
    % Kirchhoff map T(theta) of the bubble without mass transfer, whose
    % mixture coefficient is the constant kv0 (as in f_theta_of_T)
    alpha_m0 = kv0.*alpha_v + (1-kv0).*alpha_g;
    alpha_m0m1 = alpha_m0 - 1;
    om_kv0 = 1-kv0;
    % dimensionless saturated vapor pressure, Pv(T) = C1_pv*exp(C2_pv/T)
    % with Pv(1) = Pv_star
    C1_pv = Pv_star*exp(thetav);
    C2_pv = -thetav;
    Rvg_ratio = Rv_star / Rg_star;
    inv_alpha = 1/alpha_g;
    Rva_diff = Rv_star - Rg_star;
    % bubble-wall flux balance coefficients: the wall gradient rows of the
    % medium (first entry at the wall) and of the bubble (first entry at
    % the wall, then inward), with the node indices of the other entries
    if scheme == 0
        % second-order one-sided stencils
        coeff = [-1.5 , 2 ,-0.5 ];
        deltaY = 1/(Nt-1);
        deltaYm = -2/(Mt-1);
        grad_Tm_coeff = 2*chi*iota/deltaYm*coeff;
        grad_Trans_coeff = -coeff*chi/deltaY;
        grad_C_coeff = -coeff*Fom*L_heat_star/deltaY;
        iwall_m = [2 3];
        iwall_b = [Nt-1 Nt-2];
    else
        % Chebyshev: the wall rows couple every node
        grad_Tm_coeff = 2*chi*iota*full(dmatrix_Tm(1,:));
        grad_Trans_coeff = chi*full(dmatrix_Tkv(end,end:-1:1));
        grad_C_coeff = Fom*L_heat_star*full(dmatrix_Tkv(end,end:-1:1));
        iwall_m = 2:Mt;
        iwall_b = Nt-1:-1:1;
    end
    
    % index management
    if masstrans
        Nc = Nt;
    else
        Nc = 0;
    end
    if bubtherm == 0
        Nt = 0;
        Mt = 0;
    end
    if medtherm == 0
        Mt = 0;
    end
    ibubtherm   = 4:(3+Nt);
    imedtherm   = (4+Nt):(3+Nt+Mt);
    imass       = (4+Nt+Mt):(3+Nt+Mt+Nc);
    
    % precomputations for viscous dissipation
    % zT = 1 - 2./(1 + (yT - 1)/Lv);
    ZZT = 0;
    % stress spectra index management
    ivisco1 = (4+Nt+Mt+Nc):(3+Nt+Mt+Nc+Nv);
    ivisco2 = (4+Nt+Mt+Nc+Nv):(3+Nt+Mt+Nc+2*Nv);
    
    % initial condition assembly
    
    % bubble temperature initial condition
    if bubtherm
        theta0 = zeros(Nt,1);
    else
        theta0 = [];
    end
    
    % medium initial temperature
    if medtherm
        Tm0 = ones(Mt,1);
    else
        Tm0 = [];
    end
    
    % mass transfer initial condition
    if masstrans
        kv0vec = kv0*ones(Nc,1);
    else
        kv0vec = [];
    end
    
    % initial condition vector
    init = [Rzero;
    Rdotzero;
    Pb_star;
    theta0;
    Tm0;
    kv0vec;
    Szero];
    
    % thermal-only wall closure: previous wall value and fzero fallback
    theta_bw_guess = -0.0001;
    foptions = optimset('TolFun',1e-12);
    % mass-transfer wall closure: fixed bracket on the wall vapor mass
    % fraction and fixed iteration counts for the bracketed solve, which a
    % Newton warm start from the previous root replaces when it converges
    kv_bracket = [1e-13, 1 - 1e-13];
    n_grid = 64;
    n_rounds = 2;
    n_polish = 3;
    kv_w_guess = NaN;
    n_warm = 6;
    warm_tol = 1e-10;
    % constants of the bubble-wall heat (and mass) flux balance
    wallc = struct('grad_Tm',grad_Tm_coeff, 'grad_Tr',grad_Trans_coeff, ...
        'grad_C',grad_C_coeff, 'iwall_m',iwall_m, 'iwall_b',iwall_b, ...
        'alpha_g',alpha_g, 'alpha_v',alpha_v, ...
        'C1',C1_pv, 'C2',C2_pv, 'Rvg',Rvg_ratio, 'Rva',Rva_diff, ...
        'Rg',Rg_star);
    
    % solver start
    if scheme == 0
        scheme_name = 'nodal, finite difference';
    else
        scheme_name = 'nodal, Chebyshev';
    end
    f_display(radial, bubtherm, medtherm, masstrans, stress, spectral,...
        nu_model, eps3, Pv_star, Re8, De, Ca, LAM, scheme_name);
    bubble = @SVBDODE;
    % Jacobian for the stiff integrators: analytic where implemented,
    % otherwise the sparsity pattern
    if jacobian == 2 && bubtherm
        jp = struct('ith',ibubtherm, 'itm',imedtherm, 'ikv',imass, ...
            'iz',(4+Nt+Mt+Nc):numel(init), 'medtherm',medtherm, ...
            'masstrans',masstrans, 'D1',dmatrix_Tkv, 'D2',ddmatrix_Tkv, ...
            'Dm1',dmatrix_Tm, 'Dm2',ddmatrix_Tm, 'y',y, 'xi',xi, 'yT',yT, ...
            'yT2',yT2, 'yT3',yT3, 'Lt',Lt, 'Foh',Foh, 'chi',chi, ...
            'kappa',kappa, 'kapover',kapover, 'alpha_g',alpha_g, ...
            'beta_g',beta_g, 'alpha_v',alpha_v, 'beta_v',beta_v, 'kv0',kv0, ...
            'Fom',Fom, 'Rv',Rv_star, 'Rg',Rg_star, 'Rva',Rva_diff, ...
            'Rvg',Rvg_ratio, 'C1_pv',C1_pv, 'C2_pv',C2_pv, 'wall',wallc, ...
            'closure',@f_closure, 'taugradu',@f_taugradu, ...
            'theta_bw',@f_theta_bw_of_X);
        Sp = @(t,X) f_jacobian_nodal(t,X,jp);
    elseif jacobian >= 1
        Sp = f_jpattern(numel(init), ibubtherm, imedtherm, imass, ...
            (4+Nt+Mt+Nc):numel(init), scheme);
    else
        Sp = [];
    end
    [t,X] = f_odesolve(bubble, init, method, divisions, tspan, Sp, rtol, atol);
    
    % extract result
    R    = X(:,1);
    Rdot = X(:,2);
    P    = X(:,3);
    if bubtherm
        theta = X(:,ibubtherm);
        if medtherm
            Tm = X(:,imedtherm);
            % the wall entries of the state are not integrated (their time
            % derivatives are zero), so the wall state is recovered from the
            % wall closure at every output time, as the right-hand side does
            % (the loop index is named it_out because a nested function shares
            % any variable name that the parent also uses)
            kv_w = zeros(numel(t),1);
            theta_bw_guess = theta(1,end);
            kv_w_guess = NaN;
            for it_out = 1:numel(t)
                if masstrans
                    [theta(it_out,end), kv_w(it_out)] = f_wall_mass_closure( ...
                        Tm(it_out,:)', theta(it_out,:)', X(it_out,imass)', ...
                        P(it_out), R(it_out), t(it_out));
                else
                    theta(it_out,end) = f_wall_theta_exact(Tm(it_out,:)', ...
                        theta(it_out,:)');
                    theta_bw_guess = theta(it_out,end);
                end
            end
        end
    end
    if masstrans
        % T with the stored kv(end), as in the right-hand side, then the
        % wall vapor concentration
        kv = X(:,imass);
        T = f_theta_of_T(theta,kv);
        if medtherm
            kv(:,end) = kv_w;
        else
            kv(:,end) = f_kv_of_T(T(:,end),P);
        end
    elseif bubtherm
        T = f_theta_of_T(theta,kv0);
    end
    if medtherm
        Tm(:,1) = T(:,end);
    end
    
    % outputs assembly
    varargout{1} = t;
    varargout{2} = R;
    varargout{3} = Rdot;
    varargout{4} = P;
    if bubtherm
        varargout{5} = T;
    else
        varargout{5} = [];
    end
    if medtherm == 1
        varargout{6} = Tm;
    else
        varargout{6} = ones(size(t,1),1);
    end
    if masstrans == 1
        varargout{7} = kv;
    else
        varargout{7} = [];
    end
    
    % solver function
    function [dXdt] = SVBDODE(t,X)
        
        % showing output
        if progdisplay
            disp(t/tspan(end));
        end
        
        % extract solution
        
        % bubble wall radius
        R = X(1);
        % bubble wall velocity
        Rdot = X(2);
        % internal bubble pressure
        P = X(3);
        % auxiliary variable for internal bubble temperature
        if bubtherm || masstrans
            theta = X(ibubtherm);
        end
        if masstrans
            kv = X(imass);
        end
        
        % empty output variables
        thetadot = [];
        Tmdot = [];
        kvdot = [];
        
        % boundary condition evaluation
        if medtherm && masstrans
            % temperature in the material
            Tm = X(imedtherm);
            % bracketed closure on the wall vapor mass fraction
            [theta_bw, kv_bw] = f_wall_mass_closure(Tm, theta, kv, P, R, t);
            theta(end) = theta_bw;
        elseif medtherm
            % temperature in the material
            Tm = X(imedtherm);
            % closed-form wall temperature, fzero only as a fallback
            theta_bw = f_wall_theta_exact(Tm);
            theta_bw_guess = theta_bw;
            theta(end) = theta_bw;
        end
        
        % bubble temperature, f_theta_of_T written out
        if masstrans
            alpha_T = kv.*alpha_v + (1-kv).*alpha_g;
            T = (alpha_T - 1 + sqrt(1+2*theta.*alpha_T)) ./ alpha_T;
        elseif bubtherm
            T = (alpha_m0m1 + sqrt(1+2*theta.*alpha_m0)) ./ alpha_m0;
        end
        
        % bubble equations of motion
        if bubtherm && masstrans
            % vapor concentration at the boundary. With the medium, the
            % converged closure value; it is set after T so that T(end), like
            % the closure, uses the mixture coefficient of the stored kv(end)
            if medtherm
                kv(end) = kv_bw;
            else
                kv(end) = f_kv_of_T(T(end),P);
            end
            
            % temperature field gradients inside the bubble
            dtheta  = dmatrix_Tkv*theta;
            ddtheta = ddmatrix_Tkv*theta;
            
            % vapor concentration gradients inside the bubble
            dkv     = dmatrix_Tkv*kv;
            ddkv    = ddmatrix_Tkv*kv;
            om_kv   = 1-kv;
            Rmix    = kv*Rv_star + om_kv*Rg_star;
            wRmix   = Rva_diff./Rmix;
            RDkv    = wRmix.*dkv;
            
            % bubble pressure equation
            Pdot = 3/R*(chi*(kappa-1)*dtheta(end)/R - kappa*P*Rdot +...
                + kappa*P*Fom*Rv_star*dkv(end)/(T(end)*R*Rmix(end)*om_kv(end)));
            
            % bubble wall velocity equation
            Uvel = (chi/R*(kappa-1).*dtheta-y*R*Pdot/3)/(kappa*P) + Fom/R*RDkv;
            
            % mixture thermal conductivity
            alpha_mix  = kv.*alpha_v + om_kv.*alpha_g;
            Kstar_g = alpha_g*T+beta_g;
            Kstar_v = alpha_v*T+beta_v;
            Kstar   = kv.*Kstar_v + om_kv.*Kstar_g;
            
            % temperature evolution inside the bubble
            nonlinear_term = (chi*ddtheta./R^2+Pdot).*(kapover*Kstar.*T/P);
            % gas velocity relative to the moving grid
            u_rel = Uvel-y*Rdot;
            advection_term = -dtheta.*u_rel/R;
            mass_diffusion = (Fom/(R^2)).*wRmix.*dkv.*dtheta;
            thetadot = advection_term + nonlinear_term + mass_diffusion;
            % temperature at bubble wall, handled by flux boundary bc
            thetadot(end) = 0;
            
            % water vapor evolution inside the bubble
            nonlinear_diffusion = ...
                dkv.*(dtheta./(sqrt(1+2*alpha_mix.*theta).*T)+RDkv);
            advection_term = u_rel/R.*dkv;
            kvdot = Fom/R^2*(ddkv - nonlinear_diffusion) - advection_term;
            % concentration at bubble wall, handled by flux boundary bc
            kvdot(end) = 0;
            
        elseif bubtherm
            % mixture thermal conductivity
            Kstar_g = alpha_g*T+beta_g;
            Kstar_v = alpha_v*T+beta_v;
            Kstar   = kv0.*Kstar_v + om_kv0.*Kstar_g;
            
            % bubble temperature gradients
            dtheta  = dmatrix_Tkv*theta;
            ddtheta = ddmatrix_Tkv*theta;
            
            % bubble pressure evolution equation
            Pdot = 3/R*(chi*(kappa-1)*dtheta(end)/R - kappa*P*Rdot);
            
            % internal bubble velocity equation
            Uvel = (chi/R*(kappa-1).*dtheta-y*R*Pdot/3)/(kappa*P);
            
            % bubble temperature evolution equation
            diffusion_term = (chi*ddtheta./R^2+Pdot).*(kapover*Kstar.*T/P);
            advection_term = -dtheta.*(Uvel-y*Rdot)./R;
            thetadot = advection_term + diffusion_term;
            thetadot(end) = 0;
            
        else
            % polytropic gas + vapor
            P = (Pb_star-Pv_star)*(1/R)^(3*kappa) + Pv_star;
            Pdot = -3*kappa*P*Rdot/R;
        end
        
        % surroundings equations of motion
        
        % viscous forces/Reynolds number
        if nu_model ~= 0
            % fnu_m: viscosity factor at the medium grid points
            [~,intfnu,dintfnu,ddintfnu,fnu_m] = f_viscosity(nu_model, ...
                Rdot,R,v_a,v_nc,v_lambda_star,yT);
        end
        
        % surrounding temperature
        if medtherm
            % boundary temperature
            Tm(1) = T(end);
            % temperature gradients
            dTm = dmatrix_Tm*Tm;
            ddTm = ddmatrix_Tm*Tm;
            % material temperature equations
            advection = xi1_sq./(Lt*R).*(Rdot./yT2.*om_yT3/2 +...
                Foh/R.*xi_adv).*dTm;
            diffusion = Foh/R^2.*xi1_p4/Lt^2.*ddTm/4;
            % stress dissipation
            [taugradu] = f_stress_dissipation(stress,spectral,Req,R,Rdot, ...
                Ca,Br,Re8,alphax,yT2,yT3,iyT3,iyT4,iyT6,X,ZZT,ivisco1, ...
                ivisco2,fnu_m,DRe);
            % surrounding temperature evolution equation
            Tmdot = advection + diffusion + taugradu;
            % sets boundary condition on temperature
            Tmdot(1) = 0;
            Tmdot(end) = 0;
        end
        
        % pressure waveform
        if wave_const
            Pf8 = ee;
            Pf8dot = 0;
        else
            [Pf8,Pf8dot] = f_pinfinity(t,pvarargin);
        end
        
        % stress equation evolution
        [S,Sdot,Z1dot,Z2dot] = f_stress(stress,X,Req,R,Ca,De,Re8, ...
            Rdot,alphax,ivisco1,ivisco2,LAM,intfnu,dintfnu,iDRe);
        
        % bubble wall evolution / acceleration
        [Rddot] = f_radial_eq(radial, P, Pdot, Pf8, Pf8dot, iWe, R, Rdot, S, ...
            Sdot, Cstar, sam, no, GAMa, nstate, nog, hugoniot_s, JdotA, ...
            ddintfnu, iDRe);
        
        % output assembly
        dXdt = [Rdot;
        Rddot;
        Pdot;
        thetadot;
        Tmdot;
        kvdot;
        Z1dot;
        Z2dot];
        
    end
    % end of solver
    
    % medium dissipation as a function of (R, Rdot)
    function tg = f_taugradu(Rv,Rdotv,X)
        fnu_loc = 0;
        if nu_model ~= 0
            [~,~,~,~,fnu_loc] = f_viscosity(nu_model,Rdotv,Rv,v_a,v_nc, ...
                v_lambda_star,yT);
        end
        tg = f_stress_dissipation(stress,spectral,Req,Rv,Rdotv, ...
            Ca,Br,Re8,alphax,yT2,yT3,iyT3,iyT4,iyT6,X,ZZT,ivisco1, ...
            ivisco2,fnu_loc,DRe);
        tg = tg.*ones(size(yT));
    end
    
    % radial acceleration and stress evolution as a function of
    % v = [R, Rdot, P, Pdot, Z], mirroring SVBDODE
    function out = f_closure(t,v,X)
        Rv = v(1);
        Rdotv = v(2);
        X(4+Nt+Mt+Nc:end) = v(5:end);
        intf_loc = 0;
        dintf_loc = 0;
        ddintf_loc = 0;
        if nu_model ~= 0
            [~,intf_loc,dintf_loc,ddintf_loc] = f_viscosity(nu_model, ...
                Rdotv,Rv,v_a,v_nc,v_lambda_star);
        end
        if wave_const
            Pf8_loc = ee;
            Pf8dot_loc = 0;
        else
            [Pf8_loc,Pf8dot_loc] = f_pinfinity(t,pvarargin);
        end
        [S_loc,Sdot_loc,Z1dot_loc,Z2dot_loc] = f_stress(stress,X,Req,Rv, ...
            Ca,De,Re8,Rdotv,alphax,ivisco1,ivisco2,LAM, ...
            intf_loc,dintf_loc,iDRe);
        Rddot_loc = f_radial_eq(radial, v(3), v(4), Pf8_loc, Pf8dot_loc, ...
            iWe, Rv, Rdotv, S_loc, Sdot_loc, Cstar, sam, no, GAMa, nstate, ...
            nog, hugoniot_s, JdotA, ddintf_loc, iDRe);
        out = vertcat(Rddot_loc, Z1dot_loc, Z2dot_loc);
    end
    
    % solver functions
    
    % temperature at the bubble wall as a function of theta
    function Tw = f_theta_of_T(theta_w,kv)
        alpha_m  = kv.*alpha_v + (1-kv).*alpha_g;
        Tw = (alpha_m - 1 + sqrt(1+2*theta_w.*alpha_m)) ./ alpha_m;
    end
    
    % concentration at the bubble wall
    function kv_w = f_kv_of_T(Tw,pressure)
        theta_var = Rvg_ratio*(pressure./(C1_pv*exp(C2_pv./Tw))-1);
        kv_w = 1./(1+theta_var);
    end
    
    % bubble-wall state with vapor mass transfer. The unknown is the wall
    % vapor mass fraction kv_w, confined to (0, 1). For a trial kv_w the
    % saturation condition gives the wall temperature and the Kirchhoff
    % transform gives theta_w in closed form (f_wall_mass_residual).
    % Warm start: Newton with the analytic derivative from the previous root
    % kv_w_guess, accepted once a step falls below warm_tol with every
    % iterate inside kv_bracket, within n_warm steps. Otherwise (first call,
    % divergence, a non-finite value) the bracketed solve, which needs no
    % initial guess: n_rounds rounds of n_grid residual evaluations, one
    % vectorized call each, narrow kv_bracket to the grid interval with the
    % sign change (the first round includes the bracket ends), and n_polish
    % Newton steps follow.
    function [theta_w, kv_w] = f_wall_mass_closure(Tm_in, th_in, kv_in, Pw, Rw, tw)
        % mixture coefficients from the stored wall value kv_in(end), not the
        % trial kv_w. Since kvdot(end) = 0 the stored value is the initial
        % kv0. Deliberate carry-over of the existing behavior, to be
        % revisited separately.
        alpha_m = kv_in(end)*alpha_v + (1 - kv_in(end))*alpha_g;
        % the transform normalizes the conductivity to 1 at T = 1:
        % K = alpha_m*T + beta_m with beta_m = 1 - alpha_m, so that
        % theta = 0.5*alpha_m*(T^2 - 1) + beta_m*(T - 1) is the exact inverse
        % of T = (alpha_m - 1 + sqrt(1 + 2*alpha_m*theta))/alpha_m
        beta_m = 1 - alpha_m;
        % contractions of the wall rows with the other nodes, independent of
        % kv_w. Order matters: Tm(iwall_m) runs outward from the wall
        % (Tm(2:3) for finite differences), while theta(iwall_b) and
        % kv(iwall_b) run inward from the wall (theta(end-1:-1:end-2) and
        % kv(end-1:-1:end-2) for finite differences), matching the order
        % of the coefficients in grad_Tm_coeff, grad_Trans_coeff, grad_C_coeff
        flux_fixed = grad_Tm_coeff(2:end)*Tm_in(iwall_m) + ...
            grad_Trans_coeff(2:end)*th_in(iwall_b);
        carrier_fixed = grad_C_coeff(2:end)*kv_in(iwall_b);
        wp = [Pw, alpha_m, beta_m, flux_fixed, carrier_fixed];
        
        % warm start from the previous root
        kv_w = NaN;
        if isfinite(kv_w_guess)
            x = kv_w_guess;
            for kw = 1:n_warm
                [res, ~, ~, dres] = f_wall_mass_residual(x, wp);
                step = res/dres;
                x = x - step;
                if ~(isfinite(x) && x > kv_bracket(1) && x < kv_bracket(2))
                    break
                end
                if abs(step) <= warm_tol
                    kv_w = x;
                    break
                end
            end
        end
        if isfinite(kv_w)
            kv_w_guess = kv_w;
            [~, theta_w] = f_wall_mass_residual(kv_w, wp);
            return
        end
        
        % bracketed solve
        lo = kv_bracket(1);
        hi = kv_bracket(2);
        for kr = 1:n_rounds
            kv_grid = linspace(lo, hi, n_grid);
            r_grid = f_wall_mass_residual(kv_grid, wp);
            if kr == 1
                r_lo = r_grid(1);
                r_hi = r_grid(end);
                if ~(isfinite(r_lo) && isfinite(r_hi)) || sign(r_lo) == sign(r_hi)
                    error(['f_imr_nodal: the mass-transfer wall balance does ' ...
                        'not change sign on kv_w in [%g, 1 - %g] at t = %g ' ...
                        '(R = %g, P = %g, stored kv(end) = %g, residuals %g ' ...
                        'and %g). This indicates a physical or indexing ' ...
                        'problem in the state, not a solver problem.'], lo, ...
                        1 - hi, tw, Rw, Pw, kv_in(end), r_lo, r_hi);
                end
            end
            ic = find(sign(r_grid(1:end-1)) ~= sign(r_grid(2:end)), 1);
            lo = kv_grid(ic);
            hi = kv_grid(ic+1);
        end
        kv_w = (lo + hi)/2;
        for kp = 1:n_polish
            [res, ~, ~, dres] = f_wall_mass_residual(kv_w, wp);
            % stay inside the final grid interval
            kv_w = min(max(kv_w - res/dres, lo), hi);
        end
        kv_w_guess = kv_w;
        [~, theta_w] = f_wall_mass_residual(kv_w, wp);
    end
    
    % mass-transfer wall balance and its derivative as functions of the wall
    % vapor mass fraction kv_w (elementwise, kv_w may be a vector),
    % wp = [P, alpha_m, beta_m, flux_fixed, carrier_fixed]. It is the balance
    % of f_wall_mass_bc multiplied through by
    % (kv_w*Rva_diff + Rg_star)*Tw*(1 - kv_w) > 0, which has the same root
    % and no reciprocal that vanishes as kv_w -> 1
    function [r, theta_w, Tw, dr] = f_wall_mass_residual(kv_w, wp)
        Pw = wp(1);
        % vapor partial pressure, and the wall temperature from inverting
        % Pv(Tw) = C1_pv*exp(C2_pv/Tw) (C2_pv < 0)
        q = kv_w*Rvg_ratio + 1 - kv_w;
        ps = Pw*kv_w*Rvg_ratio./q;
        lg = log(ps/C1_pv);
        Tw = C2_pv./lg;
        theta_w = 0.5*wp(2)*(Tw.^2 - 1) + wp(3)*(Tw - 1);
        mix = kv_w*Rva_diff + Rg_star;
        denom = mix.*Tw.*(1 - kv_w);
        thermal = grad_Tm_coeff(1)*Tw + grad_Trans_coeff(1)*theta_w + wp(4);
        r = denom.*thermal + Pw*(grad_C_coeff(1)*kv_w + wp(5));
        if nargout > 3
            % d/d(kv_w) through ps, Tw and theta_w
            dps = Pw*Rvg_ratio./q.^2;
            dTw = -C2_pv./lg.^2.*dps./ps;
            dth = (wp(2)*Tw + wp(3)).*dTw;
            ddenom = Rva_diff*Tw.*(1 - kv_w) + mix.*dTw.*(1 - kv_w) - mix.*Tw;
            dthermal = grad_Tm_coeff(1)*dTw + grad_Trans_coeff(1)*dth;
            dr = ddenom.*thermal + denom.*dthermal + Pw*grad_C_coeff(1);
        end
    end
    
    % bubble-wall theta at the state X, for the Jacobian
    function theta_bw = f_theta_bw_of_X(t, X)
        th_X = X(ibubtherm);
        Tm_X = X(imedtherm);
        if masstrans
            kv_X = X(imass);
            theta_bw = f_wall_mass_closure(Tm_X, th_X, kv_X, X(3), X(1), t);
        else
            theta_bw = f_wall_theta_exact(Tm_X, th_X);
        end
    end
    
    % temperature boundary condition at the bubble wall, thermal only:
    % residual of the heat flux balance (theta is the shared field unless
    % given)
    function theta_w = f_bubble_wall_thermal_bc(theta_bw, Tm_in, th_in)
        Tm_trans = Tm_in(iwall_m);
        T_trans = th_in(iwall_b);
        Tw_prelim = inv_alpha*(alpha_g - 1 + sqrt(1+2*theta_bw*alpha_g));
        theta_w = grad_Tm_coeff * vertcat(Tw_prelim, Tm_trans) + ...
            grad_Trans_coeff * vertcat(theta_bw, T_trans);
    end
    
    % exact bubble-wall theta for the thermal-only balance. With
    % u = sqrt(1+2*alpha_g*theta_bw), Tw = (alpha_g-1+u)/alpha_g is linear
    % and theta_bw = (u^2-1)/(2*alpha_g) quadratic in u, so the residual is
    % c2*u^2 + 2*c1*u + q = 0. Both roots are formed, u < 0 is discarded
    % (u is a square root), and of two admissible roots the one closest to
    % theta_bw_guess is kept; the sign of c2 depends on the grid spacing
    % and the conductivity ratio, so no branch is hard-coded. c2 = 0 is the
    % linear case. fzero is the fallback when no admissible root exists.
    function theta_bw = f_wall_theta_exact(Tm_in, th_in)
        if nargin < 2
            th_in = theta;
        end
        c1 = grad_Tm_coeff(1);
        c2 = grad_Trans_coeff(1);
        K = grad_Tm_coeff(2:end)*Tm_in(iwall_m) + ...
            grad_Trans_coeff(2:end)*th_in(iwall_b);
        q = 2*c1*(alpha_g - 1) - c2 + 2*alpha_g*K;
        if abs(c2) < eps
            uroots = -q/(2*c1);
        else
            disc = c1^2 - c2*q;
            if disc >= 0
                % cancellation-free pair of roots
                sq = sqrt(disc);
                if c1 >= 0
                    qq = -(c1 + sq);
                else
                    qq = -(c1 - sq);
                end
                uroots = [qq/c2, q/qq];
            else
                uroots = [];
            end
        end
        uroots = uroots(isfinite(uroots) & uroots >= 0);
        if isempty(uroots)
            try
                theta_bw = fzero(@(x) f_bubble_wall_thermal_bc(x, Tm_in, th_in), ...
                    theta_bw_guess, foptions);
            catch
                theta_bw = NaN;
            end
            if ~isfinite(theta_bw)
                % no root exists (an unphysical trial state of the
                % integrator): take the least-residual point of the
                % quadratic so the step is rejected on its error estimate
                % rather than aborting the run
                if abs(c2) < eps
                    u = 0;
                else
                    u = max(-c1/c2, 0);
                end
                theta_bw = (u^2 - 1)/(2*alpha_g);
            end
        else
            thetas = (uroots.^2 - 1)/(2*alpha_g);
            [~,kmin] = min(abs(thetas - theta_bw_guess));
            theta_bw = thetas(kmin);
        end
    end
    
end
