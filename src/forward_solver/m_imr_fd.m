% file m_imr_fd.m
% brief contains module m_imr_fd

% brief This module features a fourth- and sixth-order accurate finite
% difference solver of the PDEs involving thermal transport and
% viscoelasticity to solve Rayleigh-Plesset equations
function varargout = m_imr_fd(varargin)
    
    % problem initialization
    [eqns_opts, solve_opts, init_opts, init_stress, tspan_opts, out_opts, ...
        acos_opts, wave_opts, sigma_opts, thermal_opts, mass_opts] ...
        = f_call_params(varargin{:});
    
    % defaults
    
    % non-Newtonian viscosity parameters
    fnu             = 0;
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
    % if (stress == 4)
    %     ptt = 1;
    % else
    %     ptt = 0;
    % end
    
    % solver options
    method          = solve_opts(1);
    spectral        = solve_opts(2);
    divisions       = solve_opts(3);
    Nv              = solve_opts(4);
    Nt              = solve_opts(5);
    Mt              = solve_opts(6);
    Lv              = solve_opts(7);
    Lt              = solve_opts(8);
    
    % dimensionless initial conditions
    Rzero           = init_opts(1);
    Uzero           = init_opts(2);
    Pb_star         = init_opts(3);
    P8              = init_opts(4);
    T8              = init_opts(5);
    Pv_star         = init_opts(6);
    Req             = init_opts(7);
    alphax          = init_opts(8);
    
    % dimensionaless initial stress
    Szero           = init_stress;
    
    % time span options
    tspan           = tspan_opts;
    tfin            = tspan(end);
    
    % output options
    dimensionalout  = out_opts(1);
    progdisplay     = out_opts(2);
    
    % physical parameters
    
    % acoustic parameters
    Cstar           = acos_opts(1);
    GAMa            = acos_opts(2);
    kappa           = acos_opts(3);
    nstate          = acos_opts(4);
    hugoniot_s      = acos_opts(5);
    
    % dimensionless waveform parameters
    om              = wave_opts(1);
    ee              = wave_opts(2);
    tw              = wave_opts(3);
    dt              = wave_opts(4);
    mn              = wave_opts(5);
    wave_type       = wave_opts(6);
    if wave_type < 0
        wave_poly   = wave_opts(7);
        wave_dpoly  = wave_opts(8);
    else
        wave_poly   = [];
        wave_dpoly  = [];
    end
    
    pvarargin = [om,ee,tw,dt,mn,wave_type,wave_poly,wave_dpoly];
    
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
    zeNO            = sigma_opts(12);
    iDRe            = sigma_opts(13);
    graded          = sigma_opts(14);
    Ca1             = sigma_opts(15);
    l1              = sigma_opts(16);
    l2              = sigma_opts(17);
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
    
    % dimensionaless mass transfer
    Fom             = mass_opts(1);
    C0              = mass_opts(2);
    Rv_star         = mass_opts(3);
    Ra_star         = mass_opts(4);
    L_heat_star     = mass_opts(5);
    % mv0             = mass_opts(6);
    % ma0             = mass_opts(7);
    
    % pre_process
    
    % creates finite difference matrices
    D_Matrix_T_C = f_finite_diff_mat(Nt,1,0);
    DD_Matrix_T_C = f_finite_diff_mat(Nt,2,0);
    D_Matrix_Tm = f_finite_diff_mat(Mt,1,1);
    DD_Matrix_Tm = f_finite_diff_mat(Mt,2,1);
    
    % inside the bubble
    N = Nt-1;
    deltaY = 1/N;
    i = 1:1:N+1;
    y = ((i-1)*deltaY)';
    
    % outside the bubble
    Nm = Mt-1;
    deltaYm = -2/Nm;
    j = 1:1:Nm+1;
    xi = (1+(j-1)*deltaYm)';
    yT = ((2./(xi+1)-1)*Lt+1);
    yT2 = yT.^2;
    yT3 = yT.^3;
    iyT3 = yT.^-3;
    iyT4 = yT.^-4;
    iyT6 = yT.^-6;
    
    % precomputations move these to f_call_params
    %LDR = LAM*De/Re8;
    sam = 1 + GAMa;
    no = (nstate-1)/nstate;
    nog = (nstate-1)/2;
    kapover = (kappa-1)/kappa;
    C1_pv = 1.17e11/P8;
    C2_pv = -5200/T8;
    % boundary conditions coefficients
    coeff = [-1.5 , 2 ,-0.5 ];
    Rva_ratio = Rv_star / Ra_star;
    inv_alpha = 1/alpha_g;
    Rva_diff = Rv_star - Ra_star;
    grad_Tm_coeff = 2*chi*iota/deltaYm*coeff;
    grad_Trans_coeff = -coeff*chi/deltaY;
    grad_C_coeff = -coeff*Fom*L_heat_star/deltaY;
    
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
    if bubtherm == 0 && masstrans
        Nt = 1;
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
    cdd = preStressInt(Lv,Nv);
    % stress spectra index management
    ivisco1 = (4+Nt+Mt+Nc):(3+Nt+Mt+Nc+Nv);
    ivisco2 = (4+Nt+Mt+Nc+Nv):(3+Nt+Mt+Nc+2*Nv);
    
    % initial condition assembly
    
    % bubble temperature initial condition
    if bubtherm
        Tau0 = zeros(Nt,1);
    elseif masstrans
        Tau0 = zeros(1,1);
    else
        Tau0 = [];
    end
    
    % medium initial temperature
    if medtherm
        Tm0 = ones(Mt,1);
    else
        Tm0 = [];
    end
    
    % mass transfer initial condition
    if masstrans
        C0vec = C0*ones(Nc,1);
    else
        C0vec = [];
    end
    
    % initial condition vector
    init = [Rzero;
    Uzero;
    Pb_star;
    Tau0;
    Tm0;
    C0vec;
    Szero];
    
    guess = -0.0001;
    foptions = optimset('TolFun',1e-10);
    
    % solver start
    f_display(radial, bubtherm, medtherm, masstrans, stress, spectral,...
        nu_model, eps3, Pv_star, Re8, De, Ca, LAM, 'finite difference');
    bubble = @SVBDODE;
    [t,X] = f_odesolve(bubble, init, method, divisions, tspan, tfin);
    
    % extract result
    R    = X(:,1);
    Rdot = X(:,2);
    P    = X(:,3);
    if bubtherm
        Tau = X(:,ibubtherm);
        T = TW(Tau);
        if medtherm
            Tm = X(:,imedtherm);
        end
    end
    if masstrans
        C = X(:,imass);
    end
    
    % transform variables back into their dimensional form
    if (dimensionalout == 1)
        t = t*tc;
        R = R*Rref;
        Rdot = Rdot*uc;
        P = P*p0;
        if bubtherm
            T = T*T8;
        end
        if medtherm
            Tm = Tm*T8;
        end
    end
    
    % outputs assembly
    varargout{1} = t;
    varargout{2} = R;
    varargout{3} = Rdot;
    varargout{4} = P;
    if bubtherm
        varargout{5} = T;
    end
    if bubtherm == 1 && medtherm == 1
        varargout{6} = Tm;
    else
        varargout{6} = ((T8 - 1)*dimensionalout + 1)*ones(size(t,1),1);
    end
    if masstrans == 1
        varargout{7} = C;
    end
    
    % solver function
    function [dXdt] = SVBDODE(t,X)
        
        % showing output
        if progdisplay
            disp(t/tfin);
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
            Tau = X(ibubtherm);
        end
        if masstrans
            C = X(imass);
        end
        
        % empty output variables
        Taudot = [];
        Tmdot = [];
        Cdot = [];
        
        % boundary condition evaluation
        if medtherm && masstrans
            % temperature in the material
            Tm = X(imedtherm);
            prelim = fzero(@Boundary_full,guess,foptions);
            guess = prelim;
        elseif medtherm
            % temperature in the material
            Tm = X(imedtherm);
            prelim = fzero(@Boundary_medtherm,guess,foptions);
            guess = prelim;
        elseif masstrans
            prelim = fzero(@Boundary_mass,guess,foptions);
            guess = prelim;
        else
            prelim = 0;
        end
        
        if bubtherm || masstrans
            Tau(end) = prelim;
            T = TW(Tau);
        end
        
        % equations of motion
        
        if bubtherm && masstrans
            % vapor concentration at the boundary
            C(end) =  CW(T(end),P);
            
            % temperature field of the gas inside the bubble
            DTau  = D_Matrix_T_C*Tau;
            DDTau = DD_Matrix_T_C*Tau;
            % concentration of vapor inside the bubble
            DC  = D_Matrix_T_C*C;
            DDC = DD_Matrix_T_C*C;
            Rmix = C*Rv_star + (1-C)*Ra_star;
            RDC = (Rva_diff./Rmix).*DC;
            
            % internal pressure equation
            Pdot = 3/R*(chi*(kappa-1)*DTau(end)/R - kappa*P*Rdot +...
                + kappa*P*Fom*Rv_star*DC(end)/(T(end)*R*Rmix(end)*(1-C(end))));
            
            % internal bubble velocity
            U_vel = (chi/R*(kappa-1).*DTau-y*R*Pdot/3)/(kappa*P) + Fom/R*RDC;
            
            % calculate mixture thermal conductivity
            Kstar_g = alpha_g*T+beta_g;
            Kstar_v = alpha_v*T+beta_v;
            Kstar = C.*Kstar_v + (1-C).*Kstar_g;
            
            % temperature of the gas inside the bubble
            nonlinear_term = (chi*DDTau./R^2+Pdot).*(kapover*Kstar.*T/P);
            advection_term = -DTau.*(U_vel-y*Rdot)./R;
            mass_diffusion = (Fom/(R^2)).*(Rva_diff./Rmix).*DC.*DTau;
            Taudot = advection_term + nonlinear_term + mass_diffusion;
            Taudot(end) = 0;
            
            % vapor equations inside the bubble
            nonlinear_diffusion = DC.*(DTau./(Kstar.*T)+RDC);
            advection_term =  (U_vel-Rdot.*y)/R.*DC;
            Cdot = Fom/R^2*(DDC - nonlinear_diffusion) - advection_term;
            Cdot(end) = 0;
            
        elseif bubtherm
            % calculate mixture thermal conductivity
            Kstar = alpha_g*T+beta_g;
            
            % temperature gradients
            DTau  = D_Matrix_T_C*Tau;
            DDTau = DD_Matrix_T_C*Tau;
            
            % internal pressure equation
            Pdot = 3/R*(chi*(kappa-1)*DTau(end)/R - kappa*P*Rdot);
            
            % internal bubble velocity
            U_vel = (chi/R*(kappa-1).*DTau-y*R*Pdot/3)/(kappa*P);
            
            % temperature of the gas inside the bubble
            diffusion_term = (chi*DDTau./R^2+Pdot).*(kapover*Kstar.*T/P);
            advection_term = -DTau.*(U_vel-y*Rdot)./R;
            Taudot = advection_term + diffusion_term;
            Taudot(end) = 0;
            
        elseif masstrans
            % calculate mixture thermal conductivity
            Kstar = alpha_g*T+beta_g;
            
            % computing mass transfer at the wall
            C(end) = CW(T,P);
            % calculate mass transfer variable
            Rmix = C*Rv_star + (1-C)*Ra_star;
            % concentration of vapor inside the bubble
            DC  = D_Matrix_T_C*C;
            DDC = DD_Matrix_T_C*C;
            RDC = (Rva_diff./Rmix).*DC;
            
            % internal pressure equation
            Pdot = 3/R*(-kappa*P*Rdot + ...
                kappa*P*Fom*Rv_star*DC(end)/(T*R*Rmix(end)*(1-C(end))));
            
            % internal bubble velocity
            U_vel = (-y*R*Pdot/3)/(kappa*P) + Fom/R*RDC;
            
            % temperature of the gas inside the bubble
            Taudot = Pdot.*(kapover*Kstar.*T/P);
            
            % vapor concentration equation
            nonlinear_diffusion = RDC.*DC;
            advection =  (U_vel-Rdot.*y)/R.*DC;
            Cdot = Fom/R^2*(DDC - nonlinear_diffusion) - advection;
            Cdot(end) = 0;
        else
            % polytropic gas
            Pdot = -3*kappa*Rdot/R*P;
        end
        
        % updating the viscous forces/Reynolds number
        if nu_model ~= 0
            [fnu,intfnu,dintfnu,ddintfnu] = f_viscosity(nu_model,Rdot, ...
                R,v_a,v_nc,v_lambda_star);
        end
        
        if medtherm
            % boundary temperature
            Tm(1) = T(end);
            % temperature gradients
            DTm = D_Matrix_Tm*Tm;
            DDTm = DD_Matrix_Tm*Tm;
            % material temperature equations
            advection = (1+xi).^2./(Lt*R).*(Rdot./yT2.*(1-yT3)/2 +...
                Foh/R.*((xi+1)/(2*Lt)-1./yT)).*DTm;
            diffusion = Foh/R^2.*(xi+1).^4/Lt^2.*DDTm/4;
            % stress dissipation
            [taugradu] = f_stress_dissipation(stress,spectral,Req,R,Rdot, ...
                Ca,Br,Re8,alphax,yT2,yT3,iyT3,iyT4,iyT6,X,ZZT,ivisco1, ...
                ivisco2,fnu,DRe);
            Tmdot = advection + diffusion + taugradu;
            % sets boundary condition on temperature
            Tmdot(1) = 0;
            Tmdot(end) = 0;
        end
        
        % pressure waveform
        [Pf8,Pf8dot] = f_pinfinity(t,pvarargin);
        
        % stress equation evolution
        if graded
            [S,Sdot,Z1dot,Z2dot] = f_stress_calc_graded(stress,Req,R,Ca,Ca1,...
                Re8,Rdot,alphax,intfnu,dintfnu,iDRe,l1,l2,v_a,v_nc);
        else
            [S,Sdot,Z1dot,Z2dot] = f_stress_calc(stress,X,Req,R,Ca,De,Re8, ...
                Rdot,alphax,ivisco1,ivisco2,LAM,zeNO,cdd,intfnu,dintfnu,iDRe);
        end
        
        % bubble wall acceleration
        [Rddot] = f_radial_eq(radial, P, Pdot, Pf8, Pf8dot, iWe, R, Rdot, S, ...
            Sdot, Cstar, sam, no, GAMa, nstate, nog, hugoniot_s, JdotA, ...
            ddintfnu, iDRe);
        
        % output assembly
        dXdt = [Rdot;
        Rddot;
        Pdot;
        Taudot;
        Tmdot;
        Cdot;
        Z1dot;
        Z2dot];
        
    end
    % end of solver
    
    % functions called by solver
    
    % calculates the temperature at the bubble wall as a function of \tau
    function Tw = TW(Tauw)
        Tw = (alpha_g -1 + sqrt(1+2*Tauw*alpha_g)) / alpha_g;
    end
    
    % calculates the concentration at the bubble wall
    function Cw = CW(Tw,P)
        theta = Rv_star/Ra_star*(P./(f_pvsat(Tw*T8)/P8)-1);
        Cw = 1./(1+theta);
    end
    
    % solves temperature boundary conditions at the bubble wall
    function Tauw = Boundary_full(prelim)
        
        Tm_trans = Tm(2:3);
        T_trans = Tau(end-1:-1:end-2);
        
        Tw_prelim = inv_alpha*(alpha_g - 1 + sqrt(1+2*prelim*alpha_g));
        Pv_prelim = C1_pv*exp(C2_pv/Tw_prelim);
        demCw = (1 + (Rva_ratio) * (P ./ Pv_prelim - 1));
        Cw_prelim = 1 / demCw;
        
        C_trans = C(end-1:-1:end-2);
        Tauw = grad_Tm_coeff * [Tw_prelim; Tm_trans]  + ...
            grad_Trans_coeff * [prelim ;T_trans] + ...
            grad_C_coeff * P * ((Cw_prelim * Rva_diff + Ra_star))^-1 * ...
            (Tw_prelim * (1 - Cw_prelim))^-1 * [Cw_prelim;
        C_trans];
        
    end
    
    % solves temperature boundary conditions at the bubble wall
    function Tauw = Boundary_medtherm(prelim)
        
        Tm_trans = Tm(2:3);
        T_trans = Tau(end-1:-1:end-2);
        
        Tw_prelim = inv_alpha*(alpha_g - 1 + sqrt(1+2*prelim*alpha_g));
        
        Tauw = grad_Tm_coeff * [Tw_prelim; Tm_trans]  + ...
            grad_Trans_coeff * [prelim ;
        T_trans];
        
    end
    
    % solves temperature boundary conditions at the bubble wall
    function Tauw = Boundary_mass(prelim)
        
        % no thermal transport in the surroundings
        Tm_trans = ones(2,1);
        T_trans = zeros(2,1);
        
        Tw_prelim = inv_alpha*(alpha_g - 1 + sqrt(1+2*prelim*alpha_g));
        Pv_prelim = C1_pv*exp(C2_pv/Tw_prelim);
        demCw = (1 + (Rva_ratio) * (P ./ Pv_prelim - 1));
        Cw_prelim = 1 / demCw;
        
        C_trans = C(end-1:-1:end-2);
        Tauw = grad_Tm_coeff * [Tw_prelim; Tm_trans]  + ...
            grad_Trans_coeff * [prelim ;T_trans] + ...
            grad_C_coeff * P * ((Cw_prelim * Rva_diff + Ra_star))^-1 * ...
            (Tw_prelim * (1 - Cw_prelim))^-1 * [Cw_prelim;
        C_trans];
        
    end
    
    % Creates finite difference matrices
    % Nodes: Number of nodes
    % order: order of differentiation ( 1st derivative vs 2nd derivative)
    % Tm_check: 0 not for external temp , 1 used for external temp
    function [Diff_Matrix] = f_finite_diff_mat(Nodes,order,Tm_check)
        
        % coordinate creation
        if Tm_check == 0
            N = Nodes-1;
            deltaY = 1/N;
            K = 1:1:N+1;
            y = (K-1)*deltaY;
        elseif Tm_check == 1
            N = Nodes-1;
            deltaY = -2/N;
            K = 1:1:N+1;
            y = 1+(K-1)*deltaY;
        end
        
        Diff_Matrix = zeros(Nodes);
        
        if order == 1
            % in between
            for counter = 2:N
                Diff_Matrix(counter,counter+1) = 0.5 ;
                Diff_Matrix(counter,counter-1) = -0.5 ;
            end
            if Tm_check == 0
                Diff_Matrix(end,end) = 1.5;
                Diff_Matrix(end,end-1) = -2;
                Diff_Matrix(end,end-2) = 0.5;
            end
            Diff_Matrix = Diff_Matrix / deltaY ;
            
        elseif order == 2
            if Tm_check == 0
                % in between
                for counter = 2:N
                    Diff_Matrix(counter,counter+1) = 1 + deltaY/y(counter);
                    Diff_Matrix(counter,counter)   = -2;
                    Diff_Matrix(counter,counter-1) = 1 - deltaY/y(counter);
                end
            elseif Tm_check == 1
                for counter=2:N  %in between
                    Diff_Matrix(counter,counter+1) = 1 ;
                    Diff_Matrix(counter,counter)   = -2 ;
                    Diff_Matrix(counter,counter-1) = 1 ;
                end
            end
            if Tm_check == 0
                Diff_Matrix(1,1)= -6;
                Diff_Matrix(1,2) = 6;
            end
            Diff_Matrix = Diff_Matrix / (deltaY^2) ;
        end
        
        Diff_Matrix = sparse(Diff_Matrix);
        
    end
    
    function cdd = preStressInt(Lt,N)
        cdd = Lt*N*0;
    end
    
end
