% file f_imr_fd.m
% brief contains function f_imr_fd

% brief This function features a second-order accurate finite difference
% solver of the PDEs involving thermal transport and viscoelasticity
% coupled to the radial (Rayleigh-Plesset family) dynamics
function varargout = f_imr_fd(varargin)
    
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
    % Lv            = solve_opts(7), only used by the spectral solver
    Lt              = solve_opts(8);
    rtol            = solve_opts(9);
    atol            = solve_opts(10);
    jacobian        = solve_opts(11);
    
    % dimensionless initial conditions
    Rzero           = init_opts(1);
    Rdotzero        = init_opts(2);
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
    
    % output options
    dimensionalout  = out_opts(1);
    progdisplay     = out_opts(2);
    tref            = out_opts(3);
    Rref            = out_opts(4);
    Uref            = out_opts(5);
    
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
    kv0             = mass_opts(2);
    Rv_star         = mass_opts(3);
    Rg_star         = mass_opts(4);
    L_heat_star     = mass_opts(5);
    
    % pre_process
    
    % creates finite difference matrices
    dmatrix_Tkv  = f_finite_diff_mat(Nt,1,0);
    ddmatrix_Tkv = f_finite_diff_mat(Nt,2,0);
    dmatrix_Tm   = f_finite_diff_mat(Mt,1,medtherm);
    ddmatrix_Tm  = f_finite_diff_mat(Mt,2,medtherm);
    
    % inside the bubble
    N = Nt-1;
    deltaY = 1/N;
    i = 1:1:N+1;
    y = ((i-1)*deltaY)';
    
    % outside the bubble
    Nm = Mt-1;
    deltaYm = -2/Nm;
    j = 1:1:Nm+1;
    xi   = (1+(j-1)*deltaYm)';
    yT   = ((2./(xi+1)-1)*Lt+1);
    yT2  = yT.^2;
    yT3  = yT.^3;
    iyT3 = yT.^-3;
    iyT4 = yT.^-4;
    iyT6 = yT.^-6;
    
    % precomputations move these to f_call_params
    %LDR = LAM*De/Re8;
    sam = 1 + GAMa;
    no = (nstate-1)/nstate;
    nog = (nstate-1)/2;
    kapover = (kappa-1)/kappa;
    [~, pv_c1, pv_c2] = f_pvsat(T8);
    C1_pv = pv_c1/P8;
    C2_pv = pv_c2/T8;
    % boundary conditions coefficients
    coeff = [-1.5 , 2 ,-0.5 ];
    Rvg_ratio = Rv_star / Rg_star;
    Rva_diff = Rv_star - Rg_star;
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
    if medtherm == 0
        Mt = 0;
    end
    ibubtherm   = 4:(3+Nt);
    imedtherm   = (4+Nt):(3+Nt+Mt);
    imass       = (4+Nt+Mt):(3+Nt+Mt+Nc);
    
    % precomputations for viscous dissipation
    % zT = 1 - 2./(1 + (yT - 1)/Lv);
    ZZT = 0;
    % stress integral precomputation, only used by the spectral solver
    cdd = 0;
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
    
    theta_bw_guess = -0.0001;
    foptions = optimset('TolFun',1e-12);
    
    % solver start
    f_display(radial, bubtherm, medtherm, masstrans, stress, spectral,...
        nu_model, eps3, Pv_star, Re8, De, Ca, LAM, 'finite difference');
    bubble = @SVBDODE;
    % Jacobian for the stiff integrators: analytic where implemented,
    % otherwise the sparsity pattern
    if jacobian == 2 && bubtherm && ~masstrans
        Sp = @f_jacobian;
    elseif jacobian >= 1
        Sp = f_jpattern(numel(init), ibubtherm, imedtherm, imass, ...
            (4+Nt+Mt+Nc):numel(init));
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
        T = f_theta_of_T(theta,kv0);
        if medtherm
            Tm = X(:,imedtherm);
        end
    end
    if masstrans
        kv = X(:,imass);
        kv(:,end) = f_kv_of_T(T(:,end),P);
        T = f_theta_of_T(theta,kv);
    end
    
    % transform variables back into their dimensional form
    if (dimensionalout == 1)
        t = t*tref;
        R = R*Rref;
        Rdot = Rdot*Uref;
        P = P*P8;
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
    else
        varargout{5} = [];
    end
    if medtherm == 1
        varargout{6} = Tm;
    else
        varargout{6} = ((T8 - 1)*dimensionalout + 1)*ones(size(t,1),1);
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
            % warm-started Newton on the transcendental wall balance
            theta_bw = f_wall_theta_newton(theta_bw_guess);
            theta_bw_guess = theta_bw;
            theta(end) = theta_bw;
        elseif medtherm
            % temperature in the material
            Tm = X(imedtherm);
            % closed-form wall temperature, no root finding
            theta_bw = f_wall_theta_exact(Tm,theta);
            theta(end) = theta_bw;
        end
        
        % bubble temperature
        if masstrans
            T = f_theta_of_T(theta,kv);
        elseif bubtherm
            T = f_theta_of_T(theta,kv0);
        end
        
        % bubble equations of motion
        if bubtherm && masstrans
            % vapor concentration at the boundary
            kv(end) = f_kv_of_T(T(end),P);
            
            % temperature field gradients inside the bubble
            dtheta  = dmatrix_Tkv*theta;
            ddtheta = ddmatrix_Tkv*theta;
            
            % vapor concentration gradients inside the bubble
            dkv     = dmatrix_Tkv*kv;
            ddkv    = ddmatrix_Tkv*kv;
            Rmix    = kv*Rv_star + (1-kv)*Rg_star;
            RDkv    = (Rva_diff./Rmix).*dkv;
            
            % bubble pressure equation
            Pdot = 3/R*(chi*(kappa-1)*dtheta(end)/R - kappa*P*Rdot +...
                + kappa*P*Fom*Rv_star*dkv(end)/(T(end)*R*Rmix(end)*(1-kv(end))));
            
            % bubble wall velocity equation
            Uvel = (chi/R*(kappa-1).*dtheta-y*R*Pdot/3)/(kappa*P) + Fom/R*RDkv;
            
            % mixture thermal conductivity
            alpha_mix  = kv.*alpha_v + (1-kv).*alpha_g;
            Kstar_g = alpha_g*T+beta_g;
            Kstar_v = alpha_v*T+beta_v;
            Kstar   = kv.*Kstar_v + (1-kv).*Kstar_g;
            
            % temperature evolution inside the bubble
            nonlinear_term = (chi*ddtheta./R^2+Pdot).*(kapover*Kstar.*T/P);
            advection_term = -dtheta.*(Uvel-y*Rdot)/R;
            mass_diffusion = (Fom/(R^2)).*(Rva_diff./Rmix).*dkv.*dtheta;
            thetadot = advection_term + nonlinear_term + mass_diffusion;
            % temperature at bubble wall, handled by flux boundary bc
            thetadot(end) = 0;
            
            % water vapor evolution inside the bubble
            nonlinear_diffusion = ...
                dkv.*(dtheta./(sqrt(1+2*alpha_mix.*theta).*T)+RDkv);
            advection_term =  (Uvel-Rdot.*y)/R.*dkv;
            kvdot = Fom/R^2*(ddkv - nonlinear_diffusion) - advection_term;
            % concentration at bubble wall, handled by flux boundary bc
            kvdot(end) = 0;
            
        elseif bubtherm
            % mixture thermal conductivity
            Kstar_g = alpha_g*T+beta_g;
            Kstar_v = alpha_v*T+beta_v;
            Kstar   = kv0.*Kstar_v + (1-kv0).*Kstar_g;
            
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
            [fnu,intfnu,dintfnu,ddintfnu] = f_viscosity(nu_model,Rdot, ...
                R,v_a,v_nc,v_lambda_star);
        end
        
        % surrounding temperature
        if medtherm
            % boundary temperature
            Tm(1) = T(end);
            % temperature gradients
            dTm = dmatrix_Tm*Tm;
            ddTm = ddmatrix_Tm*Tm;
            % material temperature equations
            advection = (1+xi).^2./(Lt*R).*(Rdot./yT2.*(1-yT3)/2 +...
                Foh/R.*((xi+1)/(2*Lt)-1./yT)).*dTm;
            diffusion = Foh/R^2.*(xi+1).^4/Lt^2.*ddTm/4;
            % stress dissipation
            [taugradu] = f_stress_dissipation(stress,spectral,Req,R,Rdot, ...
                Ca,Br,Re8,alphax,yT2,yT3,iyT3,iyT4,iyT6,X,ZZT,ivisco1, ...
                ivisco2,fnu,DRe);
            % surrounding temperature evolution equation
            Tmdot = advection + diffusion + taugradu;
            % sets boundary condition on temperature
            Tmdot(1) = 0;
            Tmdot(end) = 0;
        end
        
        % pressure waveform
        [Pf8,Pf8dot] = f_pinfinity(t,pvarargin);
        
        % stress equation evolution
        [S,Sdot,Z1dot,Z2dot] = f_stress(stress,X,Req,R,Ca,De,Re8, ...
            Rdot,alphax,ivisco1,ivisco2,LAM,zeNO,cdd,intfnu,dintfnu,iDRe);
        
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
    
    % analytic Jacobian of SVBDODE for bubtherm = 1, masstrans = 0. The
    % thermal field blocks and the bubble pressure row are assembled in
    % closed form from the finite difference operators. The radial and
    % stress closure rows are chained through Pdot, and the partials of
    % that low-dimensional closure (radial model, stress model, viscosity)
    % are taken by central differences of the closure alone, which is far
    % cheaper than differencing the full right-hand side.
    function J = f_jacobian(t,X)
        
        nX = numel(X);
        jR = X(1);
        jRdot = X(2);
        jP = X(3);
        ith = ibubtherm;
        iz = (4+Nt+Mt+Nc):nX;
        nth = numel(ith);
        
        % effective bubble temperature (wall value set by the boundary
        % condition) and its sensitivity E = d(theta_eff)/dX
        jtheta = X(ith);
        E = sparse(1:nth, ith, 1, nth, nX);
        if medtherm
            jTm = X(imedtherm);
            [theta_bw, ~, dth_dK] = f_wall_theta_exact(jTm,jtheta);
            jtheta(end) = theta_bw;
            % last row: theta_bw depends on Tm(2:3), theta(end-1), theta(end-2)
            E(nth,ith(end)) = 0;
            E(nth,ith(end-1)) = dth_dK*grad_Trans_coeff(2);
            E(nth,ith(end-2)) = dth_dK*grad_Trans_coeff(3);
            E(nth,imedtherm(2)) = dth_dK*grad_Tm_coeff(2);
            E(nth,imedtherm(3)) = dth_dK*grad_Tm_coeff(3);
        end
        
        % bubble temperature field quantities
        alpha_m = kv0*alpha_v + (1-kv0)*alpha_g;
        sq = sqrt(1+2*alpha_m*jtheta);
        jT = (alpha_m - 1 + sq)/alpha_m;
        dT = 1./sq;
        Ka = kv0*alpha_v + (1-kv0)*alpha_g;
        Kb = kv0*beta_v + (1-kv0)*beta_g;
        Kst = Ka*jT + Kb;
        dth = dmatrix_Tkv*jtheta;
        ddth = ddmatrix_Tkv*jtheta;
        
        % bubble pressure and its partials
        cP = 3*chi*(kappa-1);
        jPdot = cP*dth(end)/jR^2 - 3*kappa*jP*jRdot/jR;
        Pdot_R = -2*cP*dth(end)/jR^3 + 3*kappa*jP*jRdot/jR^2;
        Pdot_Rdot = -3*kappa*jP/jR;
        Pdot_P = -3*kappa*jRdot/jR;
        % d(Pdot)/dX, through theta_eff
        Pdot_th = (cP/jR^2)*dmatrix_Tkv(end,:);
        dPdot = Pdot_th*E;
        dPdot(1) = dPdot(1) + Pdot_R;
        dPdot(2) = dPdot(2) + Pdot_Rdot;
        dPdot(3) = dPdot(3) + Pdot_P;
        
        % bubble temperature rows
        cu = chi*(kappa-1)/kappa;
        Uv = (cu*dth/jR - y*jR*jPdot/(3*kappa))/jP;
        av = Uv - y*jRdot;
        W = kapover*Kst.*jT/jP;
        dW = kapover*(2*Ka*jT + Kb).*dT/jP;
        dfc = chi*ddth/jR^2 + jPdot;
        % coefficient of dPdot in each row
        cPd = W + dth.*y/(3*kappa*jP);
        % d(thetadot)/d(theta_eff), excluding the Pdot path
        Jtt = -spdiags(av/jR + cu*dth/(jR^2*jP), 0, nth, nth)*dmatrix_Tkv ...
            + spdiags(chi*W/jR^2, 0, nth, nth)*ddmatrix_Tkv ...
            + spdiags(dfc.*dW, 0, nth, nth);
        % direct partials with respect to R, Rdot, P
        dUv_R = -cu*dth/(jR^2*jP) - y*jPdot/(3*kappa*jP);
        th_R = -(dth/jR).*dUv_R + dth.*av/jR^2 - 2*chi*ddth.*W/jR^3;
        th_Rdot = dth.*y/jR;
        dUv_P = -cu*dth/(jR*jP^2) + y*jR*jPdot/(3*kappa*jP^2);
        th_P = -(dth/jR).*dUv_P - dfc.*W/jP;
        Jth = Jtt*E + cPd*dPdot;
        Jth(:,1) = Jth(:,1) + th_R;
        Jth(:,2) = Jth(:,2) + th_Rdot;
        Jth(:,3) = Jth(:,3) + th_P;
        % wall value held by the boundary condition
        Jth(end,:) = 0;
        
        J = sparse(nX,nX);
        J(1,2) = 1;
        J(3,:) = dPdot;
        J(ith,:) = Jth;
        
        % material temperature rows
        if medtherm
            ntm = numel(imedtherm);
            jTm(1) = jT(end);
            F = sparse(1:ntm, imedtherm, 1, ntm, nX);
            F(1,:) = dT(end)*E(nth,:);
            dTm = dmatrix_Tm*jTm;
            ddTm = ddmatrix_Tm*jTm;
            g1 = (1-yT3)./(2*yT2);
            g2 = (xi+1)/(2*Lt) - 1./yT;
            cA = (1+xi).^2/Lt;
            Am = cA.*(jRdot*g1/jR + Foh*g2/jR^2);
            Am_R = cA.*(-jRdot*g1/jR^2 - 2*Foh*g2/jR^3);
            Am_Rdot = cA.*g1/jR;
            Bm = Foh*(xi+1).^4/(4*Lt^2);
            Jmm = spdiags(Am, 0, ntm, ntm)*dmatrix_Tm + ...
                spdiags(Bm/jR^2, 0, ntm, ntm)*ddmatrix_Tm;
            % dissipation partials by central differences
            hR = eps^(1/3)*max(abs(jR),1);
            hU = eps^(1/3)*max(abs(jRdot),1);
            tg_R = (f_taugradu(jR+hR,jRdot,X) - f_taugradu(jR-hR,jRdot,X))/(2*hR);
            tg_U = (f_taugradu(jR,jRdot+hU,X) - f_taugradu(jR,jRdot-hU,X))/(2*hU);
            Jm = Jmm*F;
            Jm(:,1) = Jm(:,1) + Am_R.*dTm - 2*Bm.*ddTm/jR^3 + tg_R;
            Jm(:,2) = Jm(:,2) + Am_Rdot.*dTm + tg_U;
            Jm([1 end],:) = 0;
            J(imedtherm,:) = Jm;
        end
        
        % radial and stress closure rows: g(R, Rdot, P, Pdot, Z)
        v0 = [jR; jRdot; jP; jPdot; X(iz)];
        nv = numel(v0);
        G = zeros(1+numel(iz), nv);
        for k = 1:nv
            hk = eps^(1/3)*max(abs(v0(k)),1);
            vp = v0;
            vm = v0;
            vp(k) = vp(k) + hk;
            vm(k) = vm(k) - hk;
            G(:,k) = (f_closure(t,vp,X) - f_closure(t,vm,X))/(2*hk);
        end
        % chain rule through Pdot
        Jc = G(:,4)*dPdot;
        Jc(:,1:3) = Jc(:,1:3) + G(:,1:3);
        Jc(:,iz) = Jc(:,iz) + G(:,5:end);
        J([2 iz],:) = Jc;
        
    end
    
    % medium dissipation as a function of (R, Rdot)
    function tg = f_taugradu(Rv,Rdotv,X)
        fnu_loc = 0;
        if nu_model ~= 0
            fnu_loc = f_viscosity(nu_model,Rdotv,Rv,v_a,v_nc,v_lambda_star);
        end
        tg = f_stress_dissipation(stress,spectral,Req,Rv,Rdotv, ...
            Ca,Br,Re8,alphax,yT2,yT3,iyT3,iyT4,iyT6,X,ZZT,ivisco1, ...
            ivisco2,fnu_loc,DRe);
        tg = tg.*ones(size(yT));
    end
    
    % radial acceleration and stress evolution as a function of
    % v = [R; Rdot; P; Pdot; Z], mirroring SVBDODE
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
        [Pf8_loc,Pf8dot_loc] = f_pinfinity(t,pvarargin);
        [S_loc,Sdot_loc,Z1dot_loc,Z2dot_loc] = f_stress(stress,X,Req,Rv, ...
            Ca,De,Re8,Rdotv,alphax,ivisco1,ivisco2,LAM,zeNO,cdd, ...
            intf_loc,dintf_loc,iDRe);
        Rddot_loc = f_radial_eq(radial, v(3), v(4), Pf8_loc, Pf8dot_loc, ...
            iWe, Rv, Rdotv, S_loc, Sdot_loc, Cstar, sam, no, GAMa, nstate, ...
            nog, hugoniot_s, JdotA, ddintf_loc, iDRe);
        out = [Rddot_loc; Z1dot_loc; Z2dot_loc];
    end
    
    % solver functions
    
    % temperature at the bubble wall as a function of theta
    function Tw = f_theta_of_T(theta_w,kv)
        alpha_m  = kv.*alpha_v + (1-kv).*alpha_g;
        Tw = (alpha_m - 1 + sqrt(1+2*theta_w.*alpha_m)) ./ alpha_m;
    end
    
    % concentration at the bubble wall
    function kv_w = f_kv_of_T(Tw,pressure)
        theta_var = Rvg_ratio*(pressure./(f_pvsat(Tw*T8)/P8)-1);
        kv_w = 1./(1+theta_var);
    end
    
    % temperature boundary conditions at the bubble wall, full; the
    % optional second output is the derivative with respect to theta_bw
    function [theta_w, dtheta_w] = f_bubble_wall_full_bc(theta_bw)
        
        Tm_trans = Tm(2:3);
        T_trans = theta(end-1:-1:end-2);
        alpha_m  = kv(end)*alpha_v + (1-kv(end))*alpha_g;
        
        u_w = sqrt(1+2*theta_bw*alpha_m);
        Tw_prelim = (alpha_m - 1 + u_w)/alpha_m;
        Pv_prelim = C1_pv*exp(C2_pv/Tw_prelim);
        demCw = (1 + Rvg_ratio*(P/Pv_prelim - 1));
        kvw_prelim = 1 / demCw;
        
        kv_trans = kv(end-1:-1:end-2);
        Mw = kvw_prelim * Rva_diff + Rg_star;
        Nw = Tw_prelim * (1 - kvw_prelim);
        Gw = P / (Mw * Nw);
        Hw = grad_C_coeff * [kvw_prelim; kv_trans];
        theta_w = grad_Tm_coeff * [Tw_prelim; Tm_trans]  + ...
            grad_Trans_coeff * [theta_bw ;T_trans] + Gw * Hw;
        
        if nargout > 1
            dTw = 1/u_w;
            dPv = -Pv_prelim*C2_pv/Tw_prelim^2*dTw;
            dkvw = Rvg_ratio*P/(Pv_prelim*demCw)^2*dPv;
            dNw = dTw*(1 - kvw_prelim) - Tw_prelim*dkvw;
            dGw = -Gw*(Rva_diff*dkvw/Mw + dNw/Nw);
            dtheta_w = grad_Tm_coeff(1)*dTw + grad_Trans_coeff(1) + ...
                dGw*Hw + Gw*grad_C_coeff(1)*dkvw;
        end
        
    end
    
    % temperature boundary condition at the bubble wall, thermal only. The
    % flux balance c1*Tw + c2*theta_w + K = 0 with Tw = (ag-1+u)/ag and
    % theta_w = (u^2-1)/(2*ag), u = sqrt(1+2*ag*theta_w), is quadratic in u:
    % c2*u^2 + 2*c1*u + q = 0. The positive root (u = 1 at equilibrium) is
    % taken in cancellation-free form. dtheta_dK is d(theta_w)/dK.
    function [theta_w, u, dtheta_dK] = f_wall_theta_exact(Tm_in,theta_in)
        c1 = grad_Tm_coeff(1);
        c2 = grad_Trans_coeff(1);
        K  = grad_Tm_coeff(2:3)*Tm_in(2:3) + ...
            grad_Trans_coeff(2:3)*theta_in(end-1:-1:end-2);
        q  = 2*c1*(alpha_g - 1) - c2 + 2*alpha_g*K;
        u  = -q/(c1 + sqrt(c1^2 - c2*q));
        theta_w = (u^2 - 1)/(2*alpha_g);
        dtheta_dK = -u/(c2*u + c1);
    end
    
    % warm-started Newton iteration on the full (thermal and mass transfer)
    % wall balance, with fzero as the fallback if it does not converge
    function theta_bw = f_wall_theta_newton(theta_bw)
        converged = false;
        for iter = 1:10
            [res,dres] = f_bubble_wall_full_bc(theta_bw);
            step = res/dres;
            theta_bw = theta_bw - step;
            if ~isfinite(theta_bw)
                break;
            end
            if abs(step) <= 1e-14*max(1,abs(theta_bw))
                converged = true;
                break;
            end
        end
        if ~converged
            theta_bw = fzero(@f_bubble_wall_full_bc,theta_bw_guess,foptions);
        end
    end
    
    % finite difference matrices, second-order central in the interior and
    % second-order one-sided at the bubble wall
    % nodes: Number of nodes
    % order: order of differentiation ( 1st derivative vs 2nd derivative)
    % Tm_check: 0 not for external temp , 1 used for external temp
    function [dmatrix] = f_finite_diff_mat(nodes,order,Tm_check)
        
        % coordinate creation
        N = nodes-1;
        if Tm_check == 0
            deltaY = 1/N;
        else
            deltaY = -2/N;
        end
        K = (1:N+1)';
        y = (K-1)*deltaY;
        % interior rows 2..N; spdiags stores diagonal d of an n-by-n matrix
        % in column j of the band, so A(i,i+d) = band(i+d)
        interior = (K >= 2) & (K <= N);
        sub = zeros(nodes,1);
        dia = zeros(nodes,1);
        sup = zeros(nodes,1);
        
        if order == 1
            sup([false; interior(1:end-1)]) = 0.5;
            sub([interior(2:end); false]) = -0.5;
            dmatrix = spdiags([sub dia sup], -1:1, nodes, nodes);
            if Tm_check == 0
                dmatrix(end,end-2:end) = [0.5 -2 1.5];
            end
            dmatrix = dmatrix / deltaY;
            
        elseif order == 2
            if Tm_check == 0
                % rows i = 2..N: 1 - deltaY/y(i), -2, 1 + deltaY/y(i)
                ii = find(interior);
                sup(ii+1) = 1 + deltaY./y(ii);
                dia(ii) = -2;
                sub(ii-1) = 1 - deltaY./y(ii);
                % symmetry condition at the bubble center
                dia(1) = -6;
                sup(2) = 6;
            else
                sup([false; interior(1:end-1)]) = 1;
                dia(interior) = -2;
                sub([interior(2:end); false]) = 1;
            end
            dmatrix = spdiags([sub dia sup], -1:1, nodes, nodes);
            dmatrix = dmatrix / (deltaY^2);
        end
        
    end
    
end
