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
    Pv_star         = init_opts(4);
    Req             = init_opts(5);
    alphax          = init_opts(6);
    
    % dimensionaless initial stress
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
    thetav          = mass_opts(8);
    
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
    % dimensionless saturated vapor pressure, Pv(T) = C1_pv*exp(C2_pv/T)
    % with Pv(1) = Pv_star
    C1_pv = Pv_star*exp(thetav);
    C2_pv = -thetav;
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
    % constants of the bubble-wall heat (and mass) flux balance
    wallc = struct('grad_Tm',grad_Tm_coeff, 'grad_Tr',grad_Trans_coeff, ...
        'grad_C',grad_C_coeff, 'alpha_g',alpha_g, 'alpha_v',alpha_v, ...
        'C1',C1_pv, 'C2',C2_pv, 'Rvg',Rvg_ratio, 'Rva',Rva_diff, ...
        'Rg',Rg_star);
    
    % solver start
    f_display(radial, bubtherm, medtherm, masstrans, stress, spectral,...
        nu_model, eps3, Pv_star, Re8, De, Ca, LAM, 'finite difference');
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
        Sp = @(t,X) f_jacobian_fd(t,X,jp);
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
            % warm-started Newton on the transcendental wall balance
            theta_bw = f_wall_theta_newton([Tm(2), Tm(3), theta(end-1), ...
                theta(end-2), kv(end), kv(end-1), kv(end-2), P], theta_bw_guess);
            theta_bw_guess = theta_bw;
            theta(end) = theta_bw;
        elseif medtherm
            % temperature in the material
            Tm = X(imedtherm);
            % closed-form wall temperature, no root finding
            theta_bw = f_wall_thermal_bc(grad_Tm_coeff(2:3)*Tm(2:3) + ...
                grad_Trans_coeff(2:3)*theta(end-1:-1:end-2), ...
                grad_Tm_coeff(1), grad_Trans_coeff(1), alpha_g);
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
        [Pf8_loc,Pf8dot_loc] = f_pinfinity(t,pvarargin);
        [S_loc,Sdot_loc,Z1dot_loc,Z2dot_loc] = f_stress(stress,X,Req,Rv, ...
            Ca,De,Re8,Rdotv,alphax,ivisco1,ivisco2,LAM,zeNO,cdd, ...
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
    
    % warm-started Newton iteration on the full (thermal and mass transfer)
    % wall balance, with fzero as the fallback if it does not converge;
    % vrest = [Tm(2), Tm(3), theta(end-1), theta(end-2), kv(end),
    % kv(end-1), kv(end-2), P]
    function theta_bw = f_wall_theta_newton(vrest, theta_bw)
        theta_start = theta_bw;
        converged = false;
        for iter = 1:10
            [res,dres] = f_wall_mass_bc([theta_bw, vrest], wallc);
            step = res/dres(1);
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
            theta_bw = fzero(@(x) f_wall_mass_bc([x, vrest], wallc), ...
                theta_start, foptions);
        end
    end
    
    % bubble-wall theta at the state X (mass transfer), for the Jacobian
    function theta_bw = f_theta_bw_of_X(X)
        th_X = X(ibubtherm);
        Tm_X = X(imedtherm);
        kv_X = X(imass);
        theta_bw = f_wall_theta_newton([Tm_X(2), Tm_X(3), th_X(end-1), ...
            th_X(end-2), kv_X(end), kv_X(end-1), kv_X(end-2), X(3)], ...
            theta_bw_guess);
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
            sup(vertcat(false, interior(1:end-1))) = 0.5;
            sub(vertcat(interior(2:end), false)) = -0.5;
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
                sup(vertcat(false, interior(1:end-1))) = 1;
                dia(interior) = -2;
                sub(vertcat(interior(2:end), false)) = 1;
            end
            dmatrix = spdiags([sub dia sup], -1:1, nodes, nodes);
            dmatrix = dmatrix / (deltaY^2);
        end
        
    end
    
end
