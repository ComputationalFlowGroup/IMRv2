% file f_imr_spectral.m
% brief contains function f_imr_spectral

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function features a Chebyshev spectral collocation solver of
% the PDEs involving thermal transport and viscoelasticity coupled to the
% radial (Rayleigh-Plesset family) dynamics
function varargout = f_imr_spectral(varargin)
    
    % problem initialization
    [eqns_opts, solve_opts, init_opts, init_stress, tspan_opts, out_opts, ...
        acos_opts, wave_opts, sigma_opts, thermal_opts, ~, inputs] ...
        = f_call_params('spectral', varargin{:});
    
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
    
    % solver options
    method          = solve_opts(1);
    spectral        = solve_opts(2);
    divisions       = solve_opts(3);
    Nv              = solve_opts(4);
    Nt              = solve_opts(5);
    Mt              = solve_opts(6);
    Lv              = solve_opts(7);
    Lt              = solve_opts(8);
    rtol            = solve_opts(9);
    atol            = solve_opts(10);
    jacobian        = solve_opts(11);
    
    % dimensionless initial conditions
    Rzero           = init_opts(1);
    Uzero           = init_opts(2);
    Pb_star         = init_opts(3);
    Pv_star         = init_opts(4);
    Req             = init_opts(5);
    alphax          = init_opts(6);
    
    % dimensionless initial stress
    Szero           = init_stress;
    
    % time span options
    tspan = tspan_opts;
    
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
    % alpha_v, beta_v = thermal_opts(5:6), only used by the nodal solver
    chi             = thermal_opts(7);
    iota            = thermal_opts(8);
    
    % pre_process code
    
    % collocation point construction
    y = cos(pi*(0:Nt)'/(2*Nt));
    xi = cos(pi*(0:Mt)'/Mt);
    % collocation matrix construction
    [gA,gAI,~,~,gAPd,gAPdd] = f_dcdmtxe(Nt);
    [mA,~,~,~,mAPd,mAPdd] = f_dcdmtx(Mt);
    % combined operators from Chebyshev coefficients to the collocation
    % values of the first and second derivatives
    gDA = gAPd*gA;
    gDDA = gAPdd*gA;
    mDA = mAPd*mA;
    mDDA = mAPdd*mA;
    Q = [gA(2:end,:) zeros(Nt,Mt+1);
    zeros(Mt,Nt+1) mA(2:end,:);
    2*(0:Nt).^2 iota*(0:Mt).^2];
    Q = sparse(Q);
    % thermal boundary condition system: only its first row varies in time,
    % linearly through s = alpha_g*(T(1)-1)+1, so factor the reference
    % matrix at s = 1 once and apply a Sherman-Morrison rank-one update
    if bubtherm && medtherm
        A0 = vertcat([ones(1,Nt+1), -ones(1,Mt+1)], Q);
        [LA0,UA0,PA0,QA0] = lu(sparse(A0));
        e1 = vertcat(1, zeros(Nt+Mt+1,1));
        vrow = [zeros(1,Nt+1), ones(1,Mt+1)]';
        w1 = QA0*(UA0\(LA0\(PA0*e1)));
        vw1 = vrow'*w1;
    end
    % spectral stress representation (UCM, Oldroyd-B): tau_rr and
    % tau_thetatheta are expanded in T_n(z) - 1, n = 1..Nv, which vanish at
    % z = 1 (r -> Inf), with z = 1 - 2/(1 + (r/R - 1)/Lv) and collocation
    % points z_j = cos(pi*j/Nv), j = 1..Nv (z = -1 is the bubble wall)
    spectral_stress = (stress == 5 && spectral == 1);
    % the linear Maxwell, Jeffreys and Zener models (stress 3, 4) carry,
    % with medium heat transfer, the local polymer stress difference as a
    % field in the same basis (f_zener_field), for the medium heating
    zener_field = (stress == 3 || stress == 4) && medtherm == 1;
    if spectral_stress || zener_field
        ze = cos(pi*(1:Nv)'/Nv);
        % r/R at the stress collocation points
        yV = 2*Lv./(1-ze) - Lv + 1;
        % coefficients -> values, values -> coefficients, coefficients ->
        % dz derivative values
        [sCA,sCI,sCAd,~,~,~] = f_dcdmtx(Nv);
        sCA = sCA(2:end,2:end) - 1;
        sCI = sCI(2:end,2:end);
        sCAd = sCAd(2:end,2:end);
        % dz/d(r/R) at the collocation points
        dzdy = (1-ze).^2/(2*Lv);
        % constant factors of the stress transport, formed with the same
        % operations as in f_spectral_stress, so results are unchanged
        yV3 = yV.^3;
        yV_adv = 1./yV.^2 - yV;
        % polymer viscosity and its relaxation coefficient
        etap = (1-LAM)/Re8;
    end
    
    % precomputations
    sam = 1 + GAMa;
    no = (nstate-1)/nstate;
    nog = (nstate-1)/2;
    kapover = (kappa-1)/kappa;
    yT = 2*Lt./(1+xi) - Lt + 1;
    yT2 = yT.^2;
    yT3 = yT.^3;
    iyT3 = yT.^-3;
    iyT4 = yT.^-4;
    iyT6 = yT.^-6;
    % constant factors of the medium temperature equation, formed with the
    % same operations as in the right-hand side, so results are unchanged
    xi1_sq = (1+xi).^2;
    xi_adv = (1+xi)/(2*Lt) - 1./yT;
    yT_adv = 1./yT2 - yT;
    diff_m = Foh/4*(1+xi).^4;
    % Kirchhoff transform of K = alpha_g*T + beta_g,
    % theta = alpha_g*(T^2 - 1)/2 + beta_g*(T - 1), inverted as
    % T = (sqrt((alpha_g + beta_g)^2 + 2*alpha_g*theta) - beta_g)/alpha_g
    K1sq_g = (alpha_g + beta_g)^2;
    nn = ((-1).^(0:Nt).*(0:Nt).^2)';
    nn = sparse(nn);
    
    % precomputations for viscous dissipation: stress basis at the
    % medium collocation points
    zT = 1 - 2./(1 + (yT - 1)/Lv);
    ZZT = cos(acos(zT)*(1:Nv)) - 1;
    % stress integral weights, S = 2*int_R^Inf (tau_rr - tau_tt)/r dr =
    % 2*cdd'*(c - d) for the spectral stress coefficients c, d
    cdd = f_preStressInt(Lv,Nv);
    
    % index management: Nf spectral modes of the Zener field, and one
    % state per stress integral unless the stress 5 field is used
    Nf = Nv*zener_field;
    if spectral == 0
        Nv = 1;
    end
    if bubtherm == 0
        Nt = -1;
        Mt = -1;
    end
    if medtherm == 0
        Mt = -1;
    end
    ia = 4:(4+Nt);
    ib = (5+Nt):(5+Nt+Mt);
    if zener_field
        % Z1, then the Nf coefficients of the Zener field
        ivisco1 = 6+Nt+Mt;
        ivisco2 = [];
        izf = ivisco1 + (1:Nf);
    else
        ivisco1 = (6+Nt+Mt):(5+Nt+Mt+Nv);
        ivisco2 = (6+Nt+Mt+Nv):(5+Nt+Mt+2*Nv);
        izf = [];
    end
    
    % initial condition assembly
    
    % radius, velocity, pressure
    
    % auxiliary temperature, boundary temperature, medium temperature,
    Tau0 = zeros(Nt+1,1);
    Tm0 = ones(Mt ~= -1);
    Tm1 = zeros(Mt,1);
    
    % initial condition vector
    init = [Rzero;
    Uzero;
    Pb_star;
    Tau0;
    Tm0;
    Tm1;
    Szero;
    zeros(Nf,1)];
    % the stress states present in the state vector (Z1 and Z2, or the
    % stress 5 field coefficients), without the Zener field
    iz = (6+Nt+Mt):(numel(init)-Nf);
    
    % solver start
    f_display(radial, bubtherm, medtherm, 0, stress, spectral,...
        nu_model, eps3, Pv_star, Re8, De, Ca, LAM, 'spectral', DRe, iDRe, ...
        v_lambda_star, v_nc, v_a);
    bubble = @SVBDODE;
    % analytic Jacobian for the stiff integrators
    if jacobian == 2
        jp = struct('ia',ia, 'ib',ib, 'iz',iz, 'izf',izf, ...
            'bubtherm',bubtherm, 'medtherm',medtherm, ...
            'spectral_stress',spectral_stress, 'gA',gA, 'gDA',gDA, ...
            'gDDA',gDDA, 'gAI',gAI, 'mDA',mDA, 'mDDA',mDDA, ...
            'y',y, 'xi',xi, 'yT',yT, 'yT2',yT2, 'yT3',yT3, 'Lt',Lt, ...
            'Foh',Foh, 'Br',Br, 'chi',chi, 'kappa',kappa, 'kapover',kapover, ...
            'alpha_g',alpha_g, 'beta_g',beta_g, 'nn',nn, 'Pb_star',Pb_star, ...
            'Pv_star',Pv_star, 'ZZT',ZZT, 'cdd',cdd, 'LAM',LAM, 'Re8',Re8, ...
            'closure',@f_closure, 'closure_s',@f_closure_s, ...
            'taugradu',@f_taugradu, 'zener_field',zener_field, ...
            'zener_rhs',@f_zener_rhs);
        if bubtherm && medtherm
            jp.LA0 = LA0;
            jp.UA0 = UA0;
            jp.PA0 = PA0;
            jp.QA0 = QA0;
            jp.w1 = w1;
            jp.vw1 = vw1;
            jp.vrow = vrow;
        end
        if spectral_stress || zener_field
            jp.sCA = sCA;
            jp.sCI = sCI;
            jp.sCAd = sCAd;
            jp.dzdy = dzdy;
            jp.yV = yV;
            jp.De = De;
            jp.etap = etap;
        end
        Sp = @(t,X) f_jacobian_spectral(t,X,jp);
    else
        Sp = [];
    end
    [t,X] = f_odesolve(bubble, init, method, divisions, tspan, Sp, rtol, atol);
    
    % post processing
    
    % debug check of the Zener field: its stress integral 2*cdd'*e must
    % equal Z1/R^3, the state that drives R(t), up to the resolution of the
    % Nf modes
    if zener_field
        Sfield = 2*X(:,izf)*cdd;
        Sstate = X(:,ivisco1)./X(:,1).^3;
        fprintf(['Zener field check: max |2*cdd''*e - Z1/R^3| = %.3e ' ...
            '(max |Z1/R^3| = %.3e, %d modes)\n'], max(abs(Sfield - Sstate)), ...
            max(abs(Sstate)), Nf);
    end
    
    % extract result
    R = X(:,1);
    Rdot = X(:,2);
    P = X(:,3);
    % extracting the Chebyshev coefficients
    a = X(:,ia)';
    b = X(:,ib)';
    if bubtherm
        T = (sqrt(K1sq_g + 2*alpha_g*gA*a) - beta_g)/alpha_g;
        if medtherm
            Tm = mA*b;
        end
    end
    
    % outputs
    varargout{1} = t;
    varargout{2} = R;
    varargout{3} = Rdot;
    varargout{4} = P;
    if bubtherm == 1
        varargout{5} = T;
    else
        varargout{5} = [];
    end
    if medtherm == 1
        varargout{6} = Tm;
    else
        varargout{6} = ones(size(t,1),1);
    end
    varargout{7} = [];
    % the inputs of this run (f_imr_inputs names)
    varargout{8} = inputs;
    
    % solver function
    function dXdt = SVBDODE(t,X)
        
        % showing output
        if progdisplay
            disp(t/tspan(end));
        end
        
        % extract standard inputs
        
        % bubble wall radius
        R = X(1);
        % bubble wall velocity
        Rdot = X(2);
        % internal bubble pressure
        P = X(3);
        % thermal vector in and outside the bubble
        qdot = [];
        
        % updating the viscous forces/Reynolds number
        if nu_model ~= 0
            % fnu_m: viscosity factor at the medium grid points
            [~,intfnu,dintfnu,ddintfnu,fnu_m] = f_viscosity(nu_model, ...
                Rdot,R,v_a,v_nc,v_lambda_star,yT);
        end
        
        % bubble equations of motion
        if bubtherm
            % auxiliary temperature and its first and second derivatives
            % from the Chebyshev coefficients
            theta = gA*X(ia);
            dtheta = gDA*X(ia);
            ddtheta = gDDA*X(ia);
            % temperature and thermal diffusivity fields
            T = (sqrt(K1sq_g + 2*alpha_g*theta) - beta_g)/alpha_g;
            D = kapover*(alpha_g*T.^2 + beta_g*T)/P;
            
            % bubble pressure equation
            Pdot = 3/R*((kappa-1)*chi/R*dtheta(1) - kappa*P*Rdot);
            
            % bubble wall velocity equation
            U_vel = chi/R^2*(2*D./y - kapover/P*(dtheta-dtheta(1)*y));
            
            % auxiliary temperature derivative
            thetadot = Pdot*D + U_vel.*dtheta + chi*D/R^2.*ddtheta;
            
            thetadot(end) = Pdot*D(end) - chi/R^2*(8*D(end)*sum(nn.*X(ia)) ...
                + kapover/P*dtheta(end)^2) + chi*D(end)/R^2.*ddtheta(end);
            
            % surroundings equations of motion
            
            % surrounding temperature
            if medtherm
                % material temperature equations, derivatives from the
                % Chebyshev coefficients
                advection = xi1_sq/(Lt*R).*...
                    (Foh/R*xi_adv + Rdot/2*yT_adv).*(mDA*X(ib));
                diffusion = diff_m/(Lt^2*R^2).*(mDDA*X(ib));
                % stress dissipation
                [taugradu] = f_stress_dissipation(stress,Req,R,Rdot,Ca, ...
                    Br,Re8,LAM,alphax,iDRe,yT2,yT3,iyT3,iyT4,iyT6,fnu_m, ...
                    f_polymer_field(X));
                % surrounding temperature evolution equation
                Tmdot = advection + diffusion + taugradu;
                
                % enforce boundary condition and solve A*x=b problem with
                % the Sherman-Morrison update of the prefactored system
                Tmdot(end) = 0;
                s = alpha_g*T(1) + beta_g;
                rhs = vertcat(0, thetadot(2:end), Tmdot(2:end), 0);
                y0 = QA0*(UA0\(LA0\(PA0*rhs)));
                qdot = y0 + w1*(((s-1)*(vrow'*y0)) / (1 - (s-1)*vw1));
            else
                % cold-liquid approximation
                % solve auxiliary temperature with boundary condition
                qdot = gAI*[0;
                thetadot(2:end)];
            end
        else
            % polytropic gas + vapor
            P = (Pb_star-Pv_star)*(1/R)^(3*kappa) + Pv_star;
            Pdot = -3*kappa*P*Rdot/R;
        end
        
        % pressure waveform
        if wave_const
            Pf8 = ee;
            Pf8dot = 0;
        else
            [Pf8,Pf8dot] = f_pinfinity(t,pvarargin);
        end
        
        % stress equation evolution
        if spectral_stress
            [S,Sdot,Z1dot,Z2dot] = f_spectral_stress(R,Rdot,X);
        else
            [S,Sdot,Z1dot,Z2dot] = f_stress(stress,X,Req,R,Ca,De,Re8, ...
                Rdot,alphax,ivisco1,ivisco2,LAM,intfnu,dintfnu,iDRe);
        end
        
        % bubble wall evolution / acceleration
        [Rddot] = f_radial_eq(radial, P, Pdot, Pf8, Pf8dot, iWe, R, Rdot, S, ...
            Sdot, Cstar, sam, no, GAMa, nstate, nog, hugoniot_s, JdotA, ...
            ddintfnu, iDRe);
        
        % local polymer stress of the Zener models, for the heating
        if zener_field
            edot = f_zener_rhs(R,Rdot,X);
        else
            edot = [];
        end
        
        % output assembly
        dXdt = [Rdot;
        Rddot;
        Pdot;
        qdot;
        Z1dot;
        Z2dot;
        edot];
        
    end
    % end of solver
    
    % local polymer stress difference at the medium points, for the medium
    % heating of stress 3-5 (empty otherwise): the stress 5 field
    % tau_rr - tau_tt, or the Zener field
    function sp = f_polymer_field(X)
        if spectral_stress
            sp = ZZT*(X(ivisco1) - X(ivisco2));
        elseif zener_field
            sp = ZZT*X(izf);
        else
            sp = [];
        end
    end
    
    % time derivative of the Chebyshev coefficients e of the Zener field
    % sp = sCA*e at the stress collocation points (f_zener_field), with
    % the advection -ydot*d(sp)/dy = -dz/dt*d(sp)/dz
    function edot = f_zener_rhs(R,Rdot,X)
        e = X(izf);
        dzdt = dzdy.*(Rdot/R).*yV_adv;
        spdot = f_zener_field(sCA*e, -dzdt.*(sCAd*e), yV, R, Rdot, Req, Ca, ...
            De, Re8, LAM, alphax*(stress == 4));
        edot = sCI*spdot;
    end
    
    % UCM and Oldroyd-B stress fields in the spectral representation. The
    % polymer stresses obey, along material points with a = R^2*Rdot/r^3,
    %   D(tau_rr)/Dt = -tau_rr/De - 4*a*tau_rr - 4*etap*a/De
    %   D(tau_tt)/Dt = -tau_tt/De + 2*a*tau_tt + 2*etap*a/De
    % with etap = (1-LAM)/Re8; in the moving coordinate z the material
    % derivative adds the advection dz/dt*d(tau)/dz, where
    % dz/dt = dz/d(r/R)*(Rdot/R)*((r/R)^-2 - r/R). The solvent viscosity
    % LAM/Re8 enters the stress integral directly.
    function [S,Sdot,cdot,ddot] = f_spectral_stress(R,Rdot,X)
        c = X(ivisco1);
        d = X(ivisco2);
        trr = sCA*c;
        ttt = sCA*d;
        av = (Rdot/R)./yV3;
        dzdt = dzdy.*(Rdot/R).*yV_adv;
        trrdot = -dzdt.*(sCAd*c) - trr/De - 4*av.*trr - 4*etap*av/De;
        tttdot = -dzdt.*(sCAd*d) - ttt/De + 2*av.*ttt + 2*etap*av/De;
        cdot = sCI*trrdot;
        ddot = sCI*tttdot;
        S = 2*cdd'*(c - d) - 4*LAM/Re8*Rdot/R - 6*intfnu*iDRe;
        % the Rddot parts of the solvent and non-Newtonian terms are carried
        % by JdotA and ddintfnu in the radial equation
        Sdot = 2*cdd'*(cdot - ddot) + 4*LAM/Re8*(Rdot/R)^2 - 6*dintfnu*iDRe;
    end
    
    % medium dissipation as a function of (R, Rdot), at the polymer stress
    % field of X
    function tg = f_taugradu(Rv,Rdotv,X)
        fnu_loc = 0;
        if nu_model ~= 0
            [~,~,~,~,fnu_loc] = f_viscosity(nu_model,Rdotv,Rv,v_a,v_nc, ...
                v_lambda_star,yT);
        end
        tg = f_stress_dissipation(stress,Req,Rv,Rdotv,Ca,Br,Re8,LAM, ...
            alphax,iDRe,yT2,yT3,iyT3,iyT4,iyT6,fnu_loc,f_polymer_field(X));
        tg = tg.*ones(size(yT));
    end
    
    % radial acceleration and stress evolution as a function of
    % v = [R, Rdot, P, Pdot, Z], mirroring SVBDODE
    function out = f_closure(t,v,X)
        X(iz) = v(5:end);
        [intf_loc,dintf_loc,ddintf_loc] = f_nu(v(1),v(2));
        [S_loc,Sdot_loc,Z1dot_loc,Z2dot_loc] = f_stress(stress,X,Req,v(1), ...
            Ca,De,Re8,v(2),alphax,ivisco1,ivisco2,LAM, ...
            intf_loc,dintf_loc,iDRe);
        out = vertcat(f_rddot(t,v,S_loc,Sdot_loc,ddintf_loc), Z1dot_loc, Z2dot_loc);
    end
    
    % radial acceleration as a function of v = [R, Rdot, P, Pdot, S, Sdot]
    % for the spectral stress representation, where S and Sdot are the
    % field and solvent parts; the non-Newtonian part is added here from
    % (R, Rdot)
    function out = f_closure_s(t,v)
        [intf_loc,dintf_loc,ddintf_loc] = f_nu(v(1),v(2));
        out = f_rddot(t,v,v(5) - 6*intf_loc*iDRe, v(6) - 6*dintf_loc*iDRe, ...
            ddintf_loc);
    end
    
    % non-Newtonian stress integral terms at (R, Rdot), zero when Newtonian
    function [intf_loc,dintf_loc,ddintf_loc] = f_nu(Rv,Rdotv)
        intf_loc = 0;
        dintf_loc = 0;
        ddintf_loc = 0;
        if nu_model ~= 0
            [~,intf_loc,dintf_loc,ddintf_loc] = f_viscosity(nu_model, ...
                Rdotv,Rv,v_a,v_nc,v_lambda_star);
        end
    end
    
    function Rddot_loc = f_rddot(t,v,S_loc,Sdot_loc,ddintf_loc)
        if wave_const
            Pf8_loc = ee;
            Pf8dot_loc = 0;
        else
            [Pf8_loc,Pf8dot_loc] = f_pinfinity(t,pvarargin);
        end
        Rddot_loc = f_radial_eq(radial, v(3), v(4), Pf8_loc, Pf8dot_loc, ...
            iWe, v(1), v(2), S_loc, Sdot_loc, Cstar, sam, no, GAMa, nstate, ...
            nog, hugoniot_s, JdotA, ddintf_loc, iDRe);
    end
    
end

% precomputation functions

% discrete Chebyshev derivative matrices, FCD
function [A,B,D,E,C,F] = f_dcdmtx(N)
    
    jj = 0:N;
    jj2 = jj.^2;
    theta = pi*jj'/N;
    
    % matrix for a -> p
    A = cos(theta*jj);
    
    % matrix for p -> a
    B = 2*cos(theta*jj)'/N;
    B(:,[1 N+1]) = 0.5*B(:,[1 N+1]);
    B([1 N+1],:) = 0.5*B([1 N+1],:);
    
    theta = theta(2:N);
    
    % matrix for a -> dp/dx
    D = zeros(N+1);
    D(1,:) = jj2;
    D(2:N,:) = csc(theta)*jj.*sin(theta*jj);
    D(N+1,:) = jj2.*((-1).^(jj+1));
    
    % matrix for a -> d2p/dx2
    E = zeros(N+1);
    E(1,:) = jj2.*(jj2-1)/3;
    E(2:N,:) = (cos(theta)./sin(theta).^3)*jj.*sin(theta*jj) ...
        - (csc(theta).^2)*((0:N).^2).*cos(theta*jj);
    E(N+1,:) = jj2.*(jj2-1).*(-1).^jj/3;
    
    % matrix for p -> dp/dx
    C = D*B;
    
    % matrix for p -> d2p/dx2
    F = E*B;
    
end

% even discrete Chebyshev derivative matrices, FCD
function [A,B,D,E,C,F] = f_dcdmtxe(N)
    
    jj = 2*(0:N);
    jj2 = jj.^2;
    theta = pi*jj'/(4*N);
    
    % matrix for a -> p
    A = cos(theta*jj);
    
    % matrix for p -> a
    B = 2*cos(theta*jj)'/N;
    B(:,[1 N+1]) = 0.5*B(:,[1 N+1]);
    B([1 N+1],:) = 0.5*B([1 N+1],:);
    
    theta = theta(2:N+1);
    
    % matrix for a -> dp/dx
    D = zeros(N+1);
    D(1,:) = jj2;
    D(2:N+1,:) = csc(theta)*jj.*sin(theta*jj);
    
    % matrix for a -> d2p/dx2
    E = zeros(N+1);
    E(1,:) = jj2.*(jj2-1)/3;
    E(2:N+1,:) = (cos(theta)./sin(theta).^3)*jj.*sin(theta*jj) ...
        - (csc(theta).^2)*(jj2).*cos(theta*jj);
    
    % matrix for p -> dp/dx
    C = D*B;
    
    % matrix for p -> d2p/dx2
    F = E*B;
    
end

% prestress integrals, fixed Gauss-Legendre quadrature in theta; the (cos-1)
% zero at x = 1 tames the (1-x)^-2 endpoint
function cdd = f_preStressInt(L,N)
    
    % with x = cos(theta) the integrands of
    % 2*L*int_{-1}^{1} (T_n(x) - 1)/((L*(1+x) + 1 - x)*(1 - x)) dx become
    % smooth, (cos(n*theta) - 1)*sin(theta)/(((L+1) + (L-1)*cos(theta))*
    % (1 - cos(theta))), so one Gauss-Legendre rule on [0, pi] with N + 64
    % points gives every n to round-off
    M = N + 64;
    k = (1:M-1)';
    beta = k./sqrt(4*k.^2 - 1);
    [V,D] = eig(diag(beta,1) + diag(beta,-1));
    [x,i] = sort(diag(D));
    th = pi/2*(x + 1);
    w = pi*V(1,i)'.^2;
    w = w.*sin(th)./(((L + 1) + (L - 1)*cos(th)).*(1 - cos(th)));
    cdd = 2*L*((cos((1:N)'*th') - 1)*w);
    
end
