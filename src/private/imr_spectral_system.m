% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Scheme-specific operators, closures and reconstruction for one checked case.
% Public calls and time integration remain in f_imr_spectral.
function system = imr_spectral_system(prepared)
    inputs = prepared.inputs;
    layout = prepared.layout;
    % defaults

    % non-Newtonian viscosity parameters
    fnu_m           = 0;
    intfnu          = 0;
    dintfnu         = 0;
    ddintfnu        = 0;

    % equations settings
    radial          = prepared.model.radial;
    bubtherm        = prepared.model.bubtherm;
    medtherm        = prepared.model.medtherm;
    stress          = prepared.model.stress;
    eps3            = prepared.model.eps3;

    % solver options
    spectral        = prepared.numerics.spectral;
    Nv              = prepared.numerics.Nv;
    Nt              = prepared.numerics.Nt;
    Mt              = prepared.numerics.Mt;
    Lv              = prepared.numerics.Lv;
    Lt              = prepared.numerics.Lt;
    jacobian        = prepared.numerics.jacobian;

    % dimensionless initial conditions
    Rzero           = prepared.initial.R;
    Uzero           = prepared.initial.Rdot;
    Pb_star         = prepared.initial.P;
    Pv_star         = prepared.initial.Pv;
    Req             = prepared.initial.Req;
    alphax          = prepared.stress.alphax;

    % dimensionless initial stress
    Szero           = prepared.initial.stress;

    % time span options
    tspan = prepared.time;

    % output options
    progdisplay     = prepared.output.progress;

    % physical parameters

    % acoustic parameters
    Cstar           = prepared.acoustic.Cstar;
    GAMa            = prepared.acoustic.GAMa;
    kappa           = prepared.acoustic.kappa;
    nstate          = prepared.acoustic.nstate;
    hugoniot_s      = prepared.acoustic.hugoniot_s;

    % dimensionless waveform parameters
    om              = prepared.forcing.om;
    ee              = prepared.forcing.ee;
    tw              = prepared.forcing.tw;
    dt              = prepared.forcing.dt;
    mn              = prepared.forcing.mn;
    wave_type       = prepared.forcing.wave_type;
    wave_poly       = prepared.forcing.wave_poly;
    wave_dpoly      = prepared.forcing.wave_dpoly;

    pvarargin = {[om,ee,tw,dt,mn,wave_type], wave_poly, wave_dpoly};
    % the impulse waveform is the constant pressure ee (f_pinfinity returns
    % exactly ee and 0), so the right-hand side skips the call
    wave_const = (wave_type == 0);

    % dimensionless viscoelastic
    We              = prepared.stress.We;
    Re8             = prepared.stress.Re8;
    DRe             = prepared.stress.DRe;
    v_a             = prepared.stress.v_a;
    v_nc            = prepared.stress.v_nc;
    Ca              = prepared.stress.Ca;
    LAM             = prepared.stress.LAM;
    De              = prepared.stress.De;
    JdotA           = prepared.stress.JdotA;
    nu_model        = prepared.stress.nu_model;
    v_lambda_star   = prepared.stress.v_lambda_star;
    iDRe            = prepared.stress.iDRe;
    iWe             = 1/We;

    % dimensionless thermal
    Foh             = prepared.thermal.Foh;
    Br              = prepared.thermal.Br;
    alpha_g         = prepared.thermal.alpha_g;
    beta_g          = prepared.thermal.beta_g;
    % Vapor conductivity coefficients are used by the nodal solver only.
    chi             = prepared.thermal.chi;
    iota            = prepared.thermal.iota;

    % pre_process code

    % Enabled fields own their operators; disabled fields remain empty.
    y = [];
    xi = [];
    gA = [];
    gAI = [];
    gDA = [];
    gDDA = [];
    mA = [];
    mDA = [];
    mDDA = [];
    if bubtherm
        y = cos(pi*(0:Nt)'/(2*Nt));
        [gA,gAI,~,~,gAPd,gAPdd] = f_dcdmtxe(Nt);
        gDA = gAPd*gA;
        gDDA = gAPdd*gA;
    end
    if medtherm
        xi = cos(pi*(0:Mt)'/Mt);
        [mA,~,~,~,mAPd,mAPdd] = f_dcdmtx(Mt);
        mDA = mAPd*mA;
        mDDA = mAPdd*mA;
    end
    % thermal boundary condition system: only its first row varies in time,
    % linearly through s = alpha_g*(T(1)-1)+1, so factor the reference
    % matrix at s = 1 once and apply a Sherman-Morrison rank-one update
    if bubtherm && medtherm
        Q = [gA(2:end,:) zeros(Nt,Mt+1);
        zeros(Mt,Nt+1) mA(2:end,:);
        2*(0:Nt).^2 iota*(0:Mt).^2];
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
    nn = [];
    if bubtherm
        nn = sparse(((-1).^(0:Nt).*(0:Nt).^2)');
    end

    % precomputations for viscous dissipation: stress basis at the
    % medium collocation points
    ZZT = [];
    cdd = [];
    if spectral_stress || zener_field
        cdd = f_preStressInt(Lv,Nv);
        if medtherm
            zT = 1 - 2./(1 + (yT - 1)/Lv);
            ZZT = cos(acos(zT)*(1:Nv)) - 1;
        end
    end

    % One state contract supplies every consumer, without sentinel orders.
    ia = layout.bubble;
    ib = layout.medium;
    ivisco1 = layout.stress1;
    ivisco2 = layout.stress2;
    izf = layout.polymer;
    iz = layout.stress;
    Nf = numel(izf);
    init = zeros(layout.count,1);
    init([layout.R layout.Rdot layout.P]) = [Rzero Uzero Pb_star]';
    if medtherm
        init(ib(1)) = 1;
    end
    init(iz) = Szero;

    % solver start
    f_display(radial, bubtherm, medtherm, 0, stress, spectral,...
        nu_model, eps3, Pv_star, Re8, De, Ca, LAM, 'spectral', DRe, iDRe, ...
        v_lambda_star, v_nc, v_a);
    bubble = @f_rhs;
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
    system.initial = init;
    system.rhs = bubble;
    system.jacobian = Sp;
    system.outputs = @f_outputs;
    system.operatorSizes = [size(gA,1) size(mA,1)];

    % Reconstruct boundary values and preserve the public output tuple.
    function out = f_outputs(t,X)

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
            T = f_thermal_inverse(gA*a,alpha_g,beta_g);
            if medtherm
                Tm = mA*b;
            end
        end

        % outputs
        out{1} = t;
        out{2} = R;
        out{3} = Rdot;
        out{4} = P;
        if bubtherm == 1
            out{5} = T;
        else
            out{5} = [];
        end
        if medtherm == 1
            out{6} = Tm;
        else
            out{6} = ones(size(t,1),1);
        end
        out{7} = [];
        % the inputs of this run (f_imr_inputs names)
        out{8} = inputs;

    end

    % solver function
    function dXdt = f_rhs(t,X)

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
            T = f_thermal_inverse(theta,alpha_g,beta_g);
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
            [P, Pdot] = f_polytropic_pressure(R, Rdot, ...
                Pb_star-Pv_star, Pv_star, kappa);
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
    % v = [R, Rdot, P, Pdot, Z], mirroring f_rhs
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
