% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Scheme-specific operators, closures and reconstruction for one checked case.
% Public calls and time integration share imr.solvers.run.
%
% BUILDER STRUCTURE (roughly in order):
%  1. Parameter unpacking: extract model equations (radial, bubtherm, etc.),
%     solver options (Nv/Nt/Mt/Lv/Lt: spectral grid sizes), and the initial
%     state (Rzero, Pb_star, Szero).
%  2. Spectral basis setup: Chebyshev collocation grids (yV, zV), forward/
%     inverse projection maps (sCA/sCI for stress field if present, gA/gAI
%     for thermal).
%  3. Operator setup: derivative matrices in spectral space, mapped to
%     physical grids.
%  4. Physical constants: conductivity coefficients, saturation pressure,
%     gas constants, non-dimensional parameters.
%  5. Thermal operators and coefficients (Kirchhoff transform, Sherman-
%     Morrison terms for the thermal boundary condition).
%  6. State layout and indices: maps state to R, Rdot, P, thermal/medium
%     coefficients, stress modes, polymer field modes.
%  7. Initial condition assembly: populate init from prepared values.
%  8. Jacobian setup: conditionally build jp struct (parameters for analytic
%     Jacobian) or set up finite-difference Jacobian.
%  9. System assembly: package system.initial, system.rhs, system.jacobian,
%     system.outputs.
% 10. Nested functions: f_outputs, f_rhs and all RHS helpers (f_zener_rhs,
%     f_spectral_stress, f_closure, etc.).
%
% Spectral representation: fields are represented as Chebyshev mode
% coefficients. For example, bubble temperature theta(y,t) = gA*a(t) where
% a are the coefficients and gA is the forward basis matrix. The inverse
% map gAI projects node values back to coefficients (= gA^{-1}).
function system = build_spectral(prepared)
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
    C2              = prepared.stress.C2;

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
    % the impulse waveform is the constant pressure ee (imr.physics.far_field_pressure returns
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
    nonnewtonian    = nu_model ~= 0 && iDRe ~= 0;

    % dimensionless thermal
    Foh             = prepared.thermal.Foh;
    Br              = prepared.thermal.Br;
    alpha_g         = prepared.thermal.alpha_g;
    beta_g          = prepared.thermal.beta_g;
    alpha_v         = prepared.thermal.alpha_v;
    beta_v          = prepared.thermal.beta_v;
    chi             = prepared.thermal.chi;
    iota            = prepared.thermal.iota;

    % dimensionless mass transfer: this solver has no evolving vapor field,
    % but with vapor present the bubble is a frozen gas-vapor mixture at
    % its initial fraction kv0, same effective conductivity imr.solve_nodal
    % uses without mass transfer (kv0 = 0 without vapor, so alpha_m0 and
    % beta_m0 reduce exactly to alpha_g and beta_g)
    kv0             = prepared.mass.kv0;
    alpha_m0        = kv0*alpha_v + (1-kv0)*alpha_g;
    beta_m0         = kv0*beta_v + (1-kv0)*beta_g;

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
    % linearly through s = alpha_m0*(T(1)-1)+1, so factor the reference
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
    % spectral stress representation (UCM/Oldroyd-B, PTT, Giesekus): tau_rr
    % and tau_thetatheta are expanded in T_n(z) - 1, n = 1..Nv, which vanish
    % at z = 1 (r -> Inf), with z = 1 - 2/(1 + (r/R - 1)/Lv) and collocation
    % points z_j = cos(pi*j/Nv), j = 1..Nv (z = -1 is the bubble wall).
    % Stress 6 (PTT) and 7 (Giesekus) have no closed-form reduction and run
    % only here (imr.configuration.validate_case forces spectral = 1 for them)
    spectral_stress = (stress == 5 || stress == 6 || stress == 7) && spectral == 1;
    % PTT: stress == 6 sets the Warnez & Johnsen (2015) eps2 = 1 exponential
    % closure (f_spectral_stress). Giesekus (stress == 7) instead uses their
    % eps3 as the quadratic mobility factor; eps3 = 0 for stress 5 and 6.
    ptt = (stress == 6);
    % the linear, 2nd-order Fung and Mooney-Rivlin Zener models (stress 3,
    % 4, 9, 11) carry, with medium heat transfer, the local polymer stress
    % difference as a field in the same basis (imr.physics.zener_field), for the medium heating
    zener_field = any(stress == [3 4 9 11]) && medtherm == 1;
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
    % Kirchhoff transform of K = alpha_m0*T + beta_m0, the frozen mixture
    % conductivity, theta = alpha_m0*(T^2 - 1)/2 + beta_m0*(T - 1), inverted
    % as T = (sqrt((alpha_m0 + beta_m0)^2 + 2*alpha_m0*theta) - beta_m0)/alpha_m0
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
    % The local field and the stress integral must agree through their
    % integral relation (2*int_1^Inf sp/y dy = Z1/R^3); a zero field with
    % a nonzero Z1 gives the wrong initial medium heating.
    if zener_field
        init(izf) = f_zener_field_init(Szero);
    end

    % solver start
    imr.solvers.describe_case(radial, bubtherm, medtherm, 0, stress, spectral,...
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
            'alpha_m0',alpha_m0, 'beta_m0',beta_m0, 'nn',nn, 'Pb_star',Pb_star, ...
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
        if spectral_stress
            jp.ptt = ptt;
            jp.eps3 = eps3;
        end
        Sp = @(t,X) imr.numerics.jacobian_spectral(t,X,jp);
    else
        Sp = [];
    end
    system.initial = init;
    system.rhs = bubble;
    system.output_rhs = @(t,X) f_rhs(t,X,true);
    system.jacobian = Sp;
    system.outputs = @f_outputs;
    system.operatorSizes = [size(gA,1) size(mA,1)];

    % Reconstruct boundary values and preserve the public output tuple.
    % nargout: number of outputs requested by the caller. Expensive
    % reconstructions (temperature/medium-temperature history) are only
    % computed if their outputs are actually requested.
    function out = f_outputs(t,X,nargout_req)

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
        % Default to returning all outputs if nargout not provided
        if nargin < 3
            nargout_req = 8;
        end
        % only compute temperature if requested
        if bubtherm && nargout_req >= 5
            % extracting the Chebyshev coefficients
            a = X(:,ia)';
            T = imr.physics.thermal_inverse(gA*a,alpha_m0,beta_m0);
            if medtherm && nargout_req >= 6
                b = X(:,ib)';
                Tm = mA*b;
            end
        end

        % outputs: only populate outputs that caller requested
        out = cell(1,8);
        out{1} = t;
        out{2} = R;
        out{3} = Rdot;
        out{4} = P;
        if bubtherm == 1 && nargout_req >= 5
            out{5} = T;
        else
            out{5} = [];
        end
        if medtherm == 1 && bubtherm && nargout_req >= 6
            out{6} = Tm;
        else
            out{6} = ones(size(t,1),1);
        end
        out{7} = [];
        % the inputs of this run (imr.default_inputs names)
        out{8} = inputs;

    end

    % solver function
    function [dXdt,S,Sdot] = f_rhs(t,X,quiet)

        % showing output
        if progdisplay && (nargin < 3 || ~quiet)
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
        if nonnewtonian
            % fnu_m: viscosity factor at the medium grid points
            [~,intfnu,dintfnu,ddintfnu,fnu_m] = imr.physics.viscosity(nu_model, ...
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
            T = imr.physics.thermal_inverse(theta,alpha_m0,beta_m0);
            D = kapover*(alpha_m0*T.^2 + beta_m0*T)/P;

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
                [taugradu] = imr.physics.stress_power(stress,Req,R,Rdot,Ca, ...
                    Br,Re8,LAM,alphax,C2,iDRe,yT2,yT3,iyT3,iyT4,iyT6,fnu_m, ...
                    f_polymer_field(X));
                % surrounding temperature evolution equation
                Tmdot = advection + diffusion + taugradu;

                % enforce boundary condition and solve A*x=b problem with
                % the Sherman-Morrison update of the prefactored system
                Tmdot(end) = 0;
                s = alpha_m0*T(1) + beta_m0;
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
            [P, Pdot] = imr.physics.polytropic_pressure(R, Rdot, ...
                Pb_star-Pv_star, Pv_star, kappa);
        end

        % pressure waveform
        if wave_const
            Pf8 = ee;
            Pf8dot = 0;
        else
            [Pf8,Pf8dot] = imr.physics.far_field_pressure(t,pvarargin);
        end

        % stress equation evolution
        if spectral_stress
            [S,Sdot,Z1dot,Z2dot] = f_spectral_stress(R,Rdot,X);
        else
            [S,Sdot,Z1dot,Z2dot] = imr.physics.material_stress(stress,X,Req,R,Ca,De,Re8, ...
                Rdot,alphax,ivisco1,ivisco2,LAM,intfnu,dintfnu,iDRe,C2);
        end

        % bubble wall evolution / acceleration
        [Rddot] = imr.physics.radial_acceleration(radial, P, Pdot, Pf8, Pf8dot, iWe, R, Rdot, S, ...
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

    % Initial Zener field coefficients for the initial Z1 = Szero, mirroring
    % the stopgap of imr.solve_nodal (f_zener_field_init there). From rest
    % (Z1 = 0) the field is 0. Otherwise the field takes the shape of the
    % elastic stress difference se at the initial radius, scaled so that
    % 2*int_1^Inf sp/y dy = Z1/R^3 holds, then projected onto the Nv stress
    % modes with sCI, the same values-to-coefficients map f_zener_rhs uses.
    % Without elasticity (Ca = Inf) or at R = Req, where se = 0, it takes
    % the viscous shape y^-3 instead. The basis T_n(z) - 1 (n = 1..Nv)
    % already vanishes at z = 1 (y -> Inf), so no endpoint needs forcing.
    function e0 = f_zener_field_init(Z1)
        % Initialize the Zener field coefficients from the stress integral Z1.
        % DESIGN CHOICE: among infinitely many spatial distributions that
        % satisfy the integral constraint 2*int_1^inf sp/y dy = Z1/R^3, this
        % function picks the elastic restoring profile se(y) (derived from the
        % elasticity PDE at rest). This choice is canonical because:
        %   1. The field will immediately evolve via f_zener_rhs, which solves
        %      the polymer stress PDE; starting with a shape that matches the
        %      elasticity equations minimizes transients.
        %   2. Consistency with build_nodal.m, which uses the same shape.
        %   3. No startup noise in initial heating/stress.
        % An alternative (e.g., uniform, power-law, or viscous shape) would
        % also satisfy the integral constraint but would reshape on the first
        % RHS call. The elastic shape is preferred because it is already the
        % shape the field PDE "wants" to be in.
        e0 = zeros(Nf,1);
        if Z1 == 0
            return
        end
        Rst = Req/Rzero;
        yV2 = yV.^2;
        x2 = (yV3 - 1 + Rst^3).^(2/3);
        lam2 = yV2./x2;
        lam4i = x2.^2./yV2.^2;
        if stress == 11
            lam4 = lam2.^2;
            lam2i = x2./yV2;
            se = (lam4i - lam2)/Ca - (C2/Ca)*(lam4 - lam2i);
            Se = -(5 - 4*Rst - Rst^4)/(2*Ca) - (C2/Ca)*(2/Rst - Rst^2 - 1);
        else
            a4 = alphax*any(stress == [4 9]);
            x = lam4i + 2*lam2 - 3;
            factor = 1 + a4*x;
            Se = (3*a4 - 1)*(5 - Rst^4 - 4*Rst)/(2*Ca) + (2*a4/Ca)*(27/40 + ...
                (1/8)*Rst^8 + (1/5)*Rst^5 + Rst^2 - 2/Rst);
            if stress == 9
                factor = factor + 0.5*a4^2*x.^2;
                Se = Se + (a4^2/Ca)*(-2909/180 + 12/Rst - (4/3)/Rst^3 - 6*Rst^2 + ...
                    4*Rst^3 + (9/4)*Rst^4 - (6/5)*Rst^5 + (2/3)*Rst^6 - ...
                    (3/4)*Rst^8 + (1/9)*Rst^9 + (1/12)*Rst^12 + 4*log(Rst));
            end
            % elastic stress integral 2*int se/y dy, as in imr.physics.material_stress
            se = (lam4i - lam2)/Ca.*factor;
        end
        if isfinite(Ca) && abs(Se) > 0
            sp = Z1/(Rzero^3*Se)*se;
        else
            sp = 1.5*Z1/Rzero^3./yV3;
        end
        e0 = sCI*sp;
    end

    % time derivative of the Chebyshev coefficients e of the Zener field
    % sp = sCA*e at the stress collocation points (imr.physics.zener_field), with
    % the advection -ydot*d(sp)/dy = -dz/dt*d(sp)/dz
    function edot = f_zener_rhs(R,Rdot,X)
        e = X(izf);
        dzdt = dzdy.*(Rdot/R).*yV_adv;
        spdot = imr.physics.zener_field(sCA*e, -dzdt.*(sCAd*e), yV, R, Rdot, Req, Ca, ...
            De, Re8, LAM, stress, alphax, C2);
        edot = sCI*spdot;
    end

    % UCM/Oldroyd-B, PTT and Giesekus stress fields in the spectral
    % representation (Warnez & Johnsen, "Numerical modeling of bubble
    % dynamics in viscoelastic media with relaxation," Phys. Fluids 27,
    % 063103 (2015), https://doi.org/10.1063/1.4922598, Eq. 13). The
    % polymer stresses obey, along material points with a = R^2*Rdot/r^3,
    %   D(tau_rr)/Dt = -tau_rr/De*exp(ptt*Re8*De*(tau_rr+2*tau_tt))
    %       - 4*a*tau_rr - eps3*Re8*tau_rr^2 - 4*etap*a/De
    %   D(tau_tt)/Dt = -tau_tt/De*exp(ptt*Re8*De*(tau_rr+2*tau_tt))
    %       + 2*a*tau_tt - eps3*Re8*tau_tt^2 + 2*etap*a/De
    % with etap = (1-LAM)/Re8. stress = 5 (UCM, Oldroyd-B) has ptt = 0 and
    % eps3 = 0, recovering the original linear closure exactly. stress = 6
    % (PTT) sets ptt = 1, their linear (eps2 = 1) exponential closure.
    % stress = 7 (Giesekus) keeps ptt = 0 and takes eps3 in (0, 0.5] as
    % their mobility factor. In the moving coordinate z the material
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
        % ptt and eps3 are fixed flags/parameters, not state-dependent: skip
        % each term entirely when inactive, rather than forming ptt*Re8 or
        % eps3*Re8, so an inviscid run (Re8 = Inf) with the other model
        % does not turn an inactive term's 0*Inf into NaN
        if ptt
            exptau = exp(Re8*De*(trr + 2*ttt));
        else
            exptau = 1;
        end
        if eps3
            gies_rr = eps3*Re8*trr.^2;
            gies_tt = eps3*Re8*ttt.^2;
        else
            gies_rr = 0;
            gies_tt = 0;
        end
        trrdot = -dzdt.*(sCAd*c) - exptau.*trr/De - 4*av.*trr - ...
            gies_rr - 4*etap*av/De;
        tttdot = -dzdt.*(sCAd*d) - exptau.*ttt/De + 2*av.*ttt - ...
            gies_tt + 2*etap*av/De;
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
        if nonnewtonian
            [~,~,~,~,fnu_loc] = imr.physics.viscosity(nu_model,Rdotv,Rv,v_a,v_nc, ...
                v_lambda_star,yT);
        end
        tg = imr.physics.stress_power(stress,Req,Rv,Rdotv,Ca,Br,Re8,LAM, ...
            alphax,C2,iDRe,yT2,yT3,iyT3,iyT4,iyT6,fnu_loc,f_polymer_field(X));
        tg = tg.*ones(size(yT));
    end

    % radial acceleration and stress evolution as a function of
    % v = [R, Rdot, P, Pdot, Z], mirroring f_rhs
    function out = f_closure(t,v,X)
        X(iz) = v(5:end);
        [intf_loc,dintf_loc,ddintf_loc] = f_nu(v(1),v(2));
        [S_loc,Sdot_loc,Z1dot_loc,Z2dot_loc] = imr.physics.material_stress(stress,X,Req,v(1), ...
            Ca,De,Re8,v(2),alphax,ivisco1,ivisco2,LAM, ...
            intf_loc,dintf_loc,iDRe,C2);
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
        if nonnewtonian
            [~,intf_loc,dintf_loc,ddintf_loc] = imr.physics.viscosity(nu_model, ...
                Rdotv,Rv,v_a,v_nc,v_lambda_star);
        end
    end

    function Rddot_loc = f_rddot(t,v,S_loc,Sdot_loc,ddintf_loc)
        if wave_const
            Pf8_loc = ee;
            Pf8dot_loc = 0;
        else
            [Pf8_loc,Pf8dot_loc] = imr.physics.far_field_pressure(t,pvarargin);
        end
        Rddot_loc = imr.physics.radial_acceleration(radial, v(3), v(4), Pf8_loc, Pf8dot_loc, ...
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
