% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Scheme-specific operators, closures and reconstruction for one checked case.
% Public calls and time integration share imr.solvers.run.
%
% BUILDER STRUCTURE (roughly in order):
%  1. Parameter unpacking: extract model equations (radial, bubtherm, etc.),
%     solver options (scheme, jacobian type, grid sizes), and the initial
%     state (Rzero, Pb_star, Szero).
%  2. Operator setup: discretization matrices for temperature and medium
%     temperature fields (D1/D2 for nodal derivatives, based on scheme).
%  3. Physical constants: conductivity coefficients (alpha_g, beta_g,
%     alpha_v, beta_v), saturation pressure constants, gas constants.
%  4. Thermal/mass-transfer operators and coefficients, especially the
%     Kirchhoff transform and wall boundary gradients.
%  5. State layout and index arrays: maps state vector to R, Rdot, P,
%     bubble temperature, medium temperature, vapor concentration, stress.
%  6. Initial condition assembly: populate init vector from prepared values.
%  7. Wall closure helpers: functions f_wall_theta_exact, f_wall_mass_closure
%     and their supporting machinery (fzero options, wall constant struct).
%  8. Jacobian setup: conditionally build jp (parameters for analytic Jacobian)
%     or set up finite-difference Jacobian.
%  9. System assembly: package system.initial, system.rhs, system.jacobian,
%     system.outputs into the struct returned to run.m.
% 10. Nested functions: f_outputs, f_reset_output_wall, f_rhs and all RHS
%     helpers (f_taugradu, f_closure, etc.).
%
% Variables are named to hint at their physical meaning. For example:
%   D1/D2   = first and second-derivative matrices (spatial operators)
%   dmatrix = D1 mapped to the working grid; ddmatrix = D2 mapped.
%   alpha/beta = Kirchhoff transform coefficients for thermal conductivity.
%   theta   = transformed bubble temperature (state variable).
%   T       = physical temperature (reconstructed from theta).
%   Tm      = medium temperature.
%   kv      = vapor mole fraction (dimensionless).
%   S/Z1    = stress integral (large scalar); Z1dot its time derivative.
%   cp      = polymer stress (state variable if Zener model).
%   yT*     = spatial grids and related factors (yT = y*T where T is the
%     thermal diffusion length, used for non-dimensionalization).
function system = build_nodal(prepared)
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
    masstrans       = prepared.model.masstrans;

    % solver options
    spectral        = prepared.numerics.spectral;
    Nt              = prepared.numerics.Nt;
    Mt              = prepared.numerics.Mt;
    % Lv            = prepared.numerics.Lv, only used by the spectral solver
    Lt              = prepared.numerics.Lt;
    jacobian        = prepared.numerics.jacobian;
    scheme          = prepared.numerics.scheme;

    % dimensionless initial conditions
    Rzero           = prepared.initial.R;
    Rdotzero        = prepared.initial.Rdot;
    Pb_star         = prepared.initial.P;
    Pv_star         = prepared.initial.Pv;
    Req             = prepared.initial.Req;
    alphax          = prepared.stress.alphax;
    C2              = prepared.stress.C2;

    % dimensionless initial stress
    Szero           = prepared.initial.stress;

    % time span options
    tspan           = prepared.time;

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

    % dimensionless mass transfer
    Fom             = prepared.mass.Fom;
    kv0             = prepared.mass.kv0;
    Rv_star         = prepared.mass.Rv_star;
    Rg_star         = prepared.mass.Ra_star;
    L_heat_star     = prepared.mass.L_heat_star;
    thetav          = prepared.mass.thetav;

    % pre_process

    % Build only enabled field operators. Empty grids do not enter physics.
    y = [];
    xi = [];
    dmatrix_Tkv = sparse(0,0);
    ddmatrix_Tkv = sparse(0,0);
    dmatrix_Tm = sparse(0,0);
    ddmatrix_Tm = sparse(0,0);
    if bubtherm
        [y,dmatrix_Tkv] = imr.numerics.differentiation(Nt,1,0,scheme);
        [~,ddmatrix_Tkv] = imr.numerics.differentiation(Nt,2,0,scheme);
    end
    if medtherm
        [xi,dmatrix_Tm] = imr.numerics.differentiation(Mt,1,1,scheme);
        [~,ddmatrix_Tm] = imr.numerics.differentiation(Mt,2,1,scheme);
    end

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

    % precomputations
    sam = 1 + GAMa;
    no = (nstate-1)/nstate;
    nog = (nstate-1)/2;
    kapover = (kappa-1)/kappa;
    % Kirchhoff transform of the conductivity K = alpha*T + beta,
    % theta = alpha*(T^2 - 1)/2 + beta*(T - 1), with the inverse
    % T = (sqrt((alpha + beta)^2 + 2*alpha*theta) - beta)/alpha, where the
    % square root is K itself. alpha + beta is K at T = 1, which is not 1
    % in general (the conductivities are scaled by the mean of the gas and
    % vapor values). The map of the bubble without mass transfer uses the
    % constant mixture coefficients of kv0 (as in f_theta_of_T)
    alpha_m0 = kv0.*alpha_v + (1-kv0).*alpha_g;
    beta_m0 = kv0.*beta_v + (1-kv0).*beta_g;
    % the thermal-only wall balance (masstrans = 0) uses the same frozen
    % mixture, so the wall shares one conductivity with the bulk equation
    K1sq_m0 = (alpha_m0 + beta_m0)^2;
    om_kv0 = 1-kv0;
    % dimensionless saturated vapor pressure, Pv(T) = C1_pv*exp(C2_pv/T)
    % with Pv(1) = Pv_star
    C1_pv = Pv_star*exp(thetav);
    C2_pv = -thetav;
    Rvg_ratio = Rv_star / Rg_star;
    Rva_diff = Rv_star - Rg_star;
    % bubble-wall flux balance coefficients: the wall gradient rows of the
    % medium (first entry at the wall) and of the bubble (first entry at
    % the wall, then inward), with the node indices of the other entries
    grad_Tm_coeff = [];
    grad_Trans_coeff = [];
    grad_C_coeff = [];
    iwall_m = [];
    iwall_b = [];
    if bubtherm
        if scheme == 0
            coeff = [-1.5 2 -0.5];
            deltaY = 1/(Nt-1);
            grad_Trans_coeff = -coeff*chi/deltaY;
            grad_C_coeff = -coeff*Fom*L_heat_star/deltaY;
            iwall_b = [Nt-1 Nt-2];
        else
            grad_Trans_coeff = chi*full(dmatrix_Tkv(end,end:-1:1));
            grad_C_coeff = Fom*L_heat_star*full(dmatrix_Tkv(end,end:-1:1));
            iwall_b = Nt-1:-1:1;
        end
    end
    if medtherm
        if scheme == 0
            deltaYm = -2/(Mt-1);
            grad_Tm_coeff = 2*chi*iota/deltaYm*[-1.5 2 -0.5];
            iwall_m = [2 3];
        else
            grad_Tm_coeff = 2*chi*iota*full(dmatrix_Tm(1,:));
            iwall_m = 2:Mt;
        end
    end

    % State indices are shared by initialization, RHS, Jacobian and outputs.
    Mt = numel(layout.medium);
    ibubtherm = layout.bubble;
    imedtherm = layout.medium;
    imass = layout.vapor;
    ivisco1 = layout.stress1;
    ivisco2 = layout.stress2;
    izf = layout.polymer;
    zener_field = ~isempty(izf);

    % initial condition assembly

    init = zeros(layout.count,1);
    init([layout.R layout.Rdot layout.P]) = [Rzero Rdotzero Pb_star]';
    init(layout.medium) = 1;
    init(layout.vapor) = kv0;
    init(layout.stress) = Szero;
    if zener_field
        init(layout.polymer) = f_zener_field_init(Szero);
    end

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
        'alpha_g',alpha_g, 'alpha_v',alpha_v, 'beta_g',beta_g, ...
        'beta_v',beta_v, 'C1',C1_pv, 'C2',C2_pv, 'Rvg',Rvg_ratio, ...
        'Rva',Rva_diff, 'Rg',Rg_star);

    % solver start
    if scheme == 0
        scheme_name = 'nodal, finite difference';
    else
        scheme_name = 'nodal, Chebyshev';
    end
    imr.solvers.describe_case(radial, bubtherm, medtherm, masstrans, stress, spectral,...
        nu_model, eps3, Pv_star, Re8, De, Ca, LAM, scheme_name, DRe, iDRe, ...
        v_lambda_star, v_nc, v_a);
    bubble = @f_rhs;
    % Jacobian for the stiff integrators: analytic where implemented,
    % batched finite differences (jacobian = 3), otherwise the sparsity
    % pattern
    % the stress states and the Zener field, for the sparsity pattern
    istress = layout.allStress;
    if jacobian == 3
        jmask = full(imr.numerics.jacobian_pattern(numel(init), ibubtherm, imedtherm, imass, ...
            istress, scheme)) ~= 0;
        Sp = @f_jacobian_fd;
    elseif jacobian == 2 && bubtherm
        jp = struct('ith',ibubtherm, 'itm',imedtherm, 'ikv',imass, ...
            'iz',[ivisco1 ivisco2], 'izf',izf, 'medtherm',medtherm, ...
            'masstrans',masstrans, 'D1',dmatrix_Tkv, 'D2',ddmatrix_Tkv, ...
            'Dm1',dmatrix_Tm, 'Dm2',ddmatrix_Tm, 'y',y, 'xi',xi, 'yT',yT, ...
            'yT2',yT2, 'yT3',yT3, 'Lt',Lt, 'Foh',Foh, 'chi',chi, ...
            'kappa',kappa, 'kapover',kapover, 'alpha_g',alpha_g, ...
            'beta_g',beta_g, 'alpha_v',alpha_v, 'beta_v',beta_v, 'kv0',kv0, ...
            'Fom',Fom, 'Rv',Rv_star, 'Rg',Rg_star, 'Rva',Rva_diff, ...
            'Rvg',Rvg_ratio, 'C1_pv',C1_pv, 'C2_pv',C2_pv, 'wall',wallc, ...
            'closure',@f_closure, 'taugradu',@f_taugradu, ...
            'theta_bw',@f_theta_bw_of_X, 'zener_field',zener_field, ...
            'zener_rhs',@f_zener_rhs, 'Br',Br, 'De',De);
        Sp = @(t,X) imr.numerics.jacobian_nodal(t,X,jp);
    elseif jacobian >= 1
        Sp = imr.numerics.jacobian_pattern(numel(init), ibubtherm, imedtherm, imass, ...
            istress, scheme);
    else
        Sp = [];
    end
    system.initial = init;
    system.rhs = bubble;
    system.output_rhs = @(t,X) f_rhs(t,X,true);
    system.reset_output_rhs = @f_reset_output_wall;
    system.jacobian = Sp;
    system.outputs = @f_outputs;
    system.operatorSizes = [size(dmatrix_Tkv,1) size(dmatrix_Tm,1)];

    % Reconstruct boundary values and preserve the public output tuple.
    % nargout: number of outputs requested by the caller. Expensive
    % reconstructions (wall resolution, temperature history, vapor history)
    % are only computed if their outputs are actually requested.
    function out = f_outputs(t,X,nargout_req)

        % extract result
        R    = X(:,1);
        Rdot = X(:,2);
        P    = X(:,3);
        % Default to returning all outputs if nargout not provided
        if nargin < 3
            nargout_req = 8;
        end
        % Skip expensive thermal/vapor reconstructions if outputs are not requested
        if bubtherm && nargout_req >= 5
            theta = X(:,ibubtherm);
            if medtherm
                Tm = X(:,imedtherm);
                % the wall entries of the state are not integrated (their time
                % derivatives are zero), so the wall state is recovered from the
                % wall closure at every output time, as the right-hand side does
                % (the loop index is named it_out because a nested function shares
                % any variable name that the parent also uses)
                kv_w = zeros(numel(t),1);
                f_reset_output_wall(X(1,:).');
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
        if bubtherm && nargout_req >= 5
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
            else
                T = f_theta_of_T(theta,kv0);
            end
        end
        if medtherm && bubtherm && nargout_req >= 6
            Tm(:,1) = T(:,end);
        end

        % outputs assembly: only populate outputs that caller requested
        out = cell(1,8);
        out{1} = t;
        out{2} = R;
        out{3} = Rdot;
        out{4} = P;
        if bubtherm && nargout_req >= 5
            out{5} = T;
        else
            out{5} = [];
        end
        if medtherm == 1 && nargout_req >= 6
            out{6} = Tm;
        else
            out{6} = ones(size(t,1),1);
        end
        if masstrans == 1 && nargout_req >= 7
            out{7} = kv;
        else
            out{7} = [];
        end
        % the inputs of this run (imr.default_inputs names)
        out{8} = inputs;

    end

    % Each history starts with the same wall-root guess, rather than the
    % final integration/output cache. Later times retain their warm starts.
    function f_reset_output_wall(X)
        if bubtherm && medtherm
            theta_bw_guess = X(ibubtherm(end));
            kv_w_guess = NaN;
        end
    end

    % solver function
    function [dXdt,S,Sdot] = f_rhs(t,X,quiet)

        % X holds one state per column, so a batch of states (the
        % finite-difference Jacobian) costs about one call; every operation
        % below is elementwise across the columns, in the same order as for
        % a single state, so one column gives the same result bit for bit

        % showing output
        if progdisplay && (nargin < 3 || ~quiet)
            disp(t/tspan(end));
        end

        % extract solution

        % bubble wall radius
        R = X(1,:);
        % bubble wall velocity
        Rdot = X(2,:);
        % internal bubble pressure
        P = X(3,:);
        % auxiliary variable for internal bubble temperature
        if bubtherm
            theta = X(ibubtherm,:);
        end
        if masstrans
            kv = X(imass,:);
        end

        % empty output variables
        thetadot = [];
        Tmdot = [];
        kvdot = [];

        % boundary condition evaluation
        if medtherm
            % temperature in the material
            Tm = X(imedtherm,:);
        end
        if medtherm && masstrans
            % wall closure on the wall vapor mass fraction
            [theta_bw, kv_bw] = f_wall_mass_closure(Tm, theta, kv, P, R, t);
            theta(end,:) = theta_bw;
        elseif medtherm
            % closed-form wall temperature, fzero only as a fallback
            theta_bw = f_wall_theta_exact(Tm, theta);
            theta_bw_guess = theta_bw(1);
            theta(end,:) = theta_bw;
        end

        % Bubble temperature uses the same pure transform as its Jacobian.
        if masstrans
            alpha_T = kv.*alpha_v + (1-kv).*alpha_g;
            beta_T = kv.*beta_v + (1-kv).*beta_g;
            T = imr.physics.thermal_inverse(theta,alpha_T,beta_T);
        elseif bubtherm
            T = imr.physics.thermal_inverse(theta,alpha_m0,beta_m0);
        end

        % bubble equations of motion
        if bubtherm && masstrans
            % vapor concentration at the boundary. With the medium, the
            % converged closure value; it is set after T so that T(end), like
            % the closure, uses the mixture coefficient of the stored kv(end)
            if medtherm
                kv(end,:) = kv_bw;
            else
                kv(end,:) = f_kv_of_T(T(end,:),P);
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
            Pdot = 3./R.*(chi*(kappa-1)*dtheta(end,:)./R - kappa*P.*Rdot +...
                + kappa*P*Fom*Rv_star.*dkv(end,:)./(T(end,:).*R.*Rmix(end,:).*om_kv(end,:)));

            % bubble wall velocity equation
            Uvel = (chi./R*(kappa-1).*dtheta-y.*R.*Pdot/3)./(kappa*P) + Fom./R.*RDkv;

            % mixture thermal conductivity
            Kstar_g = alpha_g*T+beta_g;
            Kstar_v = alpha_v*T+beta_v;
            Kstar   = kv.*Kstar_v + om_kv.*Kstar_g;

            % temperature evolution inside the bubble
            nonlinear_term = (chi*ddtheta./R.^2+Pdot).*(kapover*Kstar.*T./P);
            % gas velocity relative to the moving grid
            u_rel = Uvel-y.*Rdot;
            advection_term = -dtheta.*u_rel./R;
            mass_diffusion = (Fom./(R.^2)).*wRmix.*dkv.*dtheta;
            thetadot = advection_term + nonlinear_term + mass_diffusion;
            % temperature at bubble wall, handled by flux boundary bc
            thetadot(end,:) = 0;

            % water vapor evolution inside the bubble
            nonlinear_diffusion = ...
                dkv.*(dtheta./(Kstar.*T)+RDkv);
            advection_term = u_rel./R.*dkv;
            kvdot = Fom./R.^2.*(ddkv - nonlinear_diffusion) - advection_term;
            % concentration at bubble wall, handled by flux boundary bc
            kvdot(end,:) = 0;

        elseif bubtherm
            % mixture thermal conductivity
            Kstar_g = alpha_g*T+beta_g;
            Kstar_v = alpha_v*T+beta_v;
            Kstar   = kv0.*Kstar_v + om_kv0.*Kstar_g;

            % bubble temperature gradients
            dtheta  = dmatrix_Tkv*theta;
            ddtheta = ddmatrix_Tkv*theta;

            % bubble pressure evolution equation
            Pdot = 3./R.*(chi*(kappa-1)*dtheta(end,:)./R - kappa*P.*Rdot);

            % internal bubble velocity equation
            Uvel = (chi./R*(kappa-1).*dtheta-y.*R.*Pdot/3)./(kappa*P);

            % bubble temperature evolution equation
            diffusion_term = (chi*ddtheta./R.^2+Pdot).*(kapover*Kstar.*T./P);
            advection_term = -dtheta.*(Uvel-y.*Rdot)./R;
            thetadot = advection_term + diffusion_term;
            thetadot(end,:) = 0;

        else
            % polytropic gas + vapor
            [P, Pdot] = imr.physics.polytropic_pressure(R, Rdot, ...
                Pb_star-Pv_star, Pv_star, kappa);
        end

        % surroundings equations of motion

        % viscous forces/Reynolds number
        if nonnewtonian
            % fnu_m: viscosity factor at the medium grid points
            [~,intfnu,dintfnu,ddintfnu,fnu_m] = imr.physics.viscosity(nu_model, ...
                Rdot,R,v_a,v_nc,v_lambda_star,yT);
        end

        % surrounding temperature
        if medtherm
            % boundary temperature
            Tm(1,:) = T(end,:);
            % temperature gradients
            dTm = dmatrix_Tm*Tm;
            ddTm = ddmatrix_Tm*Tm;
            % material temperature equations
            advection = xi1_sq./(Lt*R).*(Rdot./yT2.*om_yT3/2 +...
                Foh./R.*xi_adv).*dTm;
            diffusion = Foh./R.^2.*xi1_p4/Lt^2.*ddTm/4;
            % stress dissipation, with the local polymer stress of the
            % Zener models
            if zener_field
                sp = X(izf,:);
            else
                sp = [];
            end
            [taugradu] = imr.physics.stress_power(stress,Req,R,Rdot,Ca,Br, ...
                Re8,LAM,alphax,C2,iDRe,yT2,yT3,iyT3,iyT4,iyT6,fnu_m,sp);
            % surrounding temperature evolution equation
            Tmdot = advection + diffusion + taugradu;
            % sets boundary condition on temperature
            Tmdot(1,:) = 0;
            Tmdot(end,:) = 0;
        end

        % pressure waveform
        if wave_const
            Pf8 = ee;
            Pf8dot = 0;
        else
            [Pf8,Pf8dot] = imr.physics.far_field_pressure(t,pvarargin);
        end

        % stress equation evolution
        [S,Sdot,Z1dot,Z2dot] = imr.physics.material_stress(stress,X,Req,R,Ca,De,Re8, ...
            Rdot,alphax,ivisco1,ivisco2,LAM,intfnu,dintfnu,iDRe,C2);

        % bubble wall evolution / acceleration
        [Rddot] = imr.physics.radial_acceleration(radial, P, Pdot, Pf8, Pf8dot, iWe, R, Rdot, S, ...
            Sdot, Cstar, sam, no, GAMa, nstate, nog, hugoniot_s, JdotA, ...
            ddintfnu, iDRe);

        % local polymer stress of the Zener models, for the heating
        if zener_field
            spdot = f_zener_rhs(R,Rdot,sp);
        else
            spdot = [];
        end

        % output assembly
        dXdt = [Rdot;
        Rddot;
        Pdot;
        thetadot;
        Tmdot;
        kvdot;
        Z1dot;
        Z2dot;
        spdot];

    end
    % end of solver

    % finite-difference Jacobian from one batched right-hand-side call: the
    % state and its n perturbations, one per column. The step is
    % sqrt(eps)*max(|X|, 1) (the states are dimensionless and of order
    % one), rounded so that the perturbed state is exactly X + del, and
    % entries outside the sparsity pattern are set to zero (they are zero
    % up to rounding in the batched products). J is returned dense: it is
    % computed dense, and dense LU is fast at these sizes
    function J = f_jacobian_fd(t, X)
        del = sqrt(eps)*max(abs(X), 1);
        del = (X + del) - X;
        F = f_rhs(t, [X, repmat(X, 1, numel(X)) + diag(del)]);
        J = ((F(:,2:end) - F(:,1))./del.').*jmask;
    end

    % time derivative of the Zener field sp at the medium nodes
    % (imr.physics.zener_field). The advection -ydot*d(sp)/dy is the Rdot part of the
    % medium temperature advection, and the far-field node stays at 0
    function spdot = f_zener_rhs(Rv,Rdotv,sp)
        adv = xi1_sq./(Lt*Rv).*(Rdotv./yT2.*om_yT3/2).*(dmatrix_Tm*sp);
        spdot = imr.physics.zener_field(sp, adv, yT, Rv, Rdotv, Req, Ca, De, Re8, LAM, ...
            stress, alphax, C2);
        spdot(end,:) = 0;
    end

    % initial Zener field for the initial Z1 = Szero. From rest (Z1 = 0)
    % the field is 0. Otherwise (collapse, or a given stress0) it is a
    % STOPGAP: the history of the field before t = 0 is not computed, and
    % the field takes the shape of the elastic stress difference se of the
    % initial radius, scaled so that 2*int_1^Inf sp/y dy = Z1/R^3 holds.
    % Without elasticity (Ca = Inf) or at R = Req, where se = 0, it takes
    % the viscous shape y^-3, for which 2*int_1^Inf y^-4 dy = 2/3
    function sp = f_zener_field_init(Z1)
        sp = zeros(Mt,1);
        if Z1 == 0
            return
        end
        Rst = Req/Rzero;
        x2 = (yT3 - 1 + Rst^3).^(2/3);
        lam2 = yT2./x2;
        lam4i = x2.^2.*iyT4;
        if stress == 11
            lam4 = lam2.^2;
            lam2i = x2./yT2;
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
                    9*Rst + (4/3)*Rst^3 + (9/4)*Rst^4 - (6/5)*Rst^5 + (2/3)*Rst^6 - ...
                    (3/4)*Rst^8 + (1/9)*Rst^9 + (1/12)*Rst^12 + 4*log(Rst));
            end
            % elastic stress integral 2*int se/y dy, as in imr.physics.material_stress
            se = (lam4i - lam2)/Ca.*factor;
        end
        if isfinite(Ca) && abs(Se) > 0
            sp = Z1/(Rzero^3*Se)*se;
        else
            sp = 1.5*Z1/Rzero^3*iyT3;
        end
        sp(end) = 0;
    end

    % medium dissipation as a function of (R, Rdot)
    function tg = f_taugradu(Rv,Rdotv,X)
        fnu_loc = 0;
        if nonnewtonian
            [~,~,~,~,fnu_loc] = imr.physics.viscosity(nu_model,Rdotv,Rv,v_a,v_nc, ...
                v_lambda_star,yT);
        end
        if zener_field
            sp_loc = X(izf);
        else
            sp_loc = [];
        end
        tg = imr.physics.stress_power(stress,Req,Rv,Rdotv,Ca,Br,Re8,LAM, ...
            alphax,C2,iDRe,yT2,yT3,iyT3,iyT4,iyT6,fnu_loc,sp_loc);
        tg = tg.*ones(size(yT));
    end

    % radial acceleration and stress evolution as a function of
    % v = [R, Rdot, P, Pdot, Z], mirroring f_rhs
    function out = f_closure(t,v,X)
        Rv = v(1);
        Rdotv = v(2);
        X([ivisco1 ivisco2]) = v(5:end);
        intf_loc = 0;
        dintf_loc = 0;
        ddintf_loc = 0;
        if nonnewtonian
            [~,intf_loc,dintf_loc,ddintf_loc] = imr.physics.viscosity(nu_model, ...
                Rdotv,Rv,v_a,v_nc,v_lambda_star);
        end
        if wave_const
            Pf8_loc = ee;
            Pf8dot_loc = 0;
        else
            [Pf8_loc,Pf8dot_loc] = imr.physics.far_field_pressure(t,pvarargin);
        end
        [S_loc,Sdot_loc,Z1dot_loc,Z2dot_loc] = imr.physics.material_stress(stress,X,Req,Rv, ...
            Ca,De,Re8,Rdotv,alphax,ivisco1,ivisco2,LAM, ...
            intf_loc,dintf_loc,iDRe,C2);
        Rddot_loc = imr.physics.radial_acceleration(radial, v(3), v(4), Pf8_loc, Pf8dot_loc, ...
            iWe, Rv, Rdotv, S_loc, Sdot_loc, Cstar, sam, no, GAMa, nstate, ...
            nog, hugoniot_s, JdotA, ddintf_loc, iDRe);
        out = vertcat(Rddot_loc, Z1dot_loc, Z2dot_loc);
    end

    % solver functions

    % temperature as a function of theta, with the mixture coefficients of kv
    function Tw = f_theta_of_T(theta_w,kv)
        alpha_m  = kv.*alpha_v + (1-kv).*alpha_g;
        beta_m  = kv.*beta_v + (1-kv).*beta_g;
        Tw = imr.physics.thermal_inverse(theta_w,alpha_m,beta_m);
    end

    % concentration at the bubble wall
    function kv_w = f_kv_of_T(Tw,pressure)
        theta_var = Rvg_ratio*(pressure./(imr.physics.saturation_pressure(Tw,C1_pv,C2_pv))-1);
        kv_w = 1./(1+theta_var);
    end

    % bubble-wall state with vapor mass transfer. The unknown is the wall
    % vapor mass fraction kv_w, confined to (0, 1). For a trial kv_w the
    % saturation condition gives the wall temperature and the Kirchhoff
    % transform gives theta_w in closed form (f_wall_mass_residual).
    % Warm start: Newton with the analytic derivative from the previous root
    % kv_w_guess, accepted once a step falls below warm_tol with every
    % iterate inside kv_bracket, within n_warm steps. Otherwise (first call,
    % divergence, a non-finite value) the bracketed solve
    % (f_wall_mass_bracket), which needs no initial guess. The inputs may
    % hold one state per column (a batch): the warm Newton then runs on all
    % columns at once, a column that fails it takes the bracketed solve on
    % its own, and the root of the first column becomes the next guess.
    function [theta_w, kv_w] = f_wall_mass_closure(Tm_in, th_in, kv_in, Pw, Rw, tw)
        % mixture coefficients from the stored wall value kv_in(end), not the
        % trial kv_w. Since kvdot(end) = 0 the stored value is the initial
        % kv0. Deliberate carry-over of the existing behavior, to be
        % revisited separately.
        alpha_m = kv_in(end,:)*alpha_v + (1 - kv_in(end,:))*alpha_g;
        % K = alpha_m*T + beta_m, so that
        % theta = 0.5*alpha_m*(T^2 - 1) + beta_m*(T - 1) is the exact inverse
        % of T = (sqrt((alpha_m + beta_m)^2 + 2*alpha_m*theta) - beta_m)/alpha_m
        beta_m = kv_in(end,:)*beta_v + (1 - kv_in(end,:))*beta_g;
        % contractions of the wall rows with the other nodes, independent of
        % kv_w. Order matters: Tm(iwall_m) runs outward from the wall
        % (Tm(2:3) for finite differences), while theta(iwall_b) and
        % kv(iwall_b) run inward from the wall (theta(end-1:-1:end-2) and
        % kv(end-1:-1:end-2) for finite differences), matching the order
        % of the coefficients in grad_Tm_coeff, grad_Trans_coeff, grad_C_coeff
        flux_fixed = grad_Tm_coeff(2:end)*Tm_in(iwall_m,:) + ...
            grad_Trans_coeff(2:end)*th_in(iwall_b,:);
        carrier_fixed = grad_C_coeff(2:end)*kv_in(iwall_b,:);
        % one column of parameters per state
        wp = [Pw;
        alpha_m;
        beta_m;
        flux_fixed;
        carrier_fixed];
        nw = size(wp, 2);

        % warm start from the previous root, on the columns still iterating
        kv_w = NaN(1, nw);
        if isfinite(kv_w_guess)
            x = kv_w_guess*ones(1, nw);
            act = true(1, nw);
            for kw = 1:n_warm
                ia = find(act);
                [res, ~, dres] = f_wall_mass_residual(x(ia), wp(:,ia));
                step = res./dres;
                xa = x(ia) - step;
                x(ia) = xa;
                inside = isfinite(xa) & xa > kv_bracket(1) & xa < kv_bracket(2);
                done = inside & abs(step) <= warm_tol;
                kv_w(ia(done)) = xa(done);
                act(ia(done | ~inside)) = false;
                if ~any(act)
                    break
                end
            end
        end

        % bracketed solve for the columns the warm start did not settle
        for jw = find(~isfinite(kv_w))
            kv_w(jw) = f_wall_mass_bracket(wp(:,jw), Rw(jw), kv_in(end,jw), tw);
        end
        kv_w_guess = kv_w(1);
        [~, theta_w] = f_wall_mass_residual(kv_w, wp);
    end

    % bracketed solve of the wall balance for one state (one column wp):
    % n_rounds rounds of n_grid residual evaluations, one vectorized call
    % each, narrow kv_bracket to the grid interval with the sign change (the
    % first round includes the bracket ends), and n_polish Newton steps
    % follow
    function kv_w = f_wall_mass_bracket(wp, Rw, kv_end, tw)
        lo = kv_bracket(1);
        hi = kv_bracket(2);
        for kr = 1:n_rounds
            kv_grid = linspace(lo, hi, n_grid);
            r_grid = f_wall_mass_residual(kv_grid, wp);
            if kr == 1
                r_lo = r_grid(1);
                r_hi = r_grid(end);
                if ~(isfinite(r_lo) && isfinite(r_hi)) || sign(r_lo) == sign(r_hi)
                    error(['imr.solve_nodal: the mass-transfer wall balance does ' ...
                        'not change sign on kv_w in [%g, 1 - %g] at t = %g ' ...
                        '(R = %g, P = %g, stored kv(end) = %g, residuals %g ' ...
                        'and %g). This indicates a physical or indexing ' ...
                        'problem in the state, not a solver problem.'], lo, ...
                        1 - hi, tw, Rw, wp(1), kv_end, r_lo, r_hi);
                end
            end
            ic = find(sign(r_grid(1:end-1)) ~= sign(r_grid(2:end)), 1);
            lo = kv_grid(ic);
            hi = kv_grid(ic+1);
        end
        kv_w = (lo + hi)/2;
        for kp = 1:n_polish
            [res, ~, dres] = f_wall_mass_residual(kv_w, wp);
            % stay inside the final grid interval
            kv_w = min(max(kv_w - res/dres, lo), hi);
        end
    end

    % mass-transfer wall balance and its derivative as functions of the wall
    % vapor mass fraction kv_w (elementwise, kv_w may be a vector),
    % wp = [P, alpha_m, beta_m, flux_fixed, carrier_fixed] as a column, or
    % one column per entry of kv_w. It is the balance
    % of imr.physics.wall_mass_balance multiplied through by
    % (kv_w*Rva_diff + Rg_star)*Tw*(1 - kv_w) > 0, which has the same root
    % and no reciprocal that vanishes as kv_w -> 1
    function [r, theta_w, dr] = f_wall_mass_residual(kv_w, wp)
        Pw = wp(1,:);
        % vapor partial pressure, and the wall temperature from inverting
        % Pv(Tw) = C1_pv*exp(C2_pv/Tw) (C2_pv < 0)
        q = kv_w*Rvg_ratio + 1 - kv_w;
        ps = Pw.*kv_w*Rvg_ratio./q;
        lg = log(ps/C1_pv);
        Tw = C2_pv./lg;
        theta_w = 0.5*wp(2,:).*(Tw.^2 - 1) + wp(3,:).*(Tw - 1);
        mix = kv_w*Rva_diff + Rg_star;
        denom = mix.*Tw.*(1 - kv_w);
        thermal = grad_Tm_coeff(1)*Tw + grad_Trans_coeff(1)*theta_w + wp(4,:);
        r = denom.*thermal + Pw.*(grad_C_coeff(1)*kv_w + wp(5,:));
        if nargout > 2
            % d/d(kv_w) through ps, Tw and theta_w
            dps = Pw*Rvg_ratio./q.^2;
            dTw = -C2_pv./lg.^2.*dps./ps;
            dth = (wp(2,:).*Tw + wp(3,:)).*dTw;
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
    % residual of the heat flux balance. The wall uses the same frozen
    % mixture conductivity (alpha_m0, beta_m0) as the bulk equation, so a
    % run with vapor present and masstrans = 0 is thermally consistent.
    function theta_w = f_bubble_wall_thermal_bc(theta_bw, Tm_in, th_in)
        Tm_trans = Tm_in(iwall_m);
        T_trans = th_in(iwall_b);
        Tw_prelim = imr.physics.thermal_inverse(theta_bw,alpha_m0,beta_m0);
        theta_w = grad_Tm_coeff * vertcat(Tw_prelim, Tm_trans) + ...
            grad_Trans_coeff * vertcat(theta_bw, T_trans);
    end

    % Solve the wall flux balance in d = Tw-1 rather than conductivity u.
    % theta = (alpha_m0+beta_m0)*d + alpha_m0*d^2/2 gives
    % A*d^2 + B*d + C = 0. This avoids subtracting u^2-K(1)^2 and
    % dividing by a tiny alpha_m0. Discard roots with negative conductivity
    % and choose the theta closest to the previous wall value.
    function theta_bw = f_wall_theta_exact(Tm_in, th_in)
        c1 = grad_Tm_coeff(1);
        c2 = grad_Trans_coeff(1);
        K = grad_Tm_coeff(2:end)*Tm_in(iwall_m,:) + ...
            grad_Trans_coeff(2:end)*th_in(iwall_b,:);
        A = 0.5*alpha_m0*c2;
        B = c1+c2*(alpha_m0+beta_m0);
        C = c1+K;
        if A == 0
            d1 = -C/B;
            d2 = NaN(size(C));
        else
            disc = B^2-4*A*C;
            sq = sqrt(max(disc,0));
            if B >= 0
                q = -0.5*(B+sq);
            else
                q = -0.5*(B-sq);
            end
            d1 = q/A;
            d2 = C./q;
            d1(disc < 0) = NaN;
            d2(disc < 0) = NaN;
        end
        ok1 = isfinite(d1) & alpha_m0*(1+d1)+beta_m0 >= 0;
        ok2 = isfinite(d2) & alpha_m0*(1+d2)+beta_m0 >= 0;
        th1 = d1.*(alpha_m0+beta_m0+0.5*alpha_m0*d1);
        th2 = d2.*(alpha_m0+beta_m0+0.5*alpha_m0*d2);
        pick2 = ok2 & (~ok1 | abs(th2 - theta_bw_guess) < abs(th1 - theta_bw_guess));
        theta_bw = th1;
        theta_bw(pick2) = th2(pick2);
        for jw = find(~(ok1 | ok2))
            theta_bw(jw) = f_wall_theta_fallback(c1, c2, Tm_in(:,jw), th_in(:,jw));
        end
    end

    % fzero on the thermal-only balance of one state, for a state where the
    % quadratic has no admissible root
    function theta_bw = f_wall_theta_fallback(c1, c2, Tm_in, th_in)
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
            theta_bw = (u^2 - K1sq_m0)/(2*alpha_m0);
        end
    end

end
