% file f_imr_rhs_diagnostics.m
% brief evaluates the IMR RHS without time integration
function diag = f_imr_rhs_diagnostics(q, varargin)
%F_IMR_RHS_DIAGNOSTICS Evaluate radial and modal accelerations.
%
%   diag = f_imr_rhs_diagnostics(q, name, value, ...)
%
%   q = [R; Rdot; Pb_state; epnm(:); epnmd(:)]
%
%   q may also be a matrix whose columns are states. For the polytropic gas
%   branch, Pb_state is the constant gas-pressure state used in
%
%       P = (Pb_state - Pv_star)*(1/R)^(3*kappa) + Pv_star.
%
%   Passing NaN uses the Pb_state implied by f_call_params. The name/value
%   arguments are the same arguments passed to f_imr_fd.
%
%   Optional diagnostic-only names:
%       'diagnostic_time'  time at which to evaluate the pressure waveform
%       'Rddot_override'   acceleration used only in the modal equation
%       'rad_mod_override' anisotropic radial source to use
%       'ep_mod_override'  anisotropic modal source vector to use
%
%   This helper is meant for static phase portraits and force attribution.

    solver_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(solver_dir,'..','common'));

    [diagnostic_time, Rddot_override, rad_mod_override, ...
        ep_mod_override, solver_args] = ...
        strip_diagnostic_args(varargin);

    [~, eqns_opts, solve_opts, init_opts, ~, ~, ~, acos_opts, wave_opts, ...
        sigma_opts, ~, ~, pert_opts] = evalc('f_call_params(solver_args{:});');

    radial = eqns_opts(1);
    bubtherm = eqns_opts(2);
    medtherm = eqns_opts(3);
    stress = eqns_opts(4);
    masstrans = eqns_opts(6);
    perturbed = eqns_opts(7);
    pertmod = eqns_opts(8);

    if ~perturbed
        error('f_imr_rhs_diagnostics requires perturbed = 1.');
    end
    if bubtherm || medtherm || masstrans
        error(['f_imr_rhs_diagnostics currently supports the polytropic ' ...
            'gas branch only (bubtherm = medtherm = masstrans = 0).']);
    end

    Nv = solve_opts(6);

    Pb_star = init_opts(3);
    Pv_star = init_opts(6);
    Req = init_opts(7);

    Cstar = acos_opts(1);
    GAMa = acos_opts(2);
    kappa = acos_opts(3);
    nstate = acos_opts(4);
    hugoniot_s = acos_opts(5);
    sam = 1 + GAMa;
    no = (nstate - 1)/nstate;
    nog = (nstate - 1)/2;

    We = sigma_opts(1);
    Re8 = sigma_opts(2);
    v_a = sigma_opts(4);
    v_nc = sigma_opts(5);
    Ca = sigma_opts(6);
    alphax = sigma_opts(7);
    LAM = sigma_opts(8);
    De = sigma_opts(9);
    JdotA = sigma_opts(10);
    nu_model = sigma_opts(11);
    zeNO = sigma_opts(13);
    iDRe = sigma_opts(14);
    graded = sigma_opts(15);
    Ca1 = sigma_opts(16);
    l1 = sigma_opts(17);
    l2 = sigma_opts(18);
    ani1 = sigma_opts(19);
    ani2 = sigma_opts(20);
    iWe = 1/We;

    n = pert_opts.n(:);
    m = pert_opts.m(:);
    if isempty(m)
        m = zeros(size(n));
    end
    epnmeq = pert_opts.epnmeq(:);
    np = length(n);
    min_q_len = 3 + 2*np;

    % Stress-state index layout, matching f_imr_fd for the polytropic branch.
    ivisco1 = 4:(3+Nv);
    ivisco2 = (4+Nv):(3+2*Nv);
    x_len = max([3, ivisco1, ivisco2]);

    wave_type = wave_opts(6);
    if wave_type < 0
        wave_poly = wave_opts(7);
        wave_dpoly = wave_opts(8);
    else
        wave_poly = [];
        wave_dpoly = [];
    end
    pvarargin = [wave_opts(1), wave_opts(2), wave_opts(3), wave_opts(4), ...
        wave_opts(5), wave_type, wave_poly, wave_dpoly];
    [Pf8, Pf8dot] = f_pinfinity(diagnostic_time, pvarargin);

    if sum(abs([ani1 ani2])) > 0
        [chiS, M1, M2, M3, M4, M5] = f_ani_ortho([0; n], [0; m]);
    else
        chiS = [];
        M1 = [];
        M2 = [];
        M3 = [];
        M4 = [];
        M5 = [];
    end

    diag = evaluate_states(q);

    function out = evaluate_states(q_in)
        if isvector(q_in)
            out = evaluate_one(q_in(:));
            return
        end

        if size(q_in,1) < min_q_len && size(q_in,2) >= min_q_len
            q_in = q_in.';
        end
        if size(q_in,1) < min_q_len
            error('q must contain [R; Rdot; Pb_state; epnm(:); epnmd(:)].');
        end

        nstates = size(q_in,2);
        first = evaluate_one(q_in(:,1));

        scalar_fields = {'Rddot','modal_Rddot','rad_mod','P','Pdot', ...
            'Pb_state','Pf8','Pf8dot','S','Sdot'};
        vector_fields = {'epddot','ep_mod','xi','eta','elastns', ...
            'sselastns','viscns','epinertians'};
        modal_fields = fieldnames(first.modal_terms);

        out = struct();
        for sf = 1:length(scalar_fields)
            out.(scalar_fields{sf}) = zeros(1,nstates);
        end
        for vf = 1:length(vector_fields)
            out.(vector_fields{vf}) = zeros(length(first.(vector_fields{vf})),nstates);
        end
        out.Z1dot = cell(1,nstates);
        out.Z2dot = cell(1,nstates);
        out.modal_terms = struct();
        for mf = 1:length(modal_fields)
            value = first.modal_terms.(modal_fields{mf});
            out.modal_terms.(modal_fields{mf}) = zeros(numel(value),nstates);
        end

        assign_state(1, first);
        for state_idx = 2:nstates
            assign_state(state_idx, evaluate_one(q_in(:,state_idx)));
        end

        function assign_state(state_idx, one)
            for sf_idx = 1:length(scalar_fields)
                field = scalar_fields{sf_idx};
                out.(field)(state_idx) = one.(field);
            end
            for vf_idx = 1:length(vector_fields)
                field = vector_fields{vf_idx};
                out.(field)(:,state_idx) = one.(field)(:);
            end
            out.Z1dot{state_idx} = one.Z1dot;
            out.Z2dot{state_idx} = one.Z2dot;
            for mf_idx = 1:length(modal_fields)
                field = modal_fields{mf_idx};
                out.modal_terms.(field)(:,state_idx) = one.modal_terms.(field)(:);
            end
        end
    end

    function one = evaluate_one(q_state)
        if length(q_state) < min_q_len
            error('q must contain [R; Rdot; Pb_state; epnm(:); epnmd(:)].');
        end

        R = q_state(1);
        Rdot = q_state(2);
        Pb_state = q_state(3);
        if isnan(Pb_state)
            Pb_state = Pb_star;
        end
        epnm = q_state(4:(3+np));
        epnmd = q_state((4+np):(3+2*np));
        stress_state = q_state((4+2*np):end);

        P = (Pb_state - Pv_star)*(1/R)^(3*kappa) + Pv_star;
        Pdot = -3*kappa*P*Rdot/R;

        if nu_model ~= 0
            [~, intfnu, dintfnu, ddintfnu] = f_viscosity(nu_model, Rdot, ...
                R, v_a, v_nc, 0);
        else
            intfnu = 0;
            dintfnu = 0;
            ddintfnu = 0;
        end

        X = zeros(x_len,1);
        X(1:3) = [R; Rdot; Pb_state];
        if ~isempty(stress_state)
            nstress = min(length(stress_state), length(ivisco1) + length(ivisco2));
            stress_state = stress_state(1:nstress);
            if nstress > 0
                stress_idx = [ivisco1 ivisco2];
                X(stress_idx(1:nstress)) = stress_state;
            end
        end

        if graded
            [S, Sdot, Z1dot, Z2dot] = f_stress_graded(radial, stress, Req, ...
                R, Ca, Ca1, Re8, Rdot, alphax, intfnu, dintfnu, iDRe, ...
                l1, l2, v_a, v_nc);
        else
            [S, Sdot, Z1dot, Z2dot] = f_stress(stress, X, Req, R, Ca, ...
                De, Re8, Rdot, alphax, ivisco1, ivisco2, LAM, zeNO, [], ...
                intfnu, dintfnu, iDRe);
        end

        if sum(abs([ani1 ani2])) > 0
            [Ts1, Ts2, Ts3, T1, T2, T3, T4, T5] = ...
                f_ani_ortho_time_coeffs(n, m, R/Req, Req, Ca, ani1, ...
                ani2, epnmeq, epnm);
            ortho_vect = chiS(:,1)*Ts1 + chiS(:,2)*Ts2 + ...
                chiS(:,3)*Ts3 + M1*T1 + M2*T2 + M3*T3 + M4*T4 + M5*T5;
            rad_mod = ortho_vect(1);
            ep_mod = ortho_vect(2:end);
        else
            rad_mod = 0;
            ep_mod = zeros(np,1);
        end
        if ~isnan(rad_mod_override)
            rad_mod = rad_mod_override;
        end
        if ~isempty(ep_mod_override)
            if isscalar(ep_mod_override)
                ep_mod = ep_mod_override.*ones(np,1);
            else
                ep_mod = ep_mod_override(:);
            end
        end

        Rddot = f_radial_eq(radial, P, Pdot, Pf8, Pf8dot, iWe, R, Rdot, ...
            S, Sdot, Cstar, sam, no, GAMa, nstate, nog, hugoniot_s, ...
            JdotA, ddintfnu, iDRe, rad_mod);
        if isnan(Rddot_override)
            modal_Rddot = Rddot;
        else
            modal_Rddot = Rddot_override;
        end

        [epddot, modal_terms] = f_compute_perturb_coeffs(epnm, epnmd, ...
            epnmeq, ep_mod, R, Rdot, modal_Rddot, n.', Req, We, Re8, Ca, ...
            alphax, pertmod);

        one = struct();
        one.Rddot = Rddot;
        one.modal_Rddot = modal_Rddot;
        one.epddot = epddot;
        one.rad_mod = rad_mod;
        one.ep_mod = ep_mod;
        one.P = P;
        one.Pdot = Pdot;
        one.Pb_state = Pb_state;
        one.Pf8 = Pf8;
        one.Pf8dot = Pf8dot;
        one.S = S;
        one.Sdot = Sdot;
        one.Z1dot = Z1dot;
        one.Z2dot = Z2dot;
        one.xi = modal_terms.xi;
        one.eta = modal_terms.eta;
        one.elastns = modal_terms.elastns;
        one.sselastns = modal_terms.sselastns;
        one.viscns = modal_terms.viscns;
        one.epinertians = modal_terms.epinertians;
        one.modal_terms = modal_terms;
    end
end

function [diagnostic_time, Rddot_override, rad_mod_override, ...
    ep_mod_override, solver_args] = strip_diagnostic_args(args)
    if mod(length(args),2) == 1
        error('Inputs after q must be name/value pairs.');
    end

    diagnostic_time = 0;
    Rddot_override = NaN;
    rad_mod_override = NaN;
    ep_mod_override = [];
    keep = true(size(args));
    for k = 1:2:length(args)
        name = lower(args{k});
        if strcmp(name, 'diagnostic_time') || strcmp(name, 'time') || strcmp(name, 't')
            diagnostic_time = args{k+1};
            keep(k:k+1) = false;
        elseif strcmp(name, 'rddot_override') || strcmp(name, 'modal_rddot')
            Rddot_override = args{k+1};
            keep(k:k+1) = false;
        elseif strcmp(name, 'rad_mod_override')
            rad_mod_override = args{k+1};
            keep(k:k+1) = false;
        elseif strcmp(name, 'ep_mod_override')
            ep_mod_override = args{k+1};
            keep(k:k+1) = false;
        end
    end
    solver_args = args(keep);
end
