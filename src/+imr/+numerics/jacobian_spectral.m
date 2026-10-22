% file jacobian_spectral.m
% brief contains function imr.numerics.jacobian_spectral

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function computes the analytic Jacobian of the spectral
% right-hand side (imr.solve_spectral). The thermal fields, the thermal
% boundary condition solve, the bubble pressure and the spectral stress
% fields are differentiated exactly with the value/derivative pairs of
% imr.numerics.dual_operations. The boundary solve A(s)*q = rhs, with A(s) = A0 - (s-1)*e1*v'
% factored once, is differentiated as dq = A(s)\(d(rhs) + e1*(v'*q)*ds).
% The low-dimensional closure (radial model, closed-form stress models,
% viscosity, medium dissipation) is differentiated by central differences
% of the closure alone. The struct jp holds the indices, operators,
% constants and closure handles set up by imr.solve_spectral.
function J = jacobian_spectral(t, X, jp)

    n = numel(X);
    d = imr.numerics.dual_operations(n);
    iz = jp.iz;

    R = d.var(X, 1);
    Rdot = d.var(X, 2);
    R2 = d.mul(R, R);
    J = sparse(n, n);
    J(1,:) = Rdot.J;

    if jp.bubtherm
        P = d.var(X, 3);
        aC = d.var(X, jp.ia);
        Ne = numel(jp.ia);
        % auxiliary temperature and derivatives at the collocation points
        th = d.lin(jp.gA, aC);
        dth = d.lin(jp.gDA, aC);
        ddth = d.lin(jp.gDDA, aC);
        [Tv,sq] = imr.physics.thermal_inverse(th.v,jp.alpha_g,jp.beta_g);
        T = d.fun(th,Tv,1./sq);
        D = d.scale(d.div(d.add(d.scale(d.mul(T, T), jp.alpha_g), ...
            d.scale(T, jp.beta_g)), P), jp.kapover);
        dth1 = d.el(dth, 1);
        % bubble pressure
        inner = d.sub(d.scale(d.div(dth1, R), (jp.kappa - 1)*jp.chi), ...
            d.scale(d.mul(P, Rdot), jp.kappa));
        Pdot = d.scale(d.div(inner, R), 3);
        % velocity field and auxiliary temperature equation
        term1 = d.scale(D, 2./jp.y);
        term2 = d.scale(d.div(d.sub(dth, d.scale(dth1, jp.y)), P), jp.kapover);
        U = d.scale(d.div(d.sub(term1, term2), R2), jp.chi);
        thdot = d.add(d.add(d.mul(Pdot, D), d.mul(U, dth)), ...
            d.scale(d.div(d.mul(D, ddth), R2), jp.chi));
        % bubble center (y = 0) equation
        DC = d.el(D, Ne);
        dthC = d.el(dth, Ne);
        snn = d.lin(jp.nn', aC);
        inner2 = d.add(d.scale(d.mul(DC, snn), 8), ...
            d.scale(d.div(d.mul(dthC, dthC), P), jp.kapover));
        thC = d.add(d.sub(d.mul(Pdot, DC), d.scale(d.div(inner2, R2), jp.chi)), ...
            d.scale(d.div(d.mul(DC, d.el(ddth, Ne)), R2), jp.chi));
        thdot = d.set(thdot, Ne, thC);

        if jp.medtherm
            bC = d.var(X, jp.ib);
            dTm = d.lin(jp.mDA, bC);
            ddTm = d.lin(jp.mDDA, bC);
            g1 = (1./jp.yT2 - jp.yT)/2;
            g2 = (1 + jp.xi)/(2*jp.Lt) - 1./jp.yT;
            cA = (1 + jp.xi).^2/jp.Lt;
            br = d.add(d.scale(d.div(d.con(jp.Foh), R), g2), d.scale(Rdot, g1));
            adv = d.mul(d.mul(d.scale(d.div(d.con(1), R), cA), br), dTm);
            dif = d.mul(d.scale(d.div(d.con(1), R2), ...
                jp.Foh/4*(1 + jp.xi).^4/jp.Lt^2), ddTm);
            Tmdot = d.add(adv, dif);
            % stress dissipation: R and Rdot partials by central
            % differences at a fixed stress field, and the part linear in
            % the polymer stress field, -2*Br*a*sp with a = (Rdot/R)/yT^3,
            % exactly
            hR = eps^(1/3)*max(abs(R.v),1);
            hU = eps^(1/3)*max(abs(Rdot.v),1);
            tg_R = (jp.taugradu(R.v+hR,Rdot.v,X) - ...
                jp.taugradu(R.v-hR,Rdot.v,X))/(2*hR);
            tg_U = (jp.taugradu(R.v,Rdot.v+hU,X) - ...
                jp.taugradu(R.v,Rdot.v-hU,X))/(2*hU);
            tg = jp.taugradu(R.v,Rdot.v,X);
            Jtg = sparse(tg_R(:))*R.J + sparse(tg_U(:))*Rdot.J;
            if jp.spectral_stress || jp.zener_field
                dsp = -2*jp.Br*(Rdot.v/R.v)./jp.yT3.*jp.ZZT;
            end
            if jp.spectral_stress
                Nv = numel(iz)/2;
                Jtg(:,iz(1:Nv)) = Jtg(:,iz(1:Nv)) + dsp;
                Jtg(:,iz(Nv+1:end)) = Jtg(:,iz(Nv+1:end)) - dsp;
            elseif jp.zener_field
                Jtg(:,jp.izf) = Jtg(:,jp.izf) + dsp;
            end
            Tmdot = d.add(Tmdot, struct('v', tg, 'J', Jtg));
            Tmdot = d.set(Tmdot, numel(Tmdot.v), d.con(0));
            % boundary condition solve with the Sherman-Morrison update
            s = d.add(d.scale(d.el(T, 1), jp.alpha_g), d.con(jp.beta_g));
            rv = vertcat(0, thdot.v(2:end), Tmdot.v(2:end), 0);
            rJ = vertcat(sparse(1,n), thdot.J(2:end,:), Tmdot.J(2:end,:), ...
                sparse(1,n));
            q = f_ainv(rv, s.v);
            e1 = sparse(1, 1, 1, numel(rv), 1);
            qJ = f_ainv(full(rJ + e1*((jp.vrow'*q)*s.J)), s.v);
            J(jp.ia,:) = qJ(1:Ne,:);
            J(jp.ib,:) = qJ(Ne+1:end,:);
        else
            % cold-liquid approximation
            J(jp.ia,:) = jp.gAI*vertcat(sparse(1,n), thdot.J(2:end,:));
        end
    else
        % polytropic gas and vapor
        dP = jp.Pb_star - jp.Pv_star;
        P = d.fun(R, dP*R.v^(-3*jp.kappa) + jp.Pv_star, ...
            -3*jp.kappa*dP*R.v^(-3*jp.kappa-1));
        Pg = d.sub(P, d.con(jp.Pv_star));
        Pdot = d.scale(d.div(d.mul(Pg, Rdot), R), -3*jp.kappa);
    end
    J(3,:) = Pdot.J;

    % stress and radial closure
    if jp.spectral_stress
        Nv = numel(iz)/2;
        cC = d.var(X, iz(1:Nv));
        dC = d.var(X, iz(Nv+1:end));
        Rr = d.div(Rdot, R);
        av = d.scale(Rr, 1./jp.yV.^3);
        dzdt = d.scale(Rr, jp.dzdy.*(1./jp.yV.^2 - jp.yV));
        trr = d.lin(jp.sCA, cC);
        ttt = d.lin(jp.sCA, dC);
        trrdot = d.scale(d.add(d.add(d.mul(dzdt, d.lin(jp.sCAd, cC)), ...
            d.scale(trr, 1/jp.De)), d.add(d.scale(d.mul(av, trr), 4), ...
            d.scale(av, 4*jp.etap/jp.De))), -1);
        tttdot = d.add(d.add(d.scale(d.mul(dzdt, d.lin(jp.sCAd, dC)), -1), ...
            d.scale(ttt, -1/jp.De)), d.add(d.scale(d.mul(av, ttt), 2), ...
            d.scale(av, 2*jp.etap/jp.De)));
        cdot = d.lin(jp.sCI, trrdot);
        ddot = d.lin(jp.sCI, tttdot);
        S = d.sub(d.lin(2*jp.cdd', d.sub(cC, dC)), ...
            d.scale(Rr, 4*jp.LAM/jp.Re8));
        Sdot = d.add(d.lin(2*jp.cdd', d.sub(cdot, ddot)), ...
            d.scale(d.mul(Rr, Rr), 4*jp.LAM/jp.Re8));
        J(iz(1:Nv),:) = cdot.J;
        J(iz(Nv+1:end),:) = ddot.J;
        v0 = vertcat(R.v, Rdot.v, P.v, Pdot.v, S.v, Sdot.v);
        G = f_closure_fd(@(v) jp.closure_s(t,v), v0);
        J(2,:) = sparse(G)*vertcat(R.J, Rdot.J, P.J, Pdot.J, S.J, Sdot.J);
    else
        v0 = vertcat(R.v, Rdot.v, P.v, Pdot.v, X(iz));
        G = f_closure_fd(@(v) jp.closure(t,v,X), v0);
        JZ = sparse(1:numel(iz), iz, 1, numel(iz), n);
        J([2 iz],:) = sparse(G)*vertcat(R.J, Rdot.J, P.J, Pdot.J, JZ);
    end

    % Zener field rows: linear in the coefficients e, and R and Rdot
    % partials by central differences
    if jp.zener_field
        Rr = Rdot.v/R.v;
        dzdt = jp.dzdy.*Rr.*(1./jp.yV.^2 - jp.yV);
        Ae = jp.sCI*(-dzdt.*jp.sCAd - jp.sCA/jp.De - ...
            (3*Rr*(1 - jp.yV.^-3)).*jp.sCA);
        hR = eps^(1/3)*max(abs(R.v),1);
        hU = eps^(1/3)*max(abs(Rdot.v),1);
        e_R = (jp.zener_rhs(R.v+hR,Rdot.v,X) - ...
            jp.zener_rhs(R.v-hR,Rdot.v,X))/(2*hR);
        e_U = (jp.zener_rhs(R.v,Rdot.v+hU,X) - ...
            jp.zener_rhs(R.v,Rdot.v-hU,X))/(2*hU);
        Jf = sparse(e_R(:))*R.J + sparse(e_U(:))*Rdot.J;
        Jf(:,jp.izf) = Jf(:,jp.izf) + Ae;
        J(jp.izf,:) = Jf;
    end

    % A(s)\M for the thermal boundary condition system
    function Y = f_ainv(M, sv)
        Y0 = jp.QA0*(jp.UA0\(jp.LA0\(jp.PA0*M)));
        Y = Y0 + jp.w1*(((sv-1)*(jp.vrow'*Y0))/(1 - (sv-1)*jp.vw1));
    end

end

% central differences of a low-dimensional closure g(v)
function G = f_closure_fd(g, v0)
    nv = numel(v0);
    g0 = g(v0);
    G = zeros(numel(g0), nv);
    for k = 1:nv
        hk = eps^(1/3)*max(abs(v0(k)),1);
        vp = v0;
        vm = v0;
        vp(k) = vp(k) + hk;
        vm(k) = vm(k) - hk;
        G(:,k) = (g(vp) - g(vm))/(2*hk);
    end
end
