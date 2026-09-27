% file f_jacobian_nodal.m
% brief contains function f_jacobian_nodal

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function computes the analytic Jacobian of the nodal
% right-hand side (f_imr_nodal) for bubtherm = 1, with or without
% medium heat transfer and vapor mass transfer. The field equations are
% differentiated exactly: every intermediate quantity is carried as a
% pair (value v, sparse derivative J = dv/dX), so the chain rule through
% the wall boundary conditions, the bubble pressure and the transport
% equations is applied term by term. The bubble-wall temperature enters
% through the implicit function theorem on the wall balance. Only the
% low-dimensional closure (radial model, stress model, viscosity, medium
% dissipation) is differentiated by central differences of the closure
% alone. The struct jp holds the state indices, operators, constants and
% the closure handles set up by f_imr_nodal, for either differentiation
% scheme.
function J = f_jacobian_nodal(t, X, jp)
    
    n = numel(X);
    ith = jp.ith;
    itm = jp.itm;
    ikv = jp.ikv;
    iz = jp.iz;
    nth = numel(ith);
    d = f_dual_ops(n);
    
    R = d.var(X, 1);
    Rdot = d.var(X, 2);
    P = d.var(X, 3);
    
    % bubble temperature, with the wall value set by the boundary condition
    th = d.var(X, ith);
    if jp.masstrans
        kvs = d.var(X, ikv);
    end
    if jp.medtherm
        Tm = d.var(X, itm);
        if jp.masstrans
            theta_bw = jp.theta_bw(t,X);
            wc = jp.wall;
            im = wc.iwall_m;
            ib = wc.iwall_b;
            [~,dr] = f_wall_mass_bc([theta_bw, Tm.v(im)', th.v(ib)', ...
                kvs.v(end), kvs.v(ib)', P.v], wc);
            cols = [itm(im) ith(ib) ikv(end) ikv(ib) 3];
            row = sparse(1, cols, -dr(2:end)/dr(1), 1, n);
        else
            gTm = jp.wall.grad_Tm;
            gTr = jp.wall.grad_Tr;
            % same root as the solver; for c2*u^2 + 2*c1*u + q(K) = 0 the
            % implicit derivative is d(theta_w)/dK = -u/(c2*u + c1)
            theta_bw = jp.theta_bw(t,X);
            u = sqrt((jp.alpha_g + jp.beta_g)^2 + 2*jp.alpha_g*theta_bw);
            dK = -u/(gTr(1)*u + gTm(1));
            row = sparse(1, [itm(jp.wall.iwall_m) ith(jp.wall.iwall_b)], ...
                dK*[gTm(2:end) gTr(2:end)], 1, n);
        end
        th = d.set(th, nth, struct('v',theta_bw,'J',row));
    end
    
    % bubble temperature from theta and the mixture coefficients,
    % T = (s - b)/a with s = sqrt((a + b)^2 + 2*a*theta), the conductivity
    if jp.masstrans
        am = d.add(d.con(jp.alpha_g), d.scale(kvs, jp.alpha_v - jp.alpha_g));
        bm = d.add(d.con(jp.beta_g), d.scale(kvs, jp.beta_v - jp.beta_g));
        ab = am.v + bm.v;
        s = sqrt(ab.^2 + 2*am.v.*th.v);
        Tv = (s - bm.v)./am.v;
        T = struct('v', Tv, 'J', d.diag(1./s)*th.J + ...
            d.diag(((ab + th.v)./s - Tv)./am.v)*am.J + ...
            d.diag((ab./s - 1)./am.v)*bm.J);
    else
        am = jp.kv0*jp.alpha_v + (1 - jp.kv0)*jp.alpha_g;
        bm = jp.kv0*jp.beta_v + (1 - jp.kv0)*jp.beta_g;
        s = sqrt((am + bm)^2 + 2*am*th.v);
        T = struct('v', (s - bm)/am, 'J', d.diag(1./s)*th.J);
    end
    Te = d.el(T, nth);
    
    % vapor concentration, wall value saturated at (T(end), P)
    if jp.masstrans
        Pv = d.fun(Te, jp.C1_pv*exp(jp.C2_pv/Te.v), ...
            -jp.C1_pv*exp(jp.C2_pv/Te.v)*jp.C2_pv/Te.v^2);
        dem = d.add(d.con(1 - jp.Rvg), d.scale(d.div(P, Pv), jp.Rvg));
        kv = d.set(kvs, nth, d.div(d.con(1), dem));
        dkv = d.lin(jp.D1, kv);
        ddkv = d.lin(jp.D2, kv);
        Rmix = d.add(d.con(jp.Rg), d.scale(kv, jp.Rva));
        RDkv = d.scale(d.div(dkv, Rmix), jp.Rva);
    end
    
    % bubble temperature gradients
    dth = d.lin(jp.D1, th);
    ddth = d.lin(jp.D2, th);
    R2 = d.mul(R, R);
    
    % bubble pressure
    inner = d.sub(d.scale(d.div(d.el(dth, nth), R), jp.chi*(jp.kappa - 1)), ...
        d.scale(d.mul(P, Rdot), jp.kappa));
    if jp.masstrans
        den = d.mul(d.mul(Te, R), d.mul(d.el(Rmix, nth), ...
            d.add(d.con(1), d.scale(d.el(kv, nth), -1))));
        inner = d.add(inner, d.scale(d.div(d.mul(P, d.el(dkv, nth)), den), ...
            jp.kappa*jp.Fom*jp.Rv));
    end
    Pdot = d.scale(d.div(inner, R), 3);
    
    % bubble velocity field
    Uvel = d.div(d.sub(d.scale(d.div(dth, R), jp.chi*(jp.kappa - 1)), ...
        d.scale(d.mul(R, Pdot), jp.y/3)), d.scale(P, jp.kappa));
    if jp.masstrans
        Uvel = d.add(Uvel, d.scale(d.div(RDkv, R), jp.Fom));
    end
    av = d.sub(Uvel, d.scale(Rdot, jp.y));
    
    % thermal conductivity
    if jp.masstrans
        Kg = d.add(d.scale(T, jp.alpha_g), d.con(jp.beta_g));
        Kv = d.add(d.scale(T, jp.alpha_v), d.con(jp.beta_v));
        Kst = d.add(Kg, d.mul(kv, d.sub(Kv, Kg)));
    else
        Ka = jp.kv0*jp.alpha_v + (1 - jp.kv0)*jp.alpha_g;
        Kb = jp.kv0*jp.beta_v + (1 - jp.kv0)*jp.beta_g;
        Kst = d.add(d.scale(T, Ka), d.con(Kb));
    end
    
    % bubble temperature equation
    adv = d.scale(d.div(d.mul(dth, av), R), -1);
    W = d.scale(d.div(d.mul(Kst, T), P), jp.kapover);
    dfc = d.add(d.scale(d.div(ddth, R2), jp.chi), Pdot);
    thdot = d.add(adv, d.mul(dfc, W));
    if jp.masstrans
        thdot = d.add(thdot, d.scale(d.div(d.mul(dkv, dth), ...
            d.mul(Rmix, R2)), jp.Fom*jp.Rva));
    end
    
    J = sparse(n, n);
    J(1,:) = Rdot.J;
    J(3,:) = Pdot.J;
    Jth = thdot.J;
    % wall value held by the boundary condition
    Jth(end,:) = 0;
    J(ith,:) = Jth;
    
    % vapor concentration equation
    if jp.masstrans
        % the square root of the Kirchhoff inverse is the conductivity
        nld = d.mul(dkv, d.add(d.div(dth, d.mul(Kst, T)), RDkv));
        kvdot = d.sub(d.scale(d.div(d.sub(ddkv, nld), R2), jp.Fom), ...
            d.div(d.mul(av, dkv), R));
        Jkv = kvdot.J;
        % wall value held by the saturation condition
        Jkv(end,:) = 0;
        J(ikv,:) = Jkv;
    end
    
    % material temperature equation
    if jp.medtherm
        Tme = d.set(Tm, 1, Te);
        dTm = d.lin(jp.Dm1, Tme);
        ddTm = d.lin(jp.Dm2, Tme);
        g1 = (1 - jp.yT3)./(2*jp.yT2);
        g2 = (jp.xi + 1)/(2*jp.Lt) - 1./jp.yT;
        cA = (1 + jp.xi).^2/jp.Lt;
        Bm = jp.Foh*(jp.xi + 1).^4/(4*jp.Lt^2);
        bracket = d.add(d.scale(Rdot, g1), d.scale(d.div(d.con(jp.Foh), R), g2));
        advm = d.mul(d.mul(d.scale(d.div(d.con(1), R), cA), bracket), dTm);
        diffm = d.mul(d.scale(d.div(d.con(1), R2), Bm), ddTm);
        Tmdot = d.add(advm, diffm);
        % dissipation partials by central differences
        hR = eps^(1/3)*max(abs(R.v),1);
        hU = eps^(1/3)*max(abs(Rdot.v),1);
        tg_R = (jp.taugradu(R.v+hR,Rdot.v,X) - jp.taugradu(R.v-hR,Rdot.v,X))/(2*hR);
        tg_U = (jp.taugradu(R.v,Rdot.v+hU,X) - jp.taugradu(R.v,Rdot.v-hU,X))/(2*hU);
        Jm = Tmdot.J + sparse(tg_R(:))*R.J + sparse(tg_U(:))*Rdot.J;
        Jm([1 end],:) = 0;
        J(itm,:) = Jm;
    end
    
    % radial and stress closure rows g(R, Rdot, P, Pdot, Z), by central
    % differences of the closure, chained through Pdot
    v0 = vertcat(R.v, Rdot.v, P.v, Pdot.v, X(iz));
    nv = numel(v0);
    G = zeros(1 + numel(iz), nv);
    for k = 1:nv
        hk = eps^(1/3)*max(abs(v0(k)),1);
        vp = v0;
        vm = v0;
        vp(k) = vp(k) + hk;
        vm(k) = vm(k) - hk;
        G(:,k) = (jp.closure(t,vp,X) - jp.closure(t,vm,X))/(2*hk);
    end
    JZ = sparse(1:numel(iz), iz, 1, numel(iz), n);
    J([2 iz],:) = sparse(G)*vertcat(R.J, Rdot.J, P.J, Pdot.J, JZ);
    
end
