% file f_jacobian_fd.m
% brief contains function f_jacobian_fd

% brief This function computes the analytic Jacobian of the finite
% difference right-hand side (f_imr_fd) for bubtherm = 1, with or without
% medium heat transfer and vapor mass transfer. The field equations are
% differentiated exactly: every intermediate quantity is carried as a
% pair (value v, sparse derivative J = dv/dX), so the chain rule through
% the wall boundary conditions, the bubble pressure and the transport
% equations is applied term by term. The bubble-wall temperature enters
% through the implicit function theorem on the wall balance. Only the
% low-dimensional closure (radial model, stress model, viscosity, medium
% dissipation) is differentiated by central differences of the closure
% alone. The struct jp holds the state indices, operators, constants and
% the closure handles set up by f_imr_fd.
function J = f_jacobian_fd(t, X, jp)
    
    n = numel(X);
    ith = jp.ith;
    itm = jp.itm;
    ikv = jp.ikv;
    iz = jp.iz;
    nth = numel(ith);
    
    R = f_var(1);
    Rdot = f_var(2);
    P = f_var(3);
    
    % bubble temperature, with the wall value set by the boundary condition
    th = f_vars(ith);
    if jp.masstrans
        kvs = f_vars(ikv);
    end
    if jp.medtherm
        Tm = f_vars(itm);
        if jp.masstrans
            theta_bw = jp.theta_bw(X);
            wc = jp.wall;
            [~,dr] = f_wall_mass_bc([theta_bw, Tm.v(2), Tm.v(3), th.v(end-1), ...
                th.v(end-2), kvs.v(end), kvs.v(end-1), kvs.v(end-2), P.v], wc);
            cols = [itm(2) itm(3) ith(end-1) ith(end-2) ikv(end) ...
                ikv(end-1) ikv(end-2) 3];
            row = sparse(1, cols, -dr(2:9)/dr(1), 1, n);
        else
            gTm = jp.wall.grad_Tm;
            gTr = jp.wall.grad_Tr;
            K = gTm(2:3)*Tm.v(2:3) + gTr(2:3)*th.v(end-1:-1:end-2);
            [theta_bw, dK] = f_wall_thermal_bc(K, gTm(1), gTr(1), jp.alpha_g);
            row = sparse(1, [itm(2) itm(3) ith(end-1) ith(end-2)], ...
                dK*[gTm(2) gTm(3) gTr(2) gTr(3)], 1, n);
        end
        th = f_set(th, nth, struct('v',theta_bw,'J',row));
    end
    
    % bubble temperature from theta and the mixture coefficient
    if jp.masstrans
        am = f_add(f_con(jp.alpha_g), f_scale(kvs, jp.alpha_v - jp.alpha_g));
        s = sqrt(1 + 2*am.v.*th.v);
        Tv = (am.v - 1 + s)./am.v;
        T = struct('v', Tv, 'J', f_diag(1./s)*th.J + ...
            f_diag((1 + th.v./s - Tv)./am.v)*am.J);
    else
        am = jp.kv0*jp.alpha_v + (1 - jp.kv0)*jp.alpha_g;
        s = sqrt(1 + 2*am*th.v);
        T = struct('v', (am - 1 + s)/am, 'J', f_diag(1./s)*th.J);
    end
    Te = f_el(T, nth);
    
    % vapor concentration, wall value saturated at (T(end), P)
    if jp.masstrans
        Pv = f_fun(Te, jp.C1_pv*exp(jp.C2_pv/Te.v), ...
            -jp.C1_pv*exp(jp.C2_pv/Te.v)*jp.C2_pv/Te.v^2);
        dem = f_add(f_con(1 - jp.Rvg), f_scale(f_div(P, Pv), jp.Rvg));
        kv = f_set(kvs, nth, f_div(f_con(1), dem));
        dkv = f_lin(jp.D1, kv);
        ddkv = f_lin(jp.D2, kv);
        Rmix = f_add(f_con(jp.Rg), f_scale(kv, jp.Rva));
        RDkv = f_scale(f_div(dkv, Rmix), jp.Rva);
    end
    
    % bubble temperature gradients
    dth = f_lin(jp.D1, th);
    ddth = f_lin(jp.D2, th);
    R2 = f_mul(R, R);
    
    % bubble pressure
    inner = f_sub(f_scale(f_div(f_el(dth, nth), R), jp.chi*(jp.kappa - 1)), ...
        f_scale(f_mul(P, Rdot), jp.kappa));
    if jp.masstrans
        den = f_mul(f_mul(Te, R), f_mul(f_el(Rmix, nth), ...
            f_add(f_con(1), f_scale(f_el(kv, nth), -1))));
        inner = f_add(inner, f_scale(f_div(f_mul(P, f_el(dkv, nth)), den), ...
            jp.kappa*jp.Fom*jp.Rv));
    end
    Pdot = f_scale(f_div(inner, R), 3);
    
    % bubble velocity field
    Uvel = f_div(f_sub(f_scale(f_div(dth, R), jp.chi*(jp.kappa - 1)), ...
        f_scale(f_mul(R, Pdot), jp.y/3)), f_scale(P, jp.kappa));
    if jp.masstrans
        Uvel = f_add(Uvel, f_scale(f_div(RDkv, R), jp.Fom));
    end
    av = f_sub(Uvel, f_scale(Rdot, jp.y));
    
    % thermal conductivity
    if jp.masstrans
        Kg = f_add(f_scale(T, jp.alpha_g), f_con(jp.beta_g));
        Kv = f_add(f_scale(T, jp.alpha_v), f_con(jp.beta_v));
        Kst = f_add(Kg, f_mul(kv, f_sub(Kv, Kg)));
    else
        Ka = jp.kv0*jp.alpha_v + (1 - jp.kv0)*jp.alpha_g;
        Kb = jp.kv0*jp.beta_v + (1 - jp.kv0)*jp.beta_g;
        Kst = f_add(f_scale(T, Ka), f_con(Kb));
    end
    
    % bubble temperature equation
    adv = f_scale(f_div(f_mul(dth, av), R), -1);
    W = f_scale(f_div(f_mul(Kst, T), P), jp.kapover);
    dfc = f_add(f_scale(f_div(ddth, R2), jp.chi), Pdot);
    thdot = f_add(adv, f_mul(dfc, W));
    if jp.masstrans
        thdot = f_add(thdot, f_scale(f_div(f_mul(dkv, dth), ...
            f_mul(Rmix, R2)), jp.Fom*jp.Rva));
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
        amix = f_add(f_con(jp.alpha_g), f_scale(kv, jp.alpha_v - jp.alpha_g));
        q = f_add(f_con(1), f_scale(f_mul(amix, th), 2));
        sq = f_fun(q, sqrt(q.v), 0.5./sqrt(q.v));
        nld = f_mul(dkv, f_add(f_div(dth, f_mul(sq, T)), RDkv));
        kvdot = f_sub(f_scale(f_div(f_sub(ddkv, nld), R2), jp.Fom), ...
            f_div(f_mul(av, dkv), R));
        Jkv = kvdot.J;
        % wall value held by the saturation condition
        Jkv(end,:) = 0;
        J(ikv,:) = Jkv;
    end
    
    % material temperature equation
    if jp.medtherm
        Tme = f_set(Tm, 1, Te);
        dTm = f_lin(jp.Dm1, Tme);
        ddTm = f_lin(jp.Dm2, Tme);
        g1 = (1 - jp.yT3)./(2*jp.yT2);
        g2 = (jp.xi + 1)/(2*jp.Lt) - 1./jp.yT;
        cA = (1 + jp.xi).^2/jp.Lt;
        Bm = jp.Foh*(jp.xi + 1).^4/(4*jp.Lt^2);
        bracket = f_add(f_scale(Rdot, g1), f_scale(f_div(f_con(jp.Foh), R), g2));
        advm = f_mul(f_mul(f_scale(f_div(f_con(1), R), cA), bracket), dTm);
        diffm = f_mul(f_scale(f_div(f_con(1), R2), Bm), ddTm);
        Tmdot = f_add(advm, diffm);
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
    
    % value/derivative pair helpers; J has one column per state
    
    function a = f_var(k)
        a = struct('v', X(k), 'J', sparse(1, k, 1, 1, n));
    end
    
    function a = f_vars(idx)
        a = struct('v', X(idx), 'J', sparse(1:numel(idx), idx, 1, numel(idx), n));
    end
    
    function a = f_con(v)
        a = struct('v', v, 'J', sparse(numel(v), n));
    end
    
    function a = f_el(b, k)
        a = struct('v', b.v(k), 'J', b.J(k,:));
    end
    
    function a = f_set(a, k, b)
        a.v(k) = b.v;
        a.J(k,:) = b.J;
    end
    
    function a = f_lin(M, b)
        a = struct('v', M*b.v, 'J', M*b.J);
    end
    
    % elementwise function with derivative dval
    function a = f_fun(b, val, dval)
        a = struct('v', val, 'J', f_diag(dval)*b.J);
    end
    
    % multiplication by a constant scalar or vector s
    function a = f_scale(b, s)
        if numel(s) == numel(b.v)
            a = struct('v', b.v.*s, 'J', f_diag(s)*b.J);
        elseif isscalar(b.v)
            a = struct('v', b.v*s, 'J', sparse(s(:))*b.J);
        else
            a = struct('v', b.v*s, 'J', s*b.J);
        end
    end
    
    function a = f_add(b, c)
        m = max(numel(b.v), numel(c.v));
        a = struct('v', b.v + c.v, 'J', f_rows(b.J, m) + f_rows(c.J, m));
    end
    
    function a = f_sub(b, c)
        a = f_add(b, f_scale(c, -1));
    end
    
    function a = f_mul(b, c)
        if numel(b.v) == numel(c.v)
            Ja = f_diag(c.v)*b.J + f_diag(b.v)*c.J;
        elseif isscalar(b.v)
            Ja = sparse(c.v(:))*b.J + b.v*c.J;
        else
            Ja = sparse(b.v(:))*c.J + c.v*b.J;
        end
        a = struct('v', b.v.*c.v, 'J', Ja);
    end
    
    function a = f_div(b, c)
        a = f_mul(b, f_fun(c, 1./c.v, -1./c.v.^2));
    end
    
    % broadcast a scalar's derivative row to m rows
    function Jr = f_rows(Jb, m)
        if size(Jb,1) == m
            Jr = Jb;
        else
            Jr = repmat(Jb, m, 1);
        end
    end
    
    function Dg = f_diag(d)
        Dg = spdiags(d(:), 0, numel(d), numel(d));
    end
    
end
