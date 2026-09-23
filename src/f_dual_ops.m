% file f_dual_ops.m
% brief contains function f_dual_ops

% brief This function returns the operations used by the analytic
% Jacobians (f_jacobian_fd, f_jacobian_spectral). A quantity is carried as
% a struct with its value v and its sparse derivative J = dv/dX with
% respect to the n-dimensional state X, so that composing the operations
% applies the chain rule exactly. Scalars and vectors may be mixed; a
% scalar is broadcast where needed.
function d = f_dual_ops(n)
    
    d = struct('var',@f_var, 'con',@f_con, 'el',@f_el, 'set',@f_set, ...
        'lin',@f_lin, 'fun',@f_fun, 'scale',@f_scale, 'add',@f_add, ...
        'sub',@f_sub, 'mul',@f_mul, 'div',@f_div, 'diag',@f_diag);
    
    % state entries X(idx)
    function a = f_var(X, idx)
        a = struct('v', X(idx), 'J', sparse(1:numel(idx), idx, 1, numel(idx), n));
    end
    
    % constant
    function a = f_con(v)
        a = struct('v', v, 'J', sparse(numel(v), n));
    end
    
    % element(s) k
    function a = f_el(b, k)
        a = struct('v', b.v(k), 'J', b.J(k,:));
    end
    
    % replace element k of a by b
    function a = f_set(a, k, b)
        a.v(k) = b.v;
        a.J(k,:) = b.J;
    end
    
    % linear operator M*b
    function a = f_lin(M, b)
        a = struct('v', M*b.v, 'J', M*b.J);
    end
    
    % elementwise function with value val and derivative dval
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
    
    function Dg = f_diag(v)
        Dg = spdiags(v(:), 0, numel(v), numel(v));
    end
    
end
