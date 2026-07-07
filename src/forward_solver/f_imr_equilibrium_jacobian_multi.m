% file f_imr_equilibrium_jacobian_multi.m
% brief finite-difference first-order Jacobian for multi-mode equilibria
function [J, eigvals] = f_imr_equilibrium_jacobian_multi(xeq, varargin)
%F_IMR_EQUILIBRIUM_JACOBIAN_MULTI Linearize [Rdot; Rddot; epsdot; epsddot].
%
%   xeq can be [R; epsilon(:)] or [R; Rdot; epsilon(:); epsdot(:)].

    [fd_step, solver_args] = strip_jacobian_args(varargin);

    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, pert_opts] = ...
        evalc('f_call_params(solver_args{:});');
    np = length(pert_opts.n);

    xeq = xeq(:);
    if length(xeq) == np + 1
        y0 = [xeq(1); 0; xeq(2:end); zeros(np,1)];
    elseif length(xeq) == 2*np + 2
        y0 = xeq;
    else
        error(['xeq must be [R; epsilon(:)] or ' ...
            '[R; Rdot; epsilon(:); epsdot(:)].']);
    end

    nstate = length(y0);
    J = zeros(nstate);
    for k = 1:nstate
        h = fd_step*(1 + abs(y0(k)));
        yp = y0;
        ym = y0;
        yp(k) = yp(k) + h;
        ym(k) = ym(k) - h;
        J(:,k) = (rhs_full(yp) - rhs_full(ym))./(2*h);
    end

    eigvals = eig(J);

    function dy = rhs_full(y)
        R = y(1);
        Rdot = y(2);
        ep = y(3:(2+np));
        epdot = y((3+np):(2+2*np));
        D = f_imr_rhs_diagnostics([R; Rdot; NaN; ep; epdot], ...
            solver_args{:});
        dy = [Rdot; D.Rddot; epdot; D.epddot(:)];
    end
end

function [fd_step, solver_args] = strip_jacobian_args(args)
    if mod(length(args),2) == 1
        error('Inputs after xeq must be name/value pairs.');
    end

    fd_step = 1e-6;
    keep = true(size(args));
    for k = 1:2:length(args)
        name = args{k};
        if strcmpi(name, 'fd_step')
            fd_step = args{k+1};
            keep(k:k+1) = false;
        end
    end
    solver_args = args(keep);
end
