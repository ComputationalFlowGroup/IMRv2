% file f_imr_equilibrium_jacobian.m
% brief finite-difference first-order Jacobian for a one-mode equilibrium
function [J, eigvals] = f_imr_equilibrium_jacobian(xeq, varargin)
%F_IMR_EQUILIBRIUM_JACOBIAN Linearize [Rdot; Rddot; epdot; epddot].
%
%   xeq can be [R; epsilon] or [R; Rdot; epsilon; epdot].

    [fd_step, solver_args] = strip_jacobian_args(varargin);

    xeq = xeq(:);
    if length(xeq) == 2
        y0 = [xeq(1); 0; xeq(2); 0];
    elseif length(xeq) == 4
        y0 = xeq;
    else
        error('xeq must be [R; epsilon] or [R; Rdot; epsilon; epdot].');
    end

    nstate = length(y0);
    J = zeros(nstate);
    for k = 1:nstate
        h = fd_step*(1 + abs(y0(k)));
        yp = y0;
        ym = y0;
        yp(k) = yp(k) + h;
        ym(k) = ym(k) - h;
        J(:,k) = (rhs4(yp) - rhs4(ym))./(2*h);
    end

    eigvals = eig(J);

    function dy = rhs4(y)
        D = f_imr_rhs_diagnostics([y(1); y(2); NaN; y(3); y(4)], ...
            solver_args{:});
        dy = [y(2); D.Rddot; y(4); D.epddot(1)];
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
