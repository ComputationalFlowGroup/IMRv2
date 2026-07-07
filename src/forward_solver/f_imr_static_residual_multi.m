% file f_imr_static_residual_multi.m
% brief static residual for multi-mode equilibria
function F = f_imr_static_residual_multi(x, varargin)
%F_IMR_STATIC_RESIDUAL_MULTI Return [Rddot; epddot(:)].
%
%   x = [R; epsilon(:)], with all quantities nondimensional.

    x = x(:);

    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, pert_opts] = ...
        evalc('f_call_params(varargin{:});');
    np = length(pert_opts.n);

    if length(x) ~= np + 1
        error('x must contain [R; epsilon(:)] for every configured mode.');
    end

    q = [x(1); 0; NaN; x(2:end); zeros(np,1)];
    D = f_imr_rhs_diagnostics(q, varargin{:});
    F = [D.Rddot; D.epddot(:)];
end
