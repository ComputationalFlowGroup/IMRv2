% file f_imr_static_residual.m
% brief static two-equation residual for one-mode equilibria
function F = f_imr_static_residual(x, varargin)
%F_IMR_STATIC_RESIDUAL Return [Rddot; epddot] for a one-mode static state.
%
%   x = [R; epsilon], with R and epsilon nondimensional.

    x = x(:);
    if length(x) ~= 2
        error('x must be [R; epsilon].');
    end

    D = f_imr_rhs_diagnostics([x(1); 0; NaN; x(2); 0], varargin{:});
    F = [D.Rddot; D.epddot(1)];
end
