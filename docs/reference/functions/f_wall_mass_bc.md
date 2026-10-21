# `f_wall_mass_bc`

`f_wall_mass_bc(v, c)` evaluates the same balance in terms of `theta_w`, with its gradient with respect to `v = [theta_w, Tm(iwall_m), theta(iwall_b), kv(end), kv(iwall_b), P]`. `iwall_m` and `iwall_b` hold two nodes each for finite differences and every node for Chebyshev. `f_jacobian_nodal` uses the gradient at the converged root for the implicit derivative `d(theta_w)/dx = -(dr/dx)/(dr/d(theta_w))`. Scaling a residual by a positive factor does not change that ratio at a root, so the derivative is the one of the closure's residual as well.

See the [wall method](https://github.com/mrodrig6/imrv2/blob/main/method/wall-boundaries.md) for the root solve.
