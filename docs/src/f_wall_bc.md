# Bubble-wall boundary conditions

## `f_wall_thermal_bc.m`
Without mass transfer the wall heat flux balance `c1*Tw + c2*theta_w + K = 0`, with `Tw = (ag - 1 + u)/ag` and `theta_w = (u^2 - 1)/(2*ag)`, `u = sqrt(1 + 2*ag*theta_w)`, is quadratic in `u`. `f_wall_thermal_bc(K, c1, c2, alpha_g)` returns the positive root in closed form and `d(theta_w)/dK`.

## `f_wall_mass_bc.m`
With vapor mass transfer the balance adds the latent heat flux, which depends on the saturated wall vapor concentration `kvw(Tw, P)` through `Pv(Tw) = C1*exp(C2/Tw)`. The balance is transcendental in `Tw` (the unknown appears both algebraically and inside the exponential), so it has no closed-form root. `f_wall_mass_bc(v, c)` returns the residual and its gradient with respect to `v = [theta_w, Tm(2), Tm(3), theta(end-1), theta(end-2), kv(end), kv(end-1), kv(end-2), P]`. `f_imr_fd` solves it by a warm-started Newton iteration (typically 3–4 iterations), falling back to `fzero` only if Newton does not converge; `f_jacobian_fd` uses the gradient for the implicit derivative of `theta_w`.
