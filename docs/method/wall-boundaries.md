# Wall boundary conditions

The nodal solver enforces wall heat balance algebraically instead of integrating wall values. `f_wall_theta_exact` handles thermal-only cases; `f_wall_mass_closure` couples heat and vapor mass transfer.

## Thermal-only wall

The nested `f_wall_theta_exact` in `src/+imr/+solvers/build_nodal.m` solves
`c1*Tw + c2*theta_w + K = 0`. `K` collects the other nodes' wall-row contributions. Set `ag = alpha_g`, `bg = beta_g` and

$$
u=\sqrt{(a_g+b_g)^2+2a_g\theta_w},\qquad
T_w=\frac{u-b_g}{a_g},\qquad
\theta_w=\frac{u^2-(a_g+b_g)^2}{2a_g}.
$$

Here $a_g$ and $b_g$ are the code variables `ag` and `bg`; $T_w$ is wall temperature divided by ambient temperature. The Kirchhoff variable $\theta_w$ is the integral of dimensionless conductivity from ambient temperature to $T_w$, so ambient conditions give $T_w=1$ and $\theta_w=0$. The nonnegative quantity $u=a_gT_w+b_g$ is conductivity at the wall; choosing an admissible root is a physical constraint. This inverse requires positive $a_g$; unsupported zero-slope thermal inputs are rejected during validation.

Then `c2*u^2 + 2*c1*u + q = 0`, where
`q = -2*c1*bg - c2*(ag+bg)^2 + 2*ag*K`.

| Coefficient | Meaning |
|---|---|
| `c1` | Medium wall-row coefficient multiplying dimensionless wall temperature |
| `c2` | Bubble wall-row coefficient multiplying the Kirchhoff variable |
| `K` | Fixed contributions from all nonwall nodes; distinct from conductivity $u$ |
| `q` | Constant term after eliminating temperature in favor of conductivity |

All gradients use the moving, dimensionless coordinates and the flux scales in the prepared case. Uniform ambient fields give zero wall Kirchhoff variable. At zero quadratic coefficient the equation reduces to a linear solve; the algorithm below handles that limit explicitly. `imr.physics.thermal_inverse` supplies the common inverse transform and derivatives used by the wall residual and both Jacobians.

- Form both roots in a cancellation-free form and discard `u < 0`.
- If two roots remain, choose the wall `theta` closest to the previous wall value; do not hard-code the sign of `c2`.
- If `c2 = 0`, use `u = -q/(2*c1)`.
- If the discriminant is negative or no admissible root remains, try `fzero` on the residual.
- If that fails at an unphysical trial state, return the least-residual quadratic point so the integrator can reject the step.

The Jacobian uses `d(theta_w)/dK = -u/(c2*u + c1)` on the selected branch.

## Mass-transfer residual

The saturation law is `Pv(Tw) = C1_pv*exp(C2_pv/Tw)` with negative `C2_pv`. The coupled balance is transcendental, so the solver uses the bounded wall vapor mass fraction `kv_w`, not an unbounded temperature variable.

For a trial fraction, `f_wall_mass_residual` computes:

| Step | Expression |
|---|---|
| Vapor partial pressure | `ps = P*kv_w*Rvg/(kv_w*Rvg + 1 - kv_w)` |
| Saturation temperature | `Tw = C2_pv/log(ps/C1_pv)` |
| Kirchhoff temperature | `theta_w = 0.5*alpha_m*(Tw^2 - 1) + beta_m*(Tw - 1)` |
| Positive scale factor | `D = (kv_w*Rva + Rg)*Tw*(1 - kv_w)` |
| Scaled residual | `r = D*(grad_Tm(1)*Tw + grad_Tr(1)*theta_w + flux_fixed) + P*(grad_C(1)*kv_w + carrier_fixed)` |

Multiplying by positive `D` preserves the root without the reciprocal singularity at fraction 1. Fixed contributions are hoisted once per call:
`flux_fixed = grad_Tm(2:end)*Tm(iwall_m) + grad_Tr(2:end)*theta(iwall_b)` and
`carrier_fixed = grad_C(2:end)*kv(iwall_b)`.

Node order matters. Medium values run outward from the wall; bubble temperature and vapor values run inward. With finite differences, those are `Tm(2:3)`, `theta(end-1:-1:end-2)` and `kv(end-1:-1:end-2)`. Chebyshev uses every nonwall node in the corresponding order.

## Root selection and fallback

| Path | Algorithm |
|---|---|
| Warm start | Newton from previous root, analytic derivative, at most 6 steps, accepted below step size `1e-10` while every iterate stays inside the bracket |
| Bracket | `[1e-13, 1 - 1e-13]`; evaluate 64 equally spaced points and select the first sign-changing interval, twice |
| Polish | Three analytic Newton steps constrained to the final interval |
| No endpoint sign change | Raise a state/time error; do not bypass it |

The bracketed path costs two vectorized residual grids, three Newton evaluations and one final evaluation. It has fixed iteration counts. It is also used on the first call, or if warm-start iterates leave the bracket, become nonfinite or fail to converge. Batched columns can fall back independently. [Historical measurements](../validation/convergence.md#reported-wall-closure-comparisons) record earlier accuracy and cost comparisons.

## Mixture-coefficient limitation

The closure uses mixture coefficients from stored `kv(end)`, not trial `kv_w`. Because `kvdot(end) = 0`, that stored value remains the initial `kv0` during integration. This is an existing modeling limitation, not a lag updated each step; it requires separate scientific review.

Conductivity is `K = alpha_m*T + beta_m`; the transform and inverse use both coefficients as supplied. In the default case their pure-component sums are `alpha_g + beta_g = 1.0024` and `alpha_v + beta_v = 0.9976`, not exactly 1.

After convergence, the solver assigns `theta(end)`, reconstructs temperature with the stored mixture, then sets the local wall `kv(end)` to `kv_w`. Pressure, velocity, conductivity and vapor equations see the converged local wall value. Output assembly recovers wall values with the same closure.

## Implicit derivative

[`imr.physics.wall_mass_balance`](../reference/functions/f_wall_mass_bc.md) expresses the same balance in `theta_w` and returns its gradient. `imr.numerics.jacobian_nodal` uses `d(theta_w)/dx = -(dr/dx)/(dr/d(theta_w))` at the root. Positive residual scaling cancels in that ratio at a root, so the closure and derivative describe the same balance.
