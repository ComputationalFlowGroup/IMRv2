# Viscosity integrals

`imr.physics.viscosity` reduces the radial non-Newtonian stress to two integrals. One has a closed form; the other uses a fixed graded Gauss–Legendre rule.

With `gammadot = gammadot_R*(R/r)^3`, `gammadot_R = -2*Rdot/R`, and the substitution `u = R/r`, the integrals over `r` in `[R, Inf)` become integrals over `u` in `[0, 1]`, with `x_R = lambda*|gammadot_R|`:

- `I1 = int_R^Inf gammadot/r*f dr = gammadot_R*int_0^1 u^2*f(x_R*u^3) du`
- `I2 = int_R^Inf (f + x*f')/r^4 dr = R^-3*int_0^1 u^2*(f + x*f')(x_R*u^3) du = f(x_R)/(3*R^3)`

`I2` is exact: with `x = x_R*u^3`, `u^2*du = dx/(3*x_R)`, and `f + x*f' = d(x*f)/dx`, so the integral is `x_R*f(x_R)/(3*x_R)`.

`I1` uses a fixed Gauss–Legendre rule: 10 points on each of the 15 panels `[0, 2^-14], [2^-14, 2^-13], ..., [1/2, 1]`, 150 nodes, computed once. The panels shrink toward `u = 0` because at high shear rates the integrand changes on the scale `u ~ x_R^(-1/3)`. See [reported comparisons](../validation/convergence.md).
