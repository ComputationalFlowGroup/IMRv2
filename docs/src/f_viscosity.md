# `f_viscosity.m` Documentation

## Overview
`f_viscosity` computes the non-Newtonian viscosity integrals that enter the stress integral of the radial equations. The viscosity is `mu = mu8 + (mu0 - mu8)*f(x)` with `x = lambda*|gammadot|`, and each model supplies `f(x)`.

## Models (`nu_model`)

| `nu_model` | Model | `f(x)` | Parameters |
|---|---|---|---|
| 1 | Carreau | `(1 + x^2)^((n-1)/2)` | `v_nc = n` |
| 2 | Carreau–Yasuda | `(1 + x^a)^((n-1)/a)` | `v_a = a > 0`, `v_nc = n` |
| 3 | Powell–Eyring | `asinh(x)/x` | – |
| 4 | Modified Powell–Eyring | `log(1 + x)/x^n` | `v_nc = n <= 1`, `v_lambda > 0` |
| 5 | Cross | `1/(1 + x^n)` | `v_nc = n > 0` |
| 6 | Simplified Cross | `1/(1 + x)` | – |
| 7 | Modified Cross | `(1 + x^n)^(-a)` | `v_a = a > 0`, `v_nc = n > 0` |

`lambda` is the dimensionless time scale `v_lambda`. Every model reduces to the Newtonian viscosity `mu0` as `lambda -> 0` (model 4 with `n = 1`).

## Integral formulation
With `gammadot = gammadot_R*(R/r)^3`, `gammadot_R = -2*Rdot/R`, and the substitution `u = R/r`, the integrals over `r` in `[R, Inf)` become integrals over `u` in `[0, 1]`, with `x_R = lambda*|gammadot_R|`:

- `I1 = int_R^Inf gammadot/r*f dr = gammadot_R*int_0^1 u^2*f(x_R*u^3) du`
- `I2 = int_R^Inf (f + x*f')/r^4 dr = R^-3*int_0^1 u^2*(f + x*f')(x_R*u^3) du = f(x_R)/(3*R^3)`

`I2` is exact: with `x = x_R*u^3`, `u^2*du = dx/(3*x_R)`, and `f + x*f' = d(x*f)/dx`, so the integral is `x_R*f(x_R)/(3*x_R)`.

`I1` uses a fixed Gauss–Legendre rule: 10 points on each of the 15 panels `[0, 2^-14], [2^-14, 2^-13], ..., [1/2, 1]`, 150 nodes, computed once. The panels shrink toward `u = 0` because at high shear rates the integrand changes on the scale `u ~ x_R^(-1/3)`. Against an adaptive `integral` at `1e-13`, the rule is accurate to `1.3e-11` relative for all seven models, exponents 0.3–0.8, and `x_R` from `1e-4` to `1e9`.

Compared with the previous version, which called `integral` twice at `1e-8`, the outputs agree to `4e-11` and a call is about 10 times faster (58 µs against 613 µs in GNU Octave). A run with a non-Newtonian model takes about half the time.

## Outputs

| Name     | Description |
|----------|-------------|
| `f`      | Viscosity factor at the bubble wall |
| `intf`   | `-I1`; enters the stress integral as `-6*intf/DRe` |
| `dintf`  | `Rdot` part of `d(intf)/dt`, `-d(gammadot_R*R^3)/dt*I2 + h` with the Leibniz term `h = f*gammadot_R*Rdot/R` |
| `ddintf` | `Rddot` coefficient of `d(intf)/dt`, `2*R^2*I2` |
