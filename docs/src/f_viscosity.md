# `f_viscosity.m` Documentation

## Overview
`f_viscosity` computes the non-Newtonian viscosity integrals that enter the stress integral of the radial equations. The viscosity is `mu = mu8 + (mu0 - mu8)*f(x)` with `x = lambda*|gammadot|`, and each model supplies `f(x)` and `x*f'(x)`.

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
With `gammadot = gammadot_R*(R/r)^3`, `gammadot_R = -2*Rdot/R`, and the substitution `u = R/r`, the integrals over `r` in `[R, Inf)` become smooth integrals over `u` in `[0, 1]`:

- `I1 = int_R^Inf gammadot/r*f dr = gammadot_R*int_0^1 u^2*f(x_R*u^3) du`
- `I2 = int_R^Inf (f + x*f')/r^4 dr = R^-3*int_0^1 u^2*(f + x*f')(x_R*u^3) du`

with `x_R = lambda*|gammadot_R|`. They are evaluated with `integral` (tolerances `1e-8`).

## Outputs

| Name     | Description |
|----------|-------------|
| `f`      | Viscosity factor at the bubble wall |
| `intf`   | `-I1`; enters the stress integral as `-6*intf/DRe` |
| `dintf`  | `Rdot` part of `d(intf)/dt`, `-d(gammadot_R*R^3)/dt*I2 + h` with the Leibniz term `h = f*gammadot_R*Rdot/R` |
| `ddintf` | `Rddot` coefficient of `d(intf)/dt`, `2*R^2*I2` |
