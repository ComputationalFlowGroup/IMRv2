# `f_zener_field.m` Documentation

## Overview

`f_zener_field` returns the time derivative of the local polymer stress difference `sp = tau_p,rr - tau_p,thetatheta` of the linear Maxwell, Jeffreys and Zener models (stress 3 and 4). Both solvers carry `sp` as a field when medium heat transfer is on, because the heating (`f_stress_dissipation`) needs the local stress, not only its integral `Z1`.

## Equation

On the grid `y = r/R`, the field obeys

    d(sp)/dt|_y = -ydot*d(sp)/dy - (sp - sKV)/De - 3*(Rdot/R)*(1 - y^-3)*sp
    ydot = (Rdot/R)*(y^-2 - y)
    sKV  = se - 6*etap*a,   a = (Rdot/R)*y^-3,   etap = (1 - LAM)/Re8

- `se` is the neo-Hookean elastic stress difference for stress 3, or the quadratic one for stress 4 (see [`f_stress_dissipation`](f_stress_dissipation.md)).
- The caller computes the advection `-ydot*d(sp)/dy` on its own grid.
- **Boundary conditions:** none at the wall, where `ydot = 0`. `sp` is 0 at infinity.
- **Initial condition:** `sp = 0` from rest.

## Why this law

It is the local law whose integral reproduces the solver's stress state exactly: `2*int_1^Inf sp/y dy = Z1/R^3`, where `Z1` obeys `De*Z1dot + Z1 = Ze - 4*etap*R^2*Rdot` (the specification, Section 6.2, proves this). So:

- `Z1` still drives `R(t)`, unchanged.
- The field only distributes the heating in `r`.
- The total heating is fixed by `Z1`.

## In the solvers

- **`f_imr_spectral`.** `sp` is expanded in the stress 5 basis `T_n(z) - 1` with `Nv` modes, stored after `Z1`, and evaluated at the medium points with `ZZT`. After a run, the solver prints the largest mismatch of `2*cdd'*e` against `Z1/R^3` as a check.
- **`f_imr_nodal`.** `sp` is stored at the medium nodes `yT` after `Z1`: the first node is the wall, and the last, at infinity, is held at 0. The advection reuses the `Rdot` part of the medium temperature advection, `xi1_sq./(Lt*R).*(Rdot./yT2.*om_yT3/2).*(dmatrix_Tm*sp)`, which equals `-ydot*d(sp)/dy`.
- **Collapse (nodal only).** The history of the field before `t = 0` is not computed. As a stopgap, the initial field takes the shape of the elastic stress difference of the initial radius, scaled so that the integral identity holds for the initial `Z1` from `f_init_stress`. Without elasticity (`Ca = Inf`), it takes the viscous shape `y^-3` instead.

## Verification (GNU Octave)

The largest error in `2*int sp/y dy - Z1/R^3` over a run, relative to `max |Z1/R^3|` (Keller–Miksis, `Req = 0.5`):

| Mode | Stress 3 | Stress 4 |
|---|---|---|
| Nodal, finite differences, `Mt = 40, 80, 160` | `5.8e-4`, `1.4e-4`, `3.5e-5` | `1.4e-3`, `3.5e-4`, `8.9e-5` |
| Nodal, Chebyshev, `Mt = 12, 20, 30` | `2.0e-5`, `4.9e-7`, `7.7e-8` | `9.1e-5`, `7.4e-6`, `4.7e-7` |
| Spectral, `Nv = 8, 16, 32` | `1.9e-3`, `3.6e-5`, `5.9e-7` | `7.6e-3`, `6.2e-4`, `1.5e-5` |

That is second-order convergence for finite differences and spectral convergence otherwise.
