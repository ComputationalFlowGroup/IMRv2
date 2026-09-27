# `f_display.m` Documentation

## Overview
`f_display` prints a summary of the run to the command window before time integration. Both solvers call it. The summary shows:

- the solver mode;
- the radial model;
- the constitutive model;
- the viscosity model;
- the thermal, mass-transfer and vapor settings;
- how the stress is represented;
- the main dimensionless numbers (`Re8`, `De`, `Ca`, `LAM`).

## Inputs
`radial`, `bubtherm`, `medtherm`, `masstrans`, `stress`, `spectral`, `nu_model`, `eps3`, `Pv_star`, `Re8`, `De`, `Ca`, `LAM`, and `mode`, the solver mode (for example `nodal, Chebyshev`).

## Constitutive model label

The label depends on the numbers as well as on `stress`:

| `stress` | `Ca` | `Re8` or `LAM` | Label |
|---|---|---|---|
| 0 | | | no stress applied |
| 1, 2 | `Inf` | `Re8 = Inf` | inviscid fluid |
| 1, 2 | `Inf` | `Re8` finite | Newtonian fluid |
| 1, 2 | finite | `Re8 = Inf` | neo-Hookean (or quadratic neo-Hookean) elastic solid (inviscid) |
| 1, 2 | finite | `Re8` finite | neo-Hookean (or quadratic neo-Hookean) Kelvin–Voigt |
| 3, 4 | `Inf` | `LAM = 0` | linear Maxwell |
| 3, 4 | `Inf` | `LAM > 0` | linear Jeffreys |
| 3, 4 | finite | `LAM = 0` | linear Zener, with a neo-Hookean (stress 3) or quadratic neo-Hookean (stress 4) spring |
| 3, 4 | finite | `LAM > 0` | linear Zener with solvent viscosity, with the same spring |
| 5 | (always `Inf`) | `LAM = 0` | upper-convected Maxwell |
| 5 | (always `Inf`) | `LAM > 0` | Oldroyd-B |

- **Stress 3 and 4 without a spring:** with `Ca = Inf`, stress 3 and 4 are the same model, since the quadratic law only affects the spring.
- **Stress field:** with `spectral = 1`, stress 5 adds "(spectral stress field)".
- **Viscosity line:** for stress 0 it reads "none", since there is no viscous stress.

## Stress representation

- **Stress 0–2:** "closed form (no stress history)".
- **Stress 3–5:** "stress integrals (ODEs)". The exception is stress 5 with `spectral = 1`, which prints "stress field, spectral modes".
