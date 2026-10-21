# `f_display`

`f_display` prints a summary of the run to the command window before time integration. Both solvers call it. The summary shows:

- the solver mode;
- the radial model;
- the constitutive model;
- the viscosity model;
- the thermal, mass-transfer and vapor settings;
- how the stress is represented;
- the main dimensionless numbers (`Re8`, `De`, `Ca`, `LAM`).

## Inputs

`radial`, `bubtherm`, `medtherm`, `masstrans`, `stress`, `spectral`, `nu_model`, `eps3`, `Pv_star`, `Re8`, `De`, `Ca`, `LAM`, and `mode`, the solver mode (for example `nodal, Chebyshev`). Optionally, after `mode`, the non-Newtonian inputs `DRe`, `iDRe` (the factor of the non-Newtonian stress term), `v_lambda`, `v_nc` and `v_a`. Both solvers pass them.

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
- **Stress 3–5:** "stress integrals (ODEs)". The exceptions are stress 5 with `spectral = 1`, which prints "stress field, spectral modes", and stress 3–4 with `medtherm = 1`, which print "stress integral (ODE), and its local field for the medium heating" (see [`f_zener_field`](f_zener_field.md)).

## Viscosity label

The label is the model name, from `nu_model`:

| `nu_model` | Label | Parameters shown |
|---|---|---|
| 0 | Newtonian | |
| 1 | Carreau | `v_lambda`, `v_nc` |
| 2 | Carreau-Yasuda | `v_lambda`, `v_nc`, `v_a` |
| 3 | Powell-Eyring | `v_lambda` |
| 4 | Modified Powell-Eyring | `v_lambda`, `v_nc` |
| 5 | Cross | `v_lambda`, `v_nc` |
| 6 | Simplified Cross | `v_lambda` |
| 7 | Modified Cross | `v_lambda`, `v_nc`, `v_a` |

With an active model, a line after the dimensionless numbers gives `DRe` and the parameters the model uses.

The label also says when the model does not act as its name suggests:

- **Inactive with `DRe = 0`** (`mu0 = mu8`): "selected but inactive (DRe = 0, mu0 = mu8): Newtonian".
- **Inactive with `LAM = 0` for stress 3–5:** "selected but inactive (LAM = 0, no solvent share for stress 3-5): Newtonian". For these models the non-Newtonian viscosity is the solvent share `LAM`.
- **Stress 3–5, active:** ", solvent share (LAM)".
- **`v_lambda = 0`:** ", v_lambda = 0 so constant viscosity mu0". The shear rate then drops out, and the viscosity is the zero-shear value.
