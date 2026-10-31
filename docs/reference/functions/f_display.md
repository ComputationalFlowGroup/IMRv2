# `imr.solvers.describe_case`

`imr.solvers.describe_case` prints a summary of the run to the command window before time integration. Both solvers call it. The summary shows:

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
| 1, 2, 3, 4 | `Inf` | `Re8 = Inf` | inviscid fluid |
| 1, 2, 3, 4 | `Inf` | `Re8` finite | Newtonian fluid |
| 1, 2, 3, 4 | finite | `Re8 = Inf` | neo-Hookean (quadratic, 2nd-order Fung or Mooney-Rivlin) elastic solid (inviscid) |
| 1, 2, 3, 4 | finite | `Re8` finite | neo-Hookean (quadratic, 2nd-order Fung or Mooney-Rivlin) Kelvin–Voigt |
| 5, 6, 7, 8 | `Inf` | `LAM = 0` | linear Maxwell |
| 5, 6, 7, 8 | `Inf` | `LAM > 0` | linear Jeffreys |
| 5, 6, 7, 8 | finite | `LAM = 0` | linear Zener, with a neo-Hookean (5), quadratic (6), 2nd-order Fung (7) or Mooney-Rivlin (8) spring |
| 5, 6, 7, 8 | finite | `LAM > 0` | linear Zener with solvent viscosity, with the same spring |
| 9 | (always `Inf`) | `LAM = 0` | upper-convected Maxwell |
| 9 | (always `Inf`) | `LAM > 0` | Oldroyd-B |
| 10 | (always `Inf`) | | PTT (linear Phan-Thien-Tanner) |
| 11 | (always `Inf`) | | Giesekus |

- **Stress 5, 6, 7, 8 without a spring:** with `Ca = Inf`, each of these is the same model, since the spring law only affects the spring.
- **Stress field:** with `spectral = 1`, stress 9 adds "(spectral stress field)"; stress 10 and 11 always add it, since `spectral = 1` is forced for them.
- **Viscosity line:** for stress 0 it reads "none", since there is no viscous stress.
- **`eps3`:** shown as ", eps3 = <value>" whenever nonzero (stress 11, the Giesekus mobility factor).

## Stress representation

- **Stress 0–4:** "closed form (no stress history)".
- **Stress 5–11:** "stress integrals (ODEs)". The exceptions are stress 9–11 with `spectral = 1` (always on for 10, 11), which prints "stress field, spectral modes", and stress 5–8 with `medtherm = 1`, which print "stress integral (ODE), and its local field for the medium heating" (see [`imr.physics.zener_field`](f_zener_field.md)).

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
- **Inactive with `LAM = 0` for the Maxwell/Zener/UCM family (stress 5–11):** "selected but inactive (LAM = 0, no solvent share for the Maxwell/Zener/UCM family): Newtonian". For these models the non-Newtonian viscosity is the solvent share `LAM`.
- **Stress 5–11, active:** ", solvent share (LAM)".
- **`v_lambda = 0`:** ", v_lambda = 0 so constant viscosity mu0". The shear rate then drops out, and the viscosity is the zero-shear value.
