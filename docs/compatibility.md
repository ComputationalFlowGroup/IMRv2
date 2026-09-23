# Model compatibility

This page lists which model options work in each solver. `f_checkers` enforces every "rejected" entry below before a run starts, with an `INPUT ERROR` that names the conflict.

Legend: **yes** runs, **no** rejected by `f_checkers`, **planned** not implemented yet.

## Solvers

| Solver | Discretization | State |
|---|---|---|
| `f_imr_fd` | Finite differences, node values | Full capability |
| `f_imr_spectral` | Chebyshev collocation, coefficients | No mass transfer |
| `f_imr_nodal` | Node values, FD or Chebyshev operator from `f_diff_mat` | Planned. Only the operator module `f_diff_mat` exists. |

## Thermal models

Rows are the thermal configurations. The columns give the radial models that each solver accepts for them.

| Configuration | `bubtherm` | `medtherm` | `vapor` | `masstrans` | `f_imr_fd` | `f_imr_spectral` |
|---|---|---|---|---|---|---|
| Polytropic bubble | 0 | 0 | 0 | 0 | radial 1–6 | radial 1–6 |
| Polytropic bubble with vapor pressure | 0 | 0 | 1 | 0 | radial 1–6 | radial 1–6 |
| Heat transfer in the bubble | 1 | 0 | 0/1 | 0 | radial 1–4 | radial 1–4 |
| Heat transfer in the bubble and the medium | 1 | 1 | 0/1 | 0 | radial 1–4 | radial 1–4 |
| Vapor mass transfer | 1 | 0/1 | 1 | 1 | radial 1–4 | no |

Radial models: 1 Rayleigh–Plesset, 2 Keller–Miksis (pressure), 3 Keller–Miksis (enthalpy, Tait), 4 Gilmore (Tait), 5 Keller–Miksis (Mie–Grüneisen), 6 Gilmore (Mie–Grüneisen).

Always rejected:
- `medtherm = 1` without `bubtherm = 1`.
- `masstrans = 1` without `vapor = 1`, or without `bubtherm = 1`.
- `radial = 5, 6` with any of `bubtherm`, `medtherm` or `masstrans` on.

## Stress models

Every stress model works with every thermal configuration and radial model that the solver accepts.

| `stress` | Model | `f_imr_fd` | `f_imr_spectral` |
|---|---|---|---|
| 0 | None | yes | yes |
| 1 | Neo-Hookean Kelvin–Voigt | yes | yes |
| 2 | Quadratic neo-Hookean Kelvin–Voigt | yes | yes |
| 3 | Maxwell / Jeffreys / Zener (neo-Hookean) | yes | yes |
| 4 | Same as 3 (quadratic neo-Hookean) | yes | yes |
| 5, `spectral = 0` | UCM / Oldroyd-B, two stress states | yes | yes |
| 5, `spectral = 1` | UCM / Oldroyd-B, spectral stress field with `Nv` modes | no | yes |

Notes:
- `stress = 3, 4, 5` require `De > 0`.
- `stress = 5` sets `Ca = Inf`, so it has no elastic contribution.
- For `stress = 0–4`, `f_checkers` resets `spectral` to 0.
- With `medtherm = 1` and `stress >= 3`, the medium heating uses the Kelvin–Voigt closed form, which is inconsistent for history-dependent stress, and `f_checkers` warns. The spectral stress field (`stress = 5`, `spectral = 1`) computes the heating from the stress field itself.

## Other options

| Option | `f_imr_fd` | `f_imr_spectral` |
|---|---|---|
| `nu_model = 0–7` (non-Newtonian viscosity) | yes | yes |
| `wave_type = -3…3` (driving waveforms) | yes | yes |
| `eps3` (second-order strain term) | yes | yes |
| `method = 15, 23, 45` (`ode15s`, `ode23tb`, `ode45`) | yes | yes |
| `jacobian = 2` (analytic Jacobian) | yes with `bubtherm = 1`; falls back to the sparsity pattern otherwise | yes |
| `jacobian = 1` | sparsity pattern | dense numerical |
| `jacobian = 0` | dense numerical | dense numerical |
| `collapse = 1` | yes, requires `bubtherm = medtherm = masstrans = vapor = 1` | no, because it requires mass transfer |
| `stress0` (user initial stress) | yes | yes |

## Verification

Every combination in the tables was run with GNU Octave on main as of the merge of PR #4, except the two listed at the end:
- 504 runs for the thermal and stress tables: both solvers, radial 1–6, six thermal configurations and stress 0–5, with `spectral = 0, 1` for stress 5. Another 36 runs covered vapor with bubble-only heat transfer and mass transfer without medium heat transfer. Every accepted combination finished with a finite `R(t)`, and every rejected one stopped with an `INPUT ERROR`.
- The non-Newtonian viscosity, waveform, `eps3`, `method` and `jacobian` options were each run on a thermal Keller–Miksis case in both solvers.

Octave has no `ode23tb`, so `ode15s` stood in for it. Two entries could not be confirmed in Octave:
- **`collapse = 1`.** Octave's `fzero` throws an error where MATLAB's returns `NaN`, so the retry loop in `f_call_params` never runs. The equilibrium equation does have a root. It needs a MATLAB run.
- **`method = 15` with `f_imr_fd`.** Octave's `ode15s` crashes on the sparse Jacobian. This is an Octave problem, not a problem in the code.
