# `f_checkers.m` Documentation

## Overview
`f_checkers` validates a run before it starts. `f_call_params` calls it with a struct holding every model setting and dimensionless number. Incompatible or out-of-range inputs raise an error; valid but suspicious inputs raise a warning.

## Checks

### Model flags
- `radial` in 1–6, `stress` in 0–5, `nu_model` in 0–2, `wave_type` in -3–3, `method` in {15, 23, 45}, `jacobian` in 0–2.
- `bubtherm`, `medtherm`, `vapor`, `masstrans`, `collapse`, `spectral`, `progdisplay` are 0 or 1; `eps3` in [0, 0.5].

### Compatibility of the conditions
| Condition | Result |
|-----------|--------|
| `medtherm = 1` without `bubtherm` | error |
| `masstrans = 1` without `vapor` or `bubtherm` | error |
| Mie–Grüneisen (`radial = 5, 6`) with `bubtherm`, `medtherm` or `masstrans` | error: polytropic bubble only |
| `collapse = 1` without the full model with vapor | error |
| `stress = 3, 4, 5` with `De <= 0` | error |
| `stress >= 3` with `medtherm = 1` (FD dissipation) | warning `IMR:dissipation` |
| both `tvector` and `tfin` | error |

### Dimensionless numbers
- `Re`, `We`, `Ca` positive; `Inf` is allowed (inviscid, no surface tension, no elasticity).
- `De`, `LAM`, `alphax`, `kappa - 1` nonnegative; `LAM > 1` for stress models warns.
- `Cstar`, `GAMa`, `nstate`, `hugoniot_s`, `Req`, `Lv`, `Lt`, tolerances positive.
- The equilibrium gas pressure `P0_star = (1 + 1/(We*Req) - Pv)*Req^n` must be nonnegative.
- Thermal (`chi`, `Foh`, `kmstar`, conductivity coefficients) and mass-transfer (`Fom`, `lheat`, `rv`, `ra`, `pv`, `thetav`) numbers are checked when the corresponding model is on.
- Waveform: `tw > 0` for Gaussian, `om > 0` for histotripsy, measured-waveform coverage of `tvector`.
- Non-Newtonian: `DRe`, `v_lambda`, `v_nc` nonnegative, `v_a > 0` for Carreau–Yasuda, warning when `DRe = 0`.

### Dimensional-looking inputs
Warning `IMR:dimensional` when `Req < 1e-3`, the final time `< 1e-3`, `Cstar < 1`, or the Gaussian width `< 1e-3`: values this small usually mean an SI value was passed instead of a dimensionless one.
