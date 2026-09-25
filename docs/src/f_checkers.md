# `f_checkers.m` Documentation

## Overview
`f_checkers` validates a run before it starts. `f_call_params` calls it with a struct holding every model setting and dimensionless number. Incompatible or out-of-range inputs raise an error; valid but suspicious inputs raise a warning. It returns the struct with the settings tied to the stress model made consistent (see below).

## Checks

### Model flags
- `radial` in 1–6, `stress` in 0–5, `nu_model` in 0–7, `wave_type` in -3–3, `method` in {3, 15, 23}, `jacobian` in 0–3, `scheme` in 0–1.
- `bubtherm`, `medtherm`, `vapor`, `masstrans`, `collapse`, `spectral`, `progdisplay` are 0 or 1; `eps3` in [0, 0.5].

### Compatibility of the conditions
| Condition | Result |
|-----------|--------|
| `medtherm = 1` without `bubtherm` | error |
| `masstrans = 1` without `vapor` or `bubtherm` | error |
| Mie–Grüneisen (`radial = 5, 6`) with `bubtherm`, `medtherm` or `masstrans` | error: polytropic bubble only |
| Enthalpy forms (`radial = 3, 4`) with `medtherm = 1` | error: the liquid is taken at the far-field temperature |
| `collapse = 1` without `bubtherm`, `medtherm`, `vapor` and `masstrans` | error |
| `collapse = 1` with stress 5 and no `stress0` (the stress at the maximum radius is not computed for UCM and Oldroyd-B) | error |
| `stress = 3, 4, 5` with `De <= 0` | error |
| `nu_model > 0` with `stress = 0` | error: the non-Newtonian viscosity needs a viscous stress model (stress 1–5) |
| `stress >= 3` with `medtherm = 1` (FD dissipation) | warning `IMR:dissipation` |
| both `tvector` and `tfin` | error |
| `tvector` not finite and strictly increasing, or starting below 0 | error |
| `Nt < 3` with `bubtherm = 1`, or `Mt < 3` with `medtherm = 1` | error |
| `divisions` not a nonnegative integer | error |
| `vapor = 1` with `Pv >= 1` (boiling liquid) | warning `IMR:parameters` |
| `masstrans = 1` in the spectral solver | error: use `f_imr_nodal` |
| `spectral = 1` (spectral stress representation) in the nodal solver | error: `f_imr_spectral` only |
| `scheme ~= 0` in the spectral solver | error: `scheme` is for `f_imr_nodal` |

### Dimensionless numbers
- `Re`, `We`, `Ca` positive; `Inf` is allowed (inviscid, no surface tension, no elasticity).
- `De`, `LAM`, `alphax`, `kappa - 1` nonnegative; `LAM > 1` for stress models warns.
- `Cstar`, `GAMa`, `nstate`, `hugoniot_s`, `Req`, `Lv`, `Lt`, tolerances positive.
- The equilibrium gas pressure `P0_star = (1 + 1/(We*Req) - Pv)*Req^n` must be nonnegative.
- Thermal (`chi`, `Foh`, `kmstar`, conductivity coefficients) and mass-transfer (`Fom`, `lheat`, `rv`, `ra`, `pv`, `thetav`) numbers are checked when the corresponding model is on; `rv` and `ra` also with `collapse = 1`.
- Waveform: `tw > 0` for Gaussian, `om > 0` and `mn >= 0` for histotripsy, `tw >= 0` for Heaviside, and a warning when `tvector` extends past the measured waveform data.
- Non-Newtonian: `DRe`, `v_lambda`, `v_nc` nonnegative; `v_a > 0` for Carreau–Yasuda and modified Cross; `v_nc > 0` for Cross and modified Cross; `v_nc <= 1` and `v_lambda > 0` for modified Powell–Eyring; warning when `DRe = 0`.

### Dimensional-looking inputs
Warning `IMR:dimensional` when `Req < 1e-3`, the final time `< 1e-3`, `Cstar < 1`, or the Gaussian width `< 1e-3`: values this small usually mean an SI value was passed instead of a dimensionless one.

### Settings consistent with the stress model
After the checks, `f_checkers` adjusts and returns:

- `stress <= 4`: `spectral = 0`; `Nv = 1` for stress 3, 4 (one stress state) and `Nv = 0` otherwise.
- `stress = 5` (UCM, Oldroyd-B): `Nv` stress modes with the spectral representation, one otherwise, and `Ca = Inf` (no elastic contribution).

The size of a user `stress0` is checked in `f_call_params`, not here. See [`compatibility`](../compatibility.md) for the resulting grid of features that run together.
