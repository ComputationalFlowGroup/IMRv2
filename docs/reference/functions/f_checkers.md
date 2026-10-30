# `imr.configuration.validate_case`

`imr.configuration.validate_case` validates a run before it starts. `imr.configuration.prepare_case` calls it with a struct holding every model setting and dimensionless number. `imr.configuration.schema` supplies public/internal names and discrete choice/count rules; physical ranges and coupled restrictions remain here. Incompatible or out-of-range inputs raise an error; valid but suspicious inputs raise a warning. It returns the struct with the settings tied to the stress model made consistent (see below).

## Checks

### Model flags

- `radial` in 1–6, `stress` in 0–8, `nu_model` in 0–7, `wave_type` in -3–3, `method` in {3, 15, 23}, `jacobian` in 0–3, `scheme` in 0–1.
- `bubtherm`, `medtherm`, `vapor`, `masstrans`, `collapse`, `spectral`, `progdisplay` are 0 or 1; `eps3` in [0, 0.5]. Nonzero `eps3` warns unless `stress = 7` (Giesekus), the only model it enters; `stress = 7` with `eps3 = 0` warns too, since it then reduces to stress 5.

### Compatibility of the conditions

| Condition | Result |
|-----------|--------|
| `medtherm = 1` without `bubtherm` | error |
| `masstrans = 1` without `vapor` or `bubtherm` | error |
| Mie–Grüneisen (`radial = 5, 6`) with `bubtherm`, `medtherm` or `masstrans` | error: polytropic bubble only |
| Enthalpy forms (`radial = 3, 4`) with `medtherm = 1` | error: the liquid is taken at the far-field temperature |
| `collapse = 1` without `bubtherm`, `medtherm`, `vapor` and `masstrans` | error |
| `collapse = 1` with stress 5, 6 or 7 and no `stress0` (the stress at the maximum radius is not computed for UCM/Oldroyd-B, PTT or Giesekus) | error |
| `stress = 3, 4, 5, 6, 7` with `De < 1e-4` | error: the model is then its `De -> 0` limit, so the message names the simpler model to run: stress 1 for stress 3, stress 2 for stress 4, stress 1 with `Ca = Inf` for stress 5–7 (with `nu_model > 0`, `DRe/LAM` in place of `DRe`). Stress 8 has no `De` and is exempt. |
| `nu_model > 0` with `stress = 0` | error: the non-Newtonian viscosity needs a viscous stress model (any nonzero `stress`) |
| `stress = 5` with `medtherm = 1` and `spectral = 0` | error in `imr.solve_nodal`; `imr.solve_spectral` sets `spectral = 1` and warns |
| `stress = 6, 7` in `imr.solve_nodal` | error: no closed-form (two-ODE) reduction, `imr.solve_spectral` only |
| `stress = 6, 7` with `spectral = 0` in `imr.solve_spectral` | warning; sets `spectral = 1` |
| `stress = 3, 4` with `(1-LAM)/Re8 < De/Ca` (negative dashpot viscosity of the Zener branch) | warning |
| both `tvector` and `tfin` | error |
| `tvector` not finite and strictly increasing, or starting below 0 | error |
| `Nt < 3` with `bubtherm = 1`, or `Mt < 3` with `medtherm = 1` | error |
| `divisions` not a nonnegative integer | error |
| `vapor = 1` with `Pv >= 1` (boiling liquid) | warning `IMR:parameters` |
| `masstrans = 1` in the spectral solver | error: use `imr.solve_nodal` |
| `jacobian = 3` in the spectral solver | error: batched finite differences are in `imr.solve_nodal` only |
| `method = 3` (CVODES) before MATLAB R2024a, which lacks the `ode` object's `cvodesstiff` solver | error: use 15 or 23 |
| `spectral = 1` (spectral stress representation) in the nodal solver | error: `imr.solve_spectral` only |
| `scheme ~= 0` in the spectral solver | error: `scheme` is for `imr.solve_nodal` |

### Dimensionless numbers

- `Re`, `We`, `Ca` positive; `Inf` is allowed (inviscid, no surface tension, no elasticity).
- `De`, `LAM`, `alphax`, `kappa - 1` nonnegative; `LAM > 1` for the Maxwell/Zener and UCM family (stress 3–7) warns (stress 8 has no `LAM`).
- `Cstar`, `GAMa`, `nstate`, `hugoniot_s`, `Req`, `Lv`, `Lt`, tolerances positive.
- The equilibrium gas pressure `P0_star = (1 + 1/(We*Req) - Pv)*Req^n` must be nonnegative.
- Scalar parameter types/shapes are checked before derived arithmetic; values must be real and finite except the permitted infinities above. Effective numeric settings use `double`, including integer-valued inputs. Thermal bubbles require positive initial total pressure; mass transfer also requires positive initial gas pressure and a representable nonzero gas fraction.
- Thermal (`chi`, `Foh`, `kmstar`, conductivity coefficients) and mass-transfer (`Fom`, `lheat`, `rv`, `ra`, `pv`, `thetav`) numbers are checked when the corresponding model is on; `rv` and `ra` also with `collapse = 1`.
- Waveform: `tw > 0` for Gaussian, `om > 0` and `mn >= 0` for histotripsy, `tw >= 0` for Heaviside, and a warning when `tvector` extends past the measured waveform data.
- Non-Newtonian: `DRe`, `v_lambda`, `v_nc` nonnegative; `v_a > 0` for Carreau–Yasuda and modified Cross; `v_nc > 0` for Cross and modified Cross; `v_nc <= 1` and `v_lambda > 0` for modified Powell–Eyring; warning when `DRe = 0`, or with the Maxwell/Zener and UCM family (stress 3–7) when `LAM = 0` (the non-Newtonian viscosity is the solvent share `LAM`, so it vanishes; stress 8 has no `LAM` and is exempt).

### Dimensional-looking inputs

Warning `IMR:dimensional` when `Req < 1e-3`, the final time `< 1e-3`, `Cstar < 1`, or the Gaussian width `< 1e-3`: values this small usually mean an SI value was passed instead of a dimensionless one.

### Settings consistent with the stress model

After the checks, `imr.configuration.validate_case` adjusts and returns:

- `stress <= 4` or `stress = 8`: `spectral = 0`; `Nv = 1` for stress 3, 4 (one stress state) and `Nv = 0` otherwise (including stress 8, which carries no stress state at all).
- `stress = 5` (UCM, Oldroyd-B): `Nv` stress modes with the spectral representation, one otherwise, and `Ca = Inf` (no elastic contribution).
- `stress = 6, 7` (PTT, Giesekus): `spectral` forced to 1; `Nv` stress modes with the spectral representation, and `Ca = Inf` (no elastic contribution).
- `stress = 3, 4`: `Nv = 1`, except in `imr.solve_spectral` with `medtherm = 1`, where `Nv` is the number of modes of the local polymer stress field.
- `stress = 8` (linear Kelvin-Voigt) keeps `Ca` as given, unlike stress 5–7.

The type, finiteness and size of a user `stress0`, and the initial mixture conductivity, are checked in `imr.configuration.prepare_case`. `imr.solvers.build` then checks the assembled initial derivative for real, finite values before integration. See the [settings audit](../../validation/settings.md) and [`compatibility`](../compatibility.md) for coverage and limits.
