# `f_call_params.m` Documentation

## Overview
`f_call_params` sets up the simulation environment for the Inertial Microcavitation Rheometry (IMR) solver.
It loads a case file or the default case, applies the user's dimensionless overrides, validates them with `f_checkers`, and packages all simulation parameters needed for the solvers.

## Function Signature
```matlab
[eqns_opts, solve_opts, init_opts, init_stress, tspan_opts, out_opts, acos_opts, wave_opts, sigma_opts, thermal_opts, mass_opts] = f_call_params(solver, varargin)
```

---

## Inputs

The first argument names the calling solver, `'nodal'` or `'spectral'` (used by `f_checkers` for solver-specific rules). The function then accepts **dimensionless parameter-value pairs** through `varargin`; the recognized names, their definitions and defaults are listed in [`default_case`](default_case.md). Dimensional (SI) values are not accepted: an unrecognized name (for example `mu`, `g`, `r0`, `dimout`) raises an error pointing to [`toolchain/f_dim2nondim.m`](../toolchain/f_dim2nondim.md), which converts physical properties into these inputs.

- Model flags: `radial`, `bubtherm`, `medtherm`, `stress`, `eps3`, `vapor`, `masstrans`, `collapse`.
- Solver: `method`, `spectral`, `scheme`, `divisions`, `nt`, `mt`, `nv`, `lv`, `lt`, `reltol`, `abstol`, `jacobian`, `tfin`, `tvector`, `progdisplay`.
- Initial conditions: `req`, `rdotzero`, `p0star`, `stress0`.
- Acoustics: `cstar`, `gama`, `nstate`, `hugoniot_s`.
- Waveform: `wave_type`, `ee`, `om`, `tw`, `dt`, `mn`.
- Stress and viscosity: `we`, `re`, `ca`, `de`, `lam`, `alphax`, `dre`, `v_a`, `v_nc`, `v_lambda`, `nu_model`.
- Thermal: `kappa`, `foh`, `br`, `alpha_g`, `beta_g`, `alpha_v`, `beta_v`, `chi`, `kmstar`.
- Mass transfer and vapor pressure: `fom`, `lheat`, `rv`, `ra`, `pv`, `thetav`.

After reading the inputs, the function computes the derived quantities (initial gas pressure `P0_star`, initial vapor concentration, `JdotA`, ...) and calls [`f_checkers`](f_checkers.md) to validate the case.

---

## Outputs

The function returns the following grouped outputs:

- `eqns_opts`: Flags for physical models (radial, thermal, stress, mass transfer).
- `solve_opts`: Solver settings (method, divisions, spectral options).
- `init_opts`: Initial conditions in non-dimensional form (radius, velocity, pressure).
- `init_stress`: Initial stress tensor (depends on stress model).
- `tspan_opts`: Simulation time span or vector (non-dimensional).
- `out_opts`: Display control (`progdisplay`).
- `acos_opts`: Acoustic medium properties (sound speed, compressibility).
- `wave_opts`: Driving waveform parameters and optional preloaded data.
- `sigma_opts`: Viscoelastic material properties and dimensionless numbers.
- `thermal_opts`: Non-dimensional thermal transport coefficients.
- `mass_opts`: Mass transfer parameters and non-dimensional vapor concentration.

---

## Notes

- `f_call_params` must be called before running the forward solver.
- Critical inputs must be matched correctly; otherwise, the function throws errors.
- The defaults come from the function [`default_case`](default_case.md) and are copied into local variables one by one, then overridden by the name/value pairs. `casefile` is no longer accepted.
- `f_checkers.m` and `f_init_stress.m` are required for full execution.
- A user `stress0` must have 0 entries for stress 0–2, 1 for stress 3–4, and `2*Nv` for stress 5 (2 with `spectral = 0`).
- With `collapse = 1`, the equilibrium radius is recomputed from the vapor and gas content (printed as `New equilibrium radius`), and the initial wall velocity is set to zero. The recomputed value reproduces the input `Req` (checked for `Req` 0.05–0.9 and `Pv` 0.03 and 0.3), so the initial bubble pressure built from the input `Req` is consistent with it.
- With `collapse = 1` and stress 3 or 4, the initial stress comes from [`f_init_stress`](f_init_stress.md), with `alphax` for stress 4 (quadratic) and 0 for stress 3, unless `stress0` is given.

