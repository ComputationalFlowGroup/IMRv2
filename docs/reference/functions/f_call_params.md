# `imr.configuration.legacy_options`

`imr.configuration.legacy_options` preserves the legacy grouped-vector interface. It delegates preparation to `imr.configuration.prepare_case` and serializes the resulting named fields into the original output order. The forward solvers use the [named case/system contracts](../internal-contracts.md) directly.

## Function Signature

```text
[eqns_opts, solve_opts, init_opts, init_stress, tspan_opts, out_opts, acos_opts, wave_opts, sigma_opts, thermal_opts, mass_opts, inputs] = imr.configuration.legacy_options(solver, varargin)
```


## Inputs

The first argument names the calling solver, `'nodal'` or `'spectral'` (used by `imr.configuration.validate_case` for solver-specific rules). The function then accepts **dimensionless parameter-value pairs** through `varargin`; the recognized names, their definitions and defaults are listed in [`imr.configuration.defaults`](../../model/parameters.md). Dimensional (SI) values are not accepted: an unrecognized name (for example `mu`, `g`, `r0`, `dimout`) raises an error pointing to [`imr.nondimensionalize`](../../api/f_dim2nondim.md), which converts physical properties into these inputs.

- Model flags: `radial`, `bubtherm`, `medtherm`, `stress`, `eps3`, `vapor`, `masstrans`, `collapse`.
- Solver: `method` (a code, or the name `'ode23tb'`, `'ode15s'` or `'cvodes'`), `spectral`, `scheme`, `divisions`, `nt`, `mt`, `nv`, `lv`, `lt`, `reltol`, `abstol`, `jacobian`, `tfin`, `tvector`, `progdisplay`.
- Initial conditions: `req`, `rdotzero`, `p0star`, `stress0`.
- Acoustics: `cstar`, `gama`, `nstate`, `hugoniot_s`.
- Waveform: `wave_type`, `ee`, `om`, `tw`, `dt`, `mn`.
- Stress and viscosity: `we`, `re`, `ca`, `de`, `lam`, `alphax`, `dre`, `v_a`, `v_nc`, `v_lambda`, `nu_model`.
- Thermal: `kappa`, `foh`, `br`, `alpha_g`, `beta_g`, `alpha_v`, `beta_v`, `chi`, `kmstar`.
- Mass transfer and vapor pressure: `fom`, `lheat`, `rv`, `ra`, `pv`, `thetav`.

After reading the inputs, the function computes the derived quantities (initial gas pressure `P0_star`, initial vapor concentration, `JdotA`, ...) and calls [`imr.configuration.validate_case`](f_checkers.md) to validate the case.


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
- `inputs`: Requested/default input struct and final time vector; see [input-record limits](../../running/output.md).


## Notes

- This adapter remains available to existing callers and test generators. Forward solvers call `imr.configuration.prepare_case` directly.
- Critical inputs must be matched correctly; otherwise, the function throws errors.
- A struct from [`imr.default_inputs`](../../api/f_imr_inputs.md) as the first input becomes the name/value pairs of its fields that differ from the defaults, ahead of any pairs after it. The last output is the struct of the run's inputs, which the solvers return as their 8th output.
- Defaults come from [`imr.configuration.defaults`](../../model/parameters.md). `imr.configuration.schema` maps public names to internal fields and supplies unconditional discrete rules; `imr.configuration.validate_case` owns coupled physical constraints. `casefile` is no longer accepted.
- `imr.configuration.validate_case` and `imr.physics.initial_stress` are required for full execution.
- A user `stress0` must have 0 entries for stress 0–2, 1 for stress 3–4, and `2*Nv` for stress 5 (2 with `spectral = 0`).
- With `collapse = 1`, the equilibrium radius is recomputed from the vapor and gas content (printed as `New equilibrium radius`), and the initial wall velocity is set to zero. The earlier consistency measurement is preserved in [reported comparisons](../../validation/convergence.md#reported-collapse-setup-comparison).
- With `collapse = 1` and stress 3 or 4, the initial stress comes from [`imr.physics.initial_stress`](f_init_stress.md), with `alphax` for stress 4 (quadratic) and 0 for stress 3, unless `stress0` is given. Stress 0–2 carry no stress history. Stress 5 needs `stress0` (`imr.configuration.validate_case` rejects it otherwise).
