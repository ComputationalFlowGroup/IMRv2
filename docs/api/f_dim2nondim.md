# `f_dim2nondim`

`f_dim2nondim` is a side function that converts dimensional (SI) physical properties into the dimensionless inputs taken by `f_imr_nodal` and `f_imr_spectral`. The solvers never call it: they accept only dimensionless inputs.

See the [complete quick start](../getting-started/quick-start.md) for a runnable SI-to-solver example.

## Inputs (SI)

These defaults come from `src/f_dim2nondim.m`. Give inputs as name/value pairs; repeated dimensional names use the last value. Explicit `du` overrides the viscosity difference derived from `mu0` and `mu`.

| Name | Meaning | Default | Units |
|---|---|---|---|
| `r0` | Initial radius | `50e-6` | m |
| `req` | Equilibrium radius | `r0/5` | m |
| `u0` | Initial wall velocity | `0` | m/s |
| `tfin` | Final time | `Solver default if omitted` | s |
| `tvector` | Requested times | `Solver default if omitted` | s |
| `rho8` | Far-field density | `1064` | kg/m³ |
| `p8` | Far-field pressure | `101325` | Pa |
| `c8` | Far-field sound speed | `1484` | m/s |
| `gam` | Tait pressure parameter | `3049.13e5` | Pa |
| `surft` | Surface tension | `0.07` | N/m |
| `g` | Shear modulus | `1e3` | Pa |
| `mu` | High-shear viscosity | `8.3283e-4` | Pa s |
| `mu0` | Zero-shear viscosity parameter | `mu` | Pa s |
| `du` | Viscosity difference | `mu0 - mu` | Pa s |
| `lambda1` | Relaxation time | `1e-7` | s |
| `lambda2` | Retardation time | `1e-8` | s |
| `v_lambda` | Viscosity time scale | `0` | s |
| `t8` | Far-field temperature | `298.15` | K |
| `atg` | Gas conductivity slope | `5.28e-5` | W/(m K²) |
| `btg` | Gas conductivity intercept | `1.165e-2` | W/(m K) |
| `atv` | Vapor conductivity slope | `3.30e-5` | W/(m K²) |
| `btv` | Vapor conductivity intercept | `1.742e-2` | W/(m K) |
| `km` | Medium conductivity | `0.55` | W/(m K) |
| `cp` | Medium heat capacity | `4.181e3` | J/(kg K) |
| `dm` | Medium thermal diffusivity | `km/(rho8*cp)` | m²/s |
| `dmass` | Mass diffusivity | `24.2e-6` | m²/s |
| `lheat` | Latent heat | `2264.76e3` | J/kg |
| `rv` | Vapor specific gas constant | `8.3144598/18.01528e-3` | J/(kg K) |
| `ra` | Gas specific gas constant | `8.3144598/28.966e-3` | J/(kg K) |
| `pv` | Vapor pressure at t8 | `f_pvsat(t8)` | Pa |
| `p0` | Initial gas pressure | `Computed by solver` | Pa |
| `pa` | Driving amplitude | `0` | Pa |
| `omega` | Angular frequency | `0` | rad/s |
| `tw` | Pulse width | `0` | s |
| `dt` | Pulse delay | `0` | s |

Other names (model flags, grid sizes, `kappa`, `alphax`, `nstate`, `hugoniot_s`, and viscosity exponents) pass through unchanged. Unknown names are rejected later by the solver parser. Pass-through pairs are appended after converted values; explicit dimensionless overrides therefore take precedence.

## Outputs

- `varin`: dimensionless name/value pairs (see [`default_case`](../model/parameters.md) for their definitions).
- `ref`: reference scales `R0`, `Uc`, `t0`, `P8`, `T8`, `rho8`.

The saturated vapor pressure comes from [`f_pvsat`](../reference/functions/f_pvsat.md), which also supplies the vapor-pressure exponent `thetav = 5200/T8`.
