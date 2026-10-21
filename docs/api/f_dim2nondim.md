# `f_dim2nondim`

`f_dim2nondim` is a side function that converts dimensional (SI) physical properties into the dimensionless inputs taken by `f_imr_nodal` and `f_imr_spectral`. The solvers never call it: they accept only dimensionless inputs.

See the [complete quick start](../getting-started/quick-start.md) for a runnable SI-to-solver example.

## Fluid profiles

Select a MATLAB property file with `'fluid',@fluids.water` or `'fluid',@fluids.mercury`. A resolved profile struct is also accepted. Profile functions run once during preparation; the solvers receive dimensionless numbers, not material callbacks.

```matlab
[args, ref] = f_dim2nondim('fluid',@fluids.mercury, ...
    'radial',2,'stress',1,'g',0,'lambda1',0,'lambda2',0, ...
    't8',298.15,'p8',101325,'mu',0,'surft',0,'kappa',1.4, ...
    'r0',50e-6,'req',40e-6,'tvector',linspace(0,1e-6,41));
[t,R] = f_imr_nodal(args{:});
```

Explicit SI transport inputs override a profile's defaults. EoS-derived density, acoustics and caloric inputs are protected by its `derived_inputs` list. Dimensionless physical overrides such as `cstar` or `re` are rejected for explicit profiles. Model flags, numerical options, stress parameters and forcing remain ordinary inputs. Without `fluid`, the existing defaults and pass-through behavior below are retained.

The mercury profile selects its only supported radial model, pressure-form KM (`radial = 2`). Stress, elasticity and relaxation retain the converter's ordinary defaults; select the Newtonian case explicitly as in the example. It requires explicit `surft`, bubble-gas `kappa`, and `mu` when `stress > 0`. With `bubtherm = 1`, supply `atg`, `btg`, `atv`, `btv`; with `medtherm = 1`, also supply positive `km`. Thermal calculations need positive conductivities, positive `atg` and `kappa > 1`. Inactive placeholder transport scales in its profile are not physical mercury measurements.

Mercury `vapor = 1` selects ambient saturation pressure from its table; otherwise `pv = 0`. Dynamic mass transfer, collapse and radial 3–6 are rejected. The [mercury model](../model/mercury.md) explains these limits. Water retains its legacy exponential wall law and converter defaults, subject to [solver compatibility](../reference/compatibility.md).

Profile errors use `IMR:fluidMissingInput`, `IMR:fluidOverride`, `IMR:fluidRadial`, `IMR:fluidMassTransfer`, `IMR:fluidThermal` or `IMR:fluidBoiling`; invalid numeric inputs can also raise MATLAB validation errors. See [property errors](f_fluid_properties.md#errors) for EoS domain failures.

See [fluid profiles](../model/fluids.md) for data-file fields and migration from the removed mercury helpers.

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
| `pv` | Vapor pressure at t8 | Water profile saturation | Pa |
| `p0` | Initial gas pressure | `Computed by solver` | Pa |
| `pa` | Driving amplitude | `0` | Pa |
| `omega` | Angular frequency | `0` | rad/s |
| `tw` | Pulse width | `0` | s |
| `dt` | Pulse delay | `0` | s |

Without a profile, other names (model flags, grid sizes, `kappa`, `alphax`, `nstate`, `hugoniot_s`, and viscosity exponents) pass through unchanged. Unknown names are rejected later by the solver parser. Pass-through pairs are appended after converted values; explicit dimensionless overrides therefore take precedence in this legacy mode.

## Outputs

- `varin`: dimensionless name/value pairs (see [`default_case`](../model/parameters.md) for their definitions).
- `ref`: reference scales `R0`, `Uc`, `t0`, `P8`, `T8`, `rho8`.

Saturation is evaluated by [`f_fluid_properties`](f_fluid_properties.md#saturation). For the exponential model, `thetav = -c2/T8`; a saturation table sets `thetav = 0` and cannot enable the exponential wall closure.

With an explicit profile, `ref.fluid` records the resolved `profile`, unmodified evaluated `properties`, effective SI settings in `SI`, selected ambient `Pv`, original `input_pairs` with the profile snapshot, and the pressure-form KM assumption. For constant-property profiles, an SI override can differ from the original `properties`; `SI` and the top-level reference scales record the values actually used. Preserve `ref` alongside solver outputs. Replay with:

```matlab
[args, ref] = f_dim2nondim(ref.fluid.input_pairs{:});
```

The snapshot avoids depending on later edits to the profile file. The solver's eighth output alone does not reconstruct the dimensional material state.
