# `imr.nondimensionalize`

`imr.nondimensionalize` is a side function. It converts dimensional (SI) physical properties into the dimensionless inputs that `imr.solve_nodal` and `imr.solve_spectral` take. The solvers never call it. The solvers accept only dimensionless inputs.

See the [complete quick start](../getting-started/quick-start.md) for a runnable SI-to-solver example.

## Fluid profiles

Select `'fluid',@imr.fluids.water`, a new MATLAB profile function, or a resolved profile struct. The code evaluates profile data once during preparation. The solver receives dimensionless numbers, not a callback.

```matlab
imr_setup;
[args, ref] = imr.nondimensionalize('fluid',@imr.fluids.water, ...
    'radial',2,'stress',1,'g',0,'lambda1',0,'lambda2',0, ...
    't8',298.15,'p8',101325,'mu',1e-3,'surft',0.072,'kappa',1.4, ...
    'r0',50e-6,'req',40e-6,'tvector',linspace(0,1e-6,41));
[t,R] = imr.solve_nodal(args{:});
```

Explicit SI transport inputs override reference values in `profile.transport`. `profile.derived_inputs` lists protected EoS-derived inputs. The code rejects a conflicting override, and a dimensionless physical override such as `cstar` or `re`. Model, grid, numerical and forcing choices remain ordinary inputs. A call without `fluid` preserves legacy conversion and pass-through behavior.

Required active-model inputs and supported radial/phase-change choices come from the profile's `required_*` and `capabilities` fields. Gas conductivity, stress and forcing are simulation choices. They do not depend on liquid/vapor material data. The converter uses internal neutral scales only for inactive terms that have no transport data. Those scales are not physical properties. Those scales cannot satisfy an active model's explicit input requirements.

Profile errors use `IMR:fluidMissingInput`, `IMR:fluidOverride`, `IMR:fluidRadial`, `IMR:fluidMassTransfer`, `IMR:fluidThermal` or `IMR:fluidBoiling`. Invalid numeric inputs can also raise MATLAB validation errors. See [property errors](fluid_properties.md#errors), [fluid profiles](../model/fluids.md) and [defaults and units](defaults.md).

## Inputs (SI)

These defaults come from `src/+imr/nondimensionalize.m`. Give inputs as name/value pairs. If you repeat a dimensional name, the code uses the last value. An explicit `du` overrides the viscosity difference derived from `mu0` and `mu`.

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

Without a profile, other names (model flags, grid sizes, `kappa`, `alphax`, `nstate`, `hugoniot_s`, and viscosity exponents) pass through unchanged. The solver parser rejects an unknown name later. The code appends pass-through pairs after converted values. So, in this legacy mode, an explicit dimensionless override takes precedence.

## Outputs

- `varin`: dimensionless name/value pairs (see [`imr.configuration.defaults`](../model/parameters.md) for their definitions).
- `ref`: reference scales `R0`, `Uc`, `t0`, `P8`, `T8`, `rho8`.

[`imr.fluid_properties`](fluid_properties.md#saturation) evaluates saturation. For the exponential model, `thetav = -c2/T8`. A saturation table sets `thetav = 0`, and cannot enable the exponential wall closure.

With an explicit profile, `ref.fluid` records the resolved `profile`, the unmodified evaluated `properties`, the effective SI settings in `SI`, the selected ambient `Pv`, the original `input_pairs` with the profile snapshot, and the pressure-form KM assumption. For a constant-property profile, an SI override can differ from the original `properties`. `SI` and the top-level reference scales record the values actually used. Keep `ref` alongside solver outputs. Replay with:

```matlab
[args, ref] = imr.nondimensionalize(ref.fluid.input_pairs{:});
```

The snapshot avoids depending on later edits to the profile file. The solver's eighth output alone does not reconstruct the dimensional material state.
