# `f_mercury_inputs`

`f_mercury_inputs` prepares SI mercury properties and converts them to the dimensionless inputs of pressure-form Keller–Miksis. It works with either IMR solver and requires base MATLAB only.

From the checkout root:

```matlab
addpath('src');
[args, ref] = f_mercury_inputs('t8',298.15,'p8',101325, ...
    'r0',50e-6,'req',40e-6,'mu',0,'surft',0,'kappa',1.4, ...
    'vapor',1,'tvector',linspace(0,12e-6,200));
[t, R] = f_imr_nodal(args{:});
plot(t*ref.t0, R*ref.R0);
```

This example deliberately chooses zero viscosity and zero surface tension. Those are idealized limits, not measured mercury properties. The helper does not invent transport data. Its [model page](../model/mercury.md) explains what the EoS supplies and what remains fixed during the run.

## Inputs and defaults

| Name | Default or requirement | Units / role |
|---|---|---|
| `t8`, `p8` | 298.15 K, 101325 Pa | Ambient absolute state |
| `mu` | Required when `stress > 0` | High-shear viscosity, Pa s, zero allowed |
| `surft` | Required | Surface tension, N/m, zero allowed |
| `kappa` | Required | Noncondensable gas exponent, dimensionless |
| `radial` | 2; all other values rejected | KM pressure only |
| `stress`, `g` | 1, 0 Pa | Newtonian liquid by default; user constitutive choices remain explicit |
| `lambda1`, `lambda2` | Both zero | No default relaxation history |
| `vapor` | 0 | Set 1 to include constant ambient Hg saturation pressure |
| `bubtherm`, `medtherm` | Both zero | Enable only with explicit transport data below |
| `masstrans`, `collapse` | Zero only | Dynamic Hg mass transfer is not implemented |
| `atg`, `btg`, `atv`, `btv` | All required with bubble heat | Gas/vapor conductivity slopes and intercepts in SI |
| `km` | Required with medium heat | Medium conductivity, W/(m K) |
| `dm` | Derived from `km`, Hg density and fitted liquid `Cp`, or supplied | Medium diffusivity, m²/s |
| Geometry, time, forcing, gas `ra` | Same SI interface as `f_dim2nondim` | [Conversion API](f_dim2nondim.md#inputs-si) |
| Grid, integrator and Jacobian settings | Existing solver defaults | [Parameters](../model/parameters.md) |

Thermal runs need positive gas/vapor conductivities at `t8`, a positive `atg` slope for the existing inverse temperature transform, and `kappa > 1`. The converter's documented noncondensable gas constant is air unless `ra` is supplied; it is separate from the Hg-vapor constant. A constant `kappa` is still the existing model approximation, not a new variable-composition heat-capacity model.

The helper derives `rho8`, `C8`, liquid `Cp`, vapor `Rv`, the caloric enthalpy difference and ambient `Pv`. It rejects conflicting inputs such as `rho8`, `c8`, `cstar`, `pv`, `cp`, `rv`, `lheat` or dimensionless overrides. Use SI viscosity, surface tension and modulus instead of `re`, `we` and `ca`. Names are case-insensitive; repeated permitted inputs use the last value.

## Outputs and replay

`args` is a name/value cell array for `f_imr_nodal` or `f_imr_spectral`. For nodal Chebyshev, append `'scheme',1` when calling that solver. The normal [configuration](../running/configuration.md) rules still apply.

`ref` contains the ordinary reference scales plus `ref.mercury` with coefficient provenance, ambient evaluated properties, selected `Pv`, original input pairs and the fixed-ambient assumption. Preserve this structure with the solver output. Rebuild preparation with `f_mercury_inputs(ref.mercury.input_pairs{:})`; the solver's eighth output alone cannot reconstruct the original dimensional material state.

The returned `thetav = 0` is an inactive placeholder because mass transfer is disabled. Do not append options to bypass the helper's guards or treat it as a mercury wall saturation law. Ambient properties also remain fixed if a medium-temperature field is enabled.

## Failure handling

Missing required values raise `IMR:mercuryMissingInput`. Derived/unknown inputs raise `IMR:mercuryOverride`; unsupported radial models raise `IMR:mercuryRadial`. Dynamic mass transfer/collapse raises `IMR:mercuryMassTransfer`. Invalid thermal settings raise `IMR:mercuryThermal`, and ambient boiling raises `IMR:mercuryBoiling`. Property-domain errors follow the [property API](f_mercury_eos.md#errors).

Numerical examples and solver checks await the user's local MATLAB run. See [validation](../validation/mercury.md).
