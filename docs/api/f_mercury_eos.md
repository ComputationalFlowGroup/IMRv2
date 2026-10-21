# Mercury property API

These base-MATLAB functions expose the versioned liquid and vapor mercury fit. They do not read external files or call another language.

## Property evaluation

From the checkout root in MATLAB:

```matlab
addpath('src');
hg = f_mercury_eos(101325,298.15);
rho = hg.liquid.rho;
c = hg.liquid.c;
```

| Input | Meaning | Units |
|---|---|---|
| `P` | Absolute liquid/ideal-vapor evaluation pressure | Pa |
| `T` | Evaluation temperature | K |

Both arguments can be scalars or arrays compatible with MATLAB implicit expansion. All output arrays share the expanded shape. For example, a pressure row and temperature column produce a full grid. Equations and domain limits live in [the mercury model](../model/mercury.md).

| Output field | Meaning | Units |
|---|---|---|
| `P`, `T` | Expanded input arrays | Pa, K |
| `liquid.rho`, `liquid.c`, `liquid.h` | Liquid density, sound speed, caloric enthalpy | kg/m³, m/s, J/kg |
| `liquid.rho0`, `liquid.c0` | Zero-reference-pressure density and sound speed | kg/m³, m/s |
| `liquid.K0`, `liquid.K0_prime` | Temperature-specific reference modulus and pressure derivative | Pa, dimensionless |
| `vapor.rho`, `vapor.v` | Ideal-vapor density and specific volume | kg/m³, m³/kg |
| `vapor.c`, `vapor.h` | Ideal-vapor sound speed and enthalpy | m/s, J/kg |
| `latent_heat` | Difference of vapor and liquid caloric enthalpies | J/kg |

At `P = 0`, ideal-vapor density is zero and its specific volume is `Inf`; the liquid outputs remain finite. That limiting ideal-gas result is deliberate.

## Coefficient metadata

`p = f_mercury_parameters()` returns the complete fitted coefficients, diagnostic residuals, source revision and measured saturation nodes. Field meanings and units are listed in [the model page](../model/mercury.md#fitted-liquid-slices). Coefficients are bundled MATLAB literals; no fitting or optimization toolbox is used.

## Ambient saturation pressure

```matlab
addpath('src');
Pv = f_mercury_pvsat([298.15 350 400]);
```

The input is temperature in K, and the output is pressure in Pa with the same shape. The [model page](../model/mercury.md#vapor-and-saturation) defines the interpolation, data support and limitations.

## Errors

| Identifier | Meaning |
|---|---|
| `IMR:mercuryPressure` | Unsupported pressure, including tension |
| `IMR:mercuryTemperature` | Temperature outside the nominal liquid fit window |
| `IMR:mercuryDomain` | Nonpositive Murnaghan base |
| `IMR:mercurySaturation` | Temperature outside measured saturation support |
| `IMR:mercuryExtrapolation` | Warning: acoustic temperature extrapolation |

MATLAB input validation also rejects empty, nonnumeric, complex or nonfinite arguments. An incompatible array shape raises MATLAB's size error. The examples are source-checked and await local MATLAB execution; they do not show invented simulation results.
