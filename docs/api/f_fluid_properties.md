# `imr.fluid_properties`

One MATLAB evaluator serves every fluid profile. Add `src/` to the MATLAB path; no installation of the fitting repository or external runtime is needed.

```matlab
imr_setup;
p = imr.fluids.water();                 % data only, all coefficients in SI
state = imr.fluid_properties(p,101325,298.15);
state.liquid.rho                      % kg/m^3
state.liquid.c                        % m/s
sat = imr.fluid_properties(p,[],298.15); % [] explicitly selects saturation
sat.pressure                          % Pa
sat.dpressure_dT                      % Pa/K
```

Pass a profile struct or a zero-argument function handle such as `@imr.fluids.water`. The evaluator dispatches on model identifiers in that data, rather than material names. See [fluid profiles](../model/fluids.md) for the schema and [unit conversion](f_dim2nondim.md#fluid-profiles) for solver preparation.

## Bulk properties

`state = imr.fluid_properties(profile,P,T)` accepts absolute pressure in Pa and temperature in K. Scalars and compatible numeric arrays expand to the same output shape. `T` must be positive; both inputs must be finite and real and within the profile domain. At zero pressure ideal-vapor density is zero and specific volume is `Inf`.

| Field | Meaning | Units |
|---|---|---|
| `P`, `T` | Expanded inputs | Pa, K |
| `liquid.rho`, `liquid.c` | Liquid density and reported sound speed | kg/m³, m/s |
| `liquid.h` | Profile's liquid caloric slice | J/kg |
| `vapor.rho`, `vapor.v` | Ideal-vapor density and specific volume | kg/m³, m³/kg |
| `vapor.c`, `vapor.h` | Ideal-vapor sound speed and caloric slice | m/s, J/kg |
| `latent_heat` | Vapor minus liquid enthalpy on the profile convention | J/kg |

The Murnaghan family also reports `liquid.rho0`, `c0`, `K0` and `K0_prime` (kg/m³, m/s, Pa and dimensionless). These are temperature-specific zero-pressure values. See the [Murnaghan coefficient convention](../model/fluids.md#murnaghan-coefficient-convention).

## Saturation

`sat = imr.fluid_properties(profile,[],T)` evaluates **only saturation**, using its separate temperature range. It returns `model`, `pressure` (Pa), `dpressure_dT` (Pa/K) and `c2` (K). The exponential model returns its negative inverse-temperature coefficient as `c2`; a table returns zero because it does not supply the exponential wall law.

For the exponential model,

$$
P_{\mathrm{sat}}(T)=c_1\exp(c_2/T),\qquad
\frac{dP_{\mathrm{sat}}}{dT}=-\frac{c_2}{T^2}P_{\mathrm{sat}}(T).
$$

Here $T$ is dimensional temperature, $c_1$ is in Pa and $c_2$ is in K. The legacy water profile uses $c_1=1.17\times10^{11}\,\mathrm{Pa}$ and $c_2=-5200\,\mathrm{K}$. The fit is attributed to Preston (2004); an unrestricted positive input domain is not an accuracy guarantee.

For the table model, $\log P_{\mathrm{sat}}$ is linear in $1/T$ between nodes. Values match the data nodes, and extrapolation is rejected. Its derivative is piecewise smooth: at an interior node the evaluator returns the derivative from higher temperature, and at the highest temperature it uses the derivative from below. Table interpolation supports ambient vapor pressure only; it cannot be substituted for the analytically inverted wall law.

## Errors

Domain failures use `IMR:fluidPressure`, `IMR:fluidTemperature`, `IMR:fluidSaturation` or `IMR:fluidDomain`. Unknown model identifiers use `IMR:fluidModel`. Sound speed outside a profile's measured temperature interval warns with `IMR:fluidExtrapolation`; that interval can be narrower than the accepted property domain. Invalid numeric types/shapes can also raise MATLAB validation or array-expansion errors.
