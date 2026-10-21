# Mercury liquid and vapor properties

The MATLAB profile `fluids.mercury()` supplies fitted SI liquid and ideal-vapor coefficients to the shared `f_fluid_properties` evaluator. Selecting it in `f_dim2nondim` uses ambient liquid density and sound speed to prepare the pressure-form Keller–Miksis equation (`radial = 2`) for either solver.

## What enters the radial equation

At ambient absolute pressure `P8` and temperature `T8`, the converter evaluates `rho8 = rho(P8,T8)` and `C8 = c(P8,T8)`. `f_dim2nondim` then computes consistent velocity, time and acoustic scales. The standard KM pressure equation holds these ambient liquid properties fixed during a run. It does not evaluate pressure-dependent density or sound speed at the moving bubble wall.

There is no Python implementation or runtime call in IMRv2. Coefficients, saturation data, examples and tests are MATLAB files bundled with this repository. The original fitting repository is a provenance reference, not an installation dependency.

## Preparing a run

Run `imr_setup` once. Use `f_dim2nondim('fluid',@fluids.mercury,...)`; the profile selects its only supported radial model, pressure-form KM. Stress and gas choices remain explicit:

```matlab
[args, ref] = f_dim2nondim('fluid',@fluids.mercury, ...
    'radial',2,'stress',1,'g',0,'lambda1',0,'lambda2',0, ...
    'mu',0,'surft',0,'kappa',1.4,'vapor',1, ...
    'req',40e-6,'tvector',linspace(0,1e-6,41));
[t,R] = f_imr_nodal(args{:});
```

Zero viscosity and surface tension above select an idealized demonstration. For physical simulations provide measured values at the selected conditions. `surft` and bubble-gas `kappa` are required; `mu` is required when `stress > 0`. Thermal runs require explicit `atg`, `btg`, `atv`, `btv`, positive conductivities and gas slope, and `kappa > 1`; medium heat additionally requires positive `km`. The profile intentionally contains no fabricated transport values.

For migration, `fluids.mercury()` replaces `f_mercury_parameters`, `f_fluid_properties` replaces `f_mercury_eos` and `f_mercury_pvsat`, and the SI converter replaces `f_mercury_inputs`. Use the Newtonian settings in the example to reproduce that helper's model choices. `ref.fluid.profile` and `ref.fluid.input_pairs` preserve a data snapshot for replay.

## Fitted liquid slices

`f_fluid_properties(@fluids.mercury,P,T)` evaluates

$$
\rho_0=a_\rho T+b_\rho,\quad c_0=c_{0,0}+c_{0,1}T,\quad
K_0=\rho_0c_0^2/\gamma,\quad n=1+2K_0s_1/c_0,
$$

$$
\rho(P,T)=\rho_0(1+nP/K_0)^{1/n},\qquad
c(P,T)=\sqrt{\gamma(K_0+nP)/\rho(P,T)}.
$$

Here $\rho_0$ and $c_0$ are the fitted zero-pressure density and sound speed at temperature $T$; $K_0$ is the corresponding reference bulk modulus and $n$ the Murnaghan pressure exponent. The dimensionless base $1+nP/K_0$ accounts for compression from zero pressure, using absolute pressure $P$ in Pa and $T$ in K. At $P=0$ the equations recover $\rho_0$ and $c_0$. These source conventions must be preserved together.

The liquid caloric slice is `h_l = Cp_l*T + q_l`. `K0` and `n` are recomputed at each temperature; the single tabulated pair applies only at `Tref = 412.5 K`. The quadratic diagnostic coefficient `s2` is retained as metadata and is not part of these evaluators.

| Parameter group in `fluids.mercury()` | Units |
|---|---|
| `liquid.rho0_a`, `rho0_b` | kg/(m³ K), kg/m³ |
| `liquid.c0_0`, `c0_1` | m/s, m/(s K) |
| `liquid.s1`, `s2` | m/(s Pa), m/(s Pa²) |
| `liquid.K0_at_Tref`, `K0_prime_at_Tref`, `Tref` | Pa, dimensionless, K |
| Both phases' `Cp`, `Cv` | J/(kg K) |
| Both phases' `q` | J/kg |
| Both phases' `gamma`, `Gamma` | Dimensionless |
| `vapor.Rs` | J/(kg K) |
| `saturation.T`, `saturation.P` | K, Pa |
| Fields ending in `_percent` | Percent RMS diagnostic |

`gamma = Cp/Cv` is distinct from bubble `kappa` and the Grüneisen coefficient `Gamma`. The source's liquid `Gamma = gamma - 1` is provisional and is not used by the evaluators. These slices do not implement a complete `p(rho,e)` liquid model.

## Vapor and saturation

The ideal-vapor slices are `v = Rs*T/P`, `rho = P/(Rs*T)`, `c = sqrt(gamma*Rs*T)` and `h_g = Cp_g*T + q_g`. `state.latent_heat` reports `h_g - h_l` on the source's shared caloric convention. It is available for inspection, not used to enable a dynamic phase-change calculation.

[The shared saturation evaluator](../api/f_fluid_properties.md#saturation) interpolates `log(Psat)` linearly in `1/T` between measured vapor-table nodes. This matches those nodes and gives a positive, monotone ambient interpolation. It is not a new fitted saturation law or proof of interpolation accuracy. `vapor = 1` adds that constant ambient saturation pressure to a polytropic noncondensable gas bubble. Bubble `kappa` remains an explicit gas input.

Dynamic mercury mass transfer and collapse initialization are rejected. Thermal runs require explicit conductivities because the EoS fit supplies none. See the [SI converter](../api/f_dim2nondim.md) for required transport inputs.

## Limits and enthalpy formulations

The property evaluator accepts 275–550 K and 0–6.23 GPa. Sound-speed measurements cover 296–513 K in disjoint pressure windows; outside that temperature interval it warns about extrapolation. Being inside the rectangular accepted domain does not guarantee local data coverage. Negative pressure, nonfinite inputs and values outside the accepted domain are rejected. Saturation interpolation has the narrower 298.15–550 K interval and rejects extrapolation.

The source reports 0.71249% liquid sound-speed RMS over its full dataset, but about 2% uncertainty on the saturation sound-speed slice. The property port preserves the fit rather than improving it. No mercury viscosity, surface tension, conductivity or diffusivity is inferred from these coefficients.

At fixed temperature the source density gives `dP/drho = (K0+nP)/rho`, while reported `c²` is `gamma` times that value. Treating both as one barotrope would make them inconsistent. Consequently this profile rejects radial 3–6; the existing Hugoniot-based MG models remain separate. Neither the caloric `h_l(T)` nor the tabulated reference modulus can simply replace the pressure integral needed by an enthalpy model. An acoustically consistent enthalpy extension remains deferred.

## Provenance

Coefficients and observations come from [mrodrig6/mercury_eos](https://github.com/mrodrig6/mercury_eos), commit `5d51892ee4eef7d898d33ac81a02bee52d554165`, using its canonical Murnaghan evaluator. Preserve that revision when changing coefficients.

- Vukalovich et al., *Thermophysical Properties of Mercury*, Table 12 (1971): liquid saturation density and enthalpy.
- Ayrinhac et al., *J. Chem. Phys.* 140, 244201 (2014), [doi:10.1063/1.4882695](https://doi.org/10.1063/1.4882695): 115 sound-speed observations.
- Fokin, Popov and Naurzakov, *High Temperature* 49, 832–840 (2011), [doi:10.1134/S0018151X11050075](https://doi.org/10.1134/S0018151X11050075): vapor properties and saturation pressures.
- Holman and ten Seldam, *J. Phys. Chem. Ref. Data* 23, 807–827 (1994), [doi:10.1063/1.555952](https://doi.org/10.1063/1.555952): liquid heat-capacity anchor.

See [mercury validation](../validation/mercury.md) for bundled reference values, current execution status and local MATLAB commands.
