# Fluid profiles

Material data live in `src/+fluids/`. The package functions `fluids.water()` and `fluids.mercury()` return ordinary MATLAB structs; the shared `f_fluid_properties` evaluator contains the equations. A new fluid using an existing equation family needs one data file, with no changes to the radial equation or an ODE-time file read.

```matlab
addpath('src');
[args, ref] = f_dim2nondim('fluid',@fluids.mercury, ...
    'radial',2,'stress',1,'g',0,'lambda1',0,'lambda2',0, ...
    'mu',0,'surft',0,'kappa',1.4,'req',40e-6, ...
    'tvector',linspace(0,1e-6,41));
[t,R] = f_imr_nodal(args{:});
plot(t*ref.t0,R*ref.R0);
```

Zero viscosity and surface tension above select an idealized demonstration. For a physical run, supply measured values at the selected conditions. `kappa` describes the noncondensable bubble gas and is separate from each phase's heat-capacity ratio.

## Included profiles

| Profile | Liquid properties | Saturation | Solver preparation |
|---|---|---|---|
| `@fluids.water` | Legacy constant reference values | Exponential water fit | Existing options, subject to solver compatibility |
| `@fluids.mercury` | Fitted temperature-dependent Murnaghan slices | Measured table | Pressure-form KM only; ambient vapor optional |

The water profile preserves the converter's historical values, including `rho8 = 1064 kg/m³`. It is a reference set, **not a calibrated temperature-dependent water EoS**. Its constant latent heat and caloric offsets reproduce legacy conversion; its ideal-vapor acoustic ratio is an approximate reference value, not an additional fitted property. Selecting radial 3–6 with water still selects the existing Tait/Hugoniot equations and their independent coefficients. It does not turn the constant-property profile into those equations.

Mercury coefficients and frozen reference data are unchanged. [Its model page](mercury.md) explains the source domains and why an acoustically consistent enthalpy extension remains deferred. The generic interface preserves those scientific restrictions.

## Data schema

Copy an included profile to `src/+fluids/myfluid.m`, rename its function to `myfluid`, and edit all relevant coefficients, provenance and limits. Select it with `'fluid',@fluids.myfluid`. Profiles contain data, with no fitted calculations or transport fallbacks hidden in the evaluator.

| Field | Purpose |
|---|---|
| `name`, `source`, `commit`, `units` | Identity, provenance and SI convention |
| `temperature_range`, `pressure_range` | Accepted bulk-property domain |
| `sound_temperature_range` | Measured acoustic interval; warns outside |
| `liquid.model` | `constant` or `murnaghan` |
| `liquid` coefficients | Constant `rho`, `c`, or the Murnaghan coefficients; both require `Cp`, `q` |
| `vapor.model`, `Rs`, `gamma` | `ideal`, specific gas constant and acoustic heat-capacity ratio |
| `caloric.model` | `linear_enthalpy` (both phases' `Cp`, `q`) or `constant_latent_heat` (`caloric.latent_heat`) |
| `saturation.model` | `exp_inverse_temperature` (`c1`, `c2`) or `table_logp_inverseT` (ascending row vectors `T`, `P`) |
| `saturation.temperature_range` | Saturation domain, independent of bulk-property domain |
| `defaults` | SI transport defaults; gas, stress and forcing choices belong in the simulation call |
| `derived_inputs` | SI inputs locked to evaluated EoS values; conflicting user pairs are rejected |
| `required`, `required_viscous`, `required_thermal`, `required_medium` | Inputs the user must supply, unconditionally or for the selected model |
| `capabilities.radial`, `mass_transfer`, `collapse` | Supported radial selectors and phase-change capabilities |

Set conservative capabilities backed by your model and data. For a table, use positive, strictly increasing temperature nodes, positive pressures and matching row-vector lengths. Its accepted range must stay within the nodes. Dynamic mass transfer requires the exponential law used by the current wall closure; setting a capability flag cannot create that law. General solver restrictions still apply.

A new mathematical EoS family needs an additional model branch and tests in `f_fluid_properties`. The present Murnaghan family follows the mercury source's coefficient convention; it is not a universal parameterization for every liquid.

## Migration from the mercury-specific API

| Previous call | Replacement |
|---|---|
| `f_mercury_parameters()` | `fluids.mercury()` |
| `f_mercury_eos(P,T)` | `f_fluid_properties(@fluids.mercury,P,T)` |
| `f_mercury_pvsat(T)` | `sat = f_fluid_properties(@fluids.mercury,[],T); Pv = sat.pressure;` |
| `f_mercury_inputs(...)` | `f_dim2nondim('fluid',@fluids.mercury,...)` |
| `[Pv,~,c2] = f_pvsat(T)` | `sat = f_fluid_properties(@fluids.water,[],T);` then `sat.pressure`, `sat.c2` |

The four `f_mercury_*` functions and `f_pvsat` are removed. `ref.fluid` replaces `ref.mercury`, with the resolved profile in `ref.fluid.profile`. Calls to `f_dim2nondim` without `fluid` retain the legacy conversion and dimensionless-override behavior. Explicit profiles reject dimensionless physical overrides so an EoS-derived value cannot be silently replaced after conversion.

The shared converter retains its ordinary stress and elasticity defaults for every fluid. To reproduce the old mercury helper's Newtonian setup, explicitly select `stress = 1`, `g = 0` and `lambda1 = lambda2 = 0`, as above. The radial selector defaults to the first supported branch if the ordinary `radial = 1` is unavailable; mercury therefore selects radial 2. A profile does not supply a noncondensable gas or a material stress law.
