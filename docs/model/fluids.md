# Fluid profiles

Material data live in MATLAB's `src/+fluids/` package. `fluids.water()` returns a struct of SI liquid/vapor coefficients and reference transport values; `f_fluid_properties` evaluates the selected mathematical model. The simulation call chooses the bubble gas, stress model and forcing.

```matlab
imr_setup;
[args, ref] = f_dim2nondim('fluid',@fluids.water, ...
    'radial',2,'stress',1,'g',0,'lambda1',0,'lambda2',0, ...
    'mu',1e-3,'surft',0.072,'kappa',1.4,'req',40e-6, ...
    'tvector',linspace(0,1e-6,41));
[t,R] = f_imr_nodal(args{:});
plot(t*ref.t0,R*ref.R0);
```

The water profile preserves the legacy reference set, including `rho8 = 1064 kg/m³`. It is **not a calibrated temperature-dependent water EoS**. Constant latent heat and caloric offsets reproduce legacy conversion; its ideal-vapor acoustic ratio is an approximate diagnostic. Selecting radial 3–6 still selects the separate Tait/Hugoniot equations and their coefficients.

## Data schema

Copy `water.m` to `src/+fluids/myfluid.m`, rename the function to `myfluid`, and edit the relevant data, provenance and capabilities. Select it with `'fluid',@fluids.myfluid`. Existing equation families need no new solver file or ODE-time file read.

| Field | Purpose |
|---|---|
| `name`, `source`, `commit` | Material identity and reproducible provenance |
| `temperature_range`, `pressure_range` | Accepted bulk-property domain |
| `sound_temperature_range` | Measured acoustic interval; warn outside |
| `liquid.model` | `constant` or `murnaghan` |
| `liquid` coefficients | Constant `rho`, `c`, or Murnaghan coefficients; both require `Cp`, `q` |
| `liquid.tait_pressure` | Optional SI Tait pressure shift; separate from the bulk evaluator |
| `vapor.model`, `Rs`, `gamma` | `ideal`, specific gas constant and acoustic heat-capacity ratio |
| `caloric.model` | `linear_enthalpy` (both phases' `Cp`, `q`) or `constant_latent_heat` (`caloric.latent_heat`) |
| `saturation.model` | `exp_inverse_temperature` (`c1`, `c2`) or `table_logp_inverseT` (ascending row vectors `T`, `P`) |
| `saturation.temperature_range` | Saturation domain, independent of the bulk-property domain |
| `transport` | Available SI `surft`, `mu`, `km`, `dmass`, `atv`, `btv`; each has a unit beside its assignment |
| `derived_inputs` | SI inputs protected because the EoS determines them |
| `required`, `required_viscous`, `required_thermal`, `required_medium` | Inputs callers must explicitly provide for active models |
| `capabilities.radial`, `mass_transfer`, `collapse` | Allowed radial models and phase-change capabilities |

An empty transport struct means the source supplies no such data. Missing values must not be copied from another material or disguised as profile defaults. The former `p.defaults` and unused `p.units` fields have been removed; [defaults and units](../reference/defaults.md) explains where their responsibilities went. SI units are the fixed API contract, not a string selected by the profile.

Set capabilities conservatively. Table saturation needs positive temperatures and pressures, matching row-vector lengths and strictly increasing temperatures. Its accepted domain must stay within the nodes. Dynamic mass transfer requires the exponential law analytically inverted by the wall closure; changing a capability flag cannot create that law. General solver compatibility also applies.

## Murnaghan coefficient convention

The shared evaluator supports a source convention with temperature-dependent zero-pressure density and sound speed:

$$
\rho_0=a_\rho T+b_\rho,\quad c_0=c_{0,0}+c_{0,1}T,\quad
K_0=\rho_0c_0^2/\gamma,\quad n=1+2K_0s_1/c_0.
$$

It then evaluates $\rho=\rho_0(1+nP/K_0)^{1/n}$ and $c^2=\gamma(K_0+nP)/\rho$. Profile coefficient names are `rho0_a`, `rho0_b`, `c0_0`, `c0_1`, `gamma` and `s1`. This is a particular property-family convention, not a universal calibration for all liquids. Density and reported acoustics must be checked for thermodynamic consistency before coupling a new profile to enthalpy radial models.

A genuinely new equation family needs a model branch and tests in `f_fluid_properties`. A new profile using an existing family needs only data and appropriate validation evidence.

## Optional material add-on

The full distribution includes a [mercury profile](mercury.md) with its own measured-data provenance, input requirements and restrictions. Its profile, example, test suite and fixtures are independent of the water reference path.
