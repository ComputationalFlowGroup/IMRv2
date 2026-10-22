# `imr.physics.radial_acceleration`

`imr.physics.radial_acceleration` evaluates the bubble-wall acceleration `Rddot` for the six radial models. Rayleigh–Plesset treats the liquid as incompressible, the pressure form of Keller–Miksis uses a constant sound speed, and the other four use an equation of state: Tait for radial 3 and 4, Mie–Grüneisen for radial 5 and 6.

## Major Sections

1. **Rayleigh–Plesset (Model 1)**, incompressible
2. **Keller–Miksis Pressure Form (Model 2)**, constant sound speed
3. **Keller–Miksis Enthalpy Form with Tait EoS (Model 3)**
4. **Gilmore with Tait EoS (Model 4)**
5. **Keller–Miksis with Mie–Grüneisen EoS (Model 5)**
6. **Gilmore with Mie–Grüneisen EoS (Model 6)**
7. **EOS Helper Functions**:
    - `f_mie_gruneisen_eos_scalar`: sound speed `C`, enthalpy `hB` and `hH = 1/rho` at the wall or far-field pressure
    - `f_mie_mu`: density strain `mu = rho - 1` at a pressure

## Inputs

| Name        | Description                                | Units |
|-------------|--------------------------------------------|-------|
| `radial`    | Model selector (1–6)                      | –     |
| `P`         | Internal bubble pressure (non-dim)        | –     |
| `Pdot`      | Time derivative of `P`                    | –     |
| `Pf8`       | Far-field pressure above ambient (non-dim)               | –     |
| `Pf8dot`    | Time derivative of `Pf8`                  | –     |
| `iWe`       | Inverse Weber number                      | –     |
| `R`         | Bubble radius (non-dim)                   | –     |
| `Rdot`      | Bubble wall velocity (non-dim)            | –     |
| `S`         | Stress term (non-dim)                     | –     |
| `Sdot`      | Acceleration-independent part of the stress derivative                    | –     |
| `Cstar`     | Non-dim sound speed                       | –     |
| `sam`       | `1 + GAMa`, Tait pressure at rest         | –     |
| `no`        | `(nstate - 1)/nstate`                     | –     |
| `GAMa`      | Tait parameter `GAM/P8`                   | –     |
| `nstate`    | Tait exponent                             | –     |
| `nog`       | `(nstate - 1)/2`, used by the Mie–Grüneisen EoS | –     |
| `hugoniot_s`| Hugoniot slope                            | –     |
| `JdotA`     | `Rddot` coefficient of the Newtonian viscous stress, moved to the left-hand side | –     |
| `ddintfnu`  | `Rddot` coefficient of the time derivative of the non-Newtonian stress integral | –     |
| `iDRe`      | Common non-Newtonian factor; see [definition](../../model/stress-viscosity.md#non-newtonian-viscosity) | –     |

## Outputs

| Name     | Description                             | Units |
|----------|-----------------------------------------|-------|
| `Rddot`  | Bubble wall acceleration                | –     |

## Notes

- The EOS helper functions work elementwise, on one pressure or a row of them.
- `R`, `Rdot` and the other state inputs may be row vectors, with `X` holding one state per column (the batched right-hand side of [`imr.solve_nodal`](../../method/time-integration.md#batched-right-hand-side)). Operations are elementwise across the columns. Driven enthalpy branches now evaluate far-field enthalpy and its pressure derivative with the same EoS as the wall; see the [equations](../../theory/radial-equations.md#enthalpy-forms).
- `imr.configuration.validate_case` allows radial 5 and 6 only with the polytropic bubble, and radial 3 and 4 without medium heat transfer: the enthalpy forms take the liquid at the far-field temperature.
