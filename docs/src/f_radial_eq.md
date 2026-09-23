# `f_radial_eq.m` Documentation

## Overview
`f_radial_eq` evaluates the bubble-wall acceleration `Rddot` for the six radial models. Rayleigh–Plesset treats the liquid as incompressible, the pressure form of Keller–Miksis uses a constant sound speed, and the other four use an equation of state: Tait for radial 3 and 4, Mie–Grüneisen for radial 5 and 6.

## Major Sections
1. **Rayleigh–Plesset (Model 1)**, incompressible
2. **Keller–Miksis Pressure Form (Model 2)**, constant sound speed
3. **Keller–Miksis Enthalpy Form with Tait EoS (Model 3)**
4. **Gilmore with Tait EoS (Model 4)**
5. **Keller–Miksis with Mie–Grüneisen EoS (Model 5)**
6. **Gilmore with Mie–Grüneisen EoS (Model 6)**
7. **EOS Helper Functions**:
   - `f_mie_gruneisen_eos_scalar`
   - `f_mie_rho_from_p_scalar`

## Inputs

| Name        | Description                                | Units |
|-------------|--------------------------------------------|-------|
| `radial`    | Model selector (1–6)                      | –     |
| `P`         | Internal bubble pressure (non-dim)        | –     |
| `Pdot`      | Time derivative of `P`                    | –     |
| `Pf8`       | External pressure (non-dim)               | –     |
| `Pf8dot`    | Time derivative of `Pf8`                  | –     |
| `iWe`       | Inverse Weber number                      | –     |
| `R`         | Bubble radius (non-dim)                   | –     |
| `Rdot`      | Bubble wall velocity (non-dim)            | –     |
| `S`         | Stress term (non-dim)                     | –     |
| `Sdot`      | Time derivative of `S`                    | –     |
| `Cstar`     | Non-dim sound speed                       | –     |
| `sam`       | `1 + GAMa`, Tait pressure at rest         | –     |
| `no`        | `(nstate - 1)/nstate`                     | –     |
| `GAMa`      | Tait parameter `GAM/P8`                   | –     |
| `nstate`    | Tait exponent                             | –     |
| `nog`       | `(nstate - 1)/2`, used by the Mie–Grüneisen EoS | –     |
| `hugoniot_s`| Hugoniot slope                            | –     |
| `JdotA`     | `Rddot` coefficient of the Newtonian viscous stress, moved to the left-hand side | –     |
| `ddintfnu`  | `Rddot` coefficient of the time derivative of the non-Newtonian stress integral | –     |
| `iDRe`      | `1/DRe`, the inverse non-Newtonian Reynolds number | –     |

## Outputs

| Name     | Description                             | Units |
|----------|-----------------------------------------|-------|
| `Rddot`  | Bubble wall acceleration                | –     |

## Notes
- EOS helper functions perform scalar EOS and density calculations.
- `f_checkers` allows radial 5 and 6 only with the polytropic bubble.
