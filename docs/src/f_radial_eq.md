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
   - `f_mie_gruneisen_eos_scalar`: sound speed `C`, enthalpy `hB` and `hH = 1/rho` at the wall pressure
   - `f_mie_mu`: density strain `mu = rho - 1` at a pressure

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
- Mie–Grüneisen enthalpy: with `s = hugoniot_s` and `A = p/Cstar^2 = g(mu) = mu*(1 + nog*mu)/(1 - s*mu)^2`, `h(P) = int_1^P dp/(1 + mu)` is integrated by parts, `h = [p/(1 + mu)] + Cstar^2*int g/(1 + mu)^2 dmu`, and the partial fractions of `g/(1 + mu)^2` (poles at `mu = -1` and `mu = 1/s`) give it in closed form. The differences are written so that they hold for `s = 0` and near `P = 1`. Against an adaptive `integral` at `1e-12` the closed form agrees to `2.4e-13` relative to `|P - 1|` over 396 states (`Cstar` 5–1000, `nog` 0.5–2.5, `s` 0–2, `P` from -0.5 to `1e4`). It replaces an adaptive `integral` on every right-hand-side call, and radial 5 and 6 runs take about half the time, with `R(t)` unchanged to `1e-13`.
- `f_checkers` allows radial 5 and 6 only with the polytropic bubble, and radial 3 and 4 without medium heat transfer: the enthalpy forms take the liquid at the far-field temperature.
