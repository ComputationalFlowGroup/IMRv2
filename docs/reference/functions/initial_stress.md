# `imr.physics.initial_stress`

`imr.physics.initial_stress` computes the stress at the bubble's maximum radius for a collapse run. It shoots on the initial wall velocity of a Keller–Miksis growth from `Req`, so that the bubble reaches `R = 1` at its maximum. It returns the stress at that point. `imr.configuration.prepare_case` calls this function for stress 5 and 6 with `collapse = 1`, when the input has no `stress0`. The model has no heat or mass transfer. The bubble holds gas at the equilibrium pressure and vapor at `Pv_star`.

## Inputs

| Name      | Description                            | Units   |
|-----------|----------------------------------------|---------|
| `Req`     | Equilibrium radius                     | –       |
| `Re`      | Reynolds number                        | –       |
| `Ca`      | Cauchy number                          | –       |
| `De`      | Deborah number                         | –       |
| `We`      | Weber number                           | –       |
| `CL`      | Sound speed `C8/Uc` (`Cstar`)          | –       |
| `Pv_star` | Vapor pressure `Pv/P8`; sets the bubble pressure `Pv_star + (1 - Pv_star + 1/(We*Req))*(Req/R)^3` | –       |
| `alphax`  | Strain-stiffening parameter of the quadratic model (optional, default 0) | –       |

## Outputs

| Name | Description                      | Units   |
|------|----------------------------------|---------|
| `S0` | Initial stress at bubble maximum | –       |

For the algorithm and assumptions, see the [topic guide](../../method/initial-stress.md). See [historical reports](../../verification/convergence.md) for the numerical measurements.
