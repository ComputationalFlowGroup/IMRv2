# `f_init_stress.m` Documentation

## Overview
`f_init_stress` computes the stress at the bubble's maximum radius for a collapse run. It shoots on the initial wall velocity of a Rayleigh–Plesset growth from `Req` so that the bubble reaches `R = 1` at its maximum, and returns the stress there. `f_call_params` calls it for stress 3 and 4 with `collapse = 1` when no `stress0` is given. The model has no heat or mass transfer.

## Major Sections
1. **Shooting Setup**: Configure `fminsearch` options (`sopts`).
2. **Shooting**: `fminsearch` finds the initial wall velocity for which `f_maxwell_growth_iter` reaches a maximum radius of 1.
3. **Time Integration**: `ode23tb` on `f_RP_bub_collapse_comp`, with an event at `Rdot = 0` to locate the maximum radius.
4. **Stress Extraction**: the stress state at the maximum radius.

### Internal Functions
- `f_maxwell_growth_iter`: returns the maximum radius for a trial initial wall velocity.
- `f_RP_bub_collapse_comp`: Nondimensional Rayleigh–Plesset ODE for collapse with stress evolution.

## Inputs

| Name      | Description                            | Units   |
|-----------|----------------------------------------|---------|
| `Req`     | Equilibrium radius                     | –       |
| `Re`      | Reynolds number                        | –       |
| `Ca`      | Cauchy number                          | –       |
| `De`      | Deborah number                         | –       |
| `We`      | Weber number                           | –       |
| `CL`      | Sound speed `C8/Uc` (`Cstar`)          | –       |
| `Pv_star` | Vapor pressure `Pv/P8` (currently unused) | –       |

## Outputs

| Name | Description                      | Units   |
|------|----------------------------------|---------|
| `S0` | Initial stress at bubble maximum | –       |

## Notes
- Uses `ode23tb`, so it needs MATLAB; GNU Octave has no `ode23tb`.
