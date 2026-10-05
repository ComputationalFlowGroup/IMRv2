# `f_init_stress.m` Documentation

## Overview
`f_init_stress` computes the stress at the bubble's maximum radius for a collapse run. It shoots on the initial wall velocity of a Keller–Miksis growth from `Req` so that the bubble reaches `R = 1` at its maximum, and returns the stress there. `f_call_params` calls it for stress 3 and 4 with `collapse = 1` when no `stress0` is given. The model has no heat or mass transfer; the bubble holds gas at the equilibrium pressure and vapor at `Pv_star`.

## Major Sections
1. **Bracket**: the maximum radius grows with the initial wall velocity and equals `Req` at rest. Starting from `[0, 1]`, the upper end doubles until the maximum radius exceeds 1. The growth equation is singular at `Rdot = CL`, so the bracket stays below the sound speed; if even a nearly sonic start stays below `R = 1`, the function stops with an error that asks for `stress0`.
2. **Shooting**: `fzero` on the bracket finds the initial wall velocity for which `f_maxwell_growth_iter` reaches a maximum radius of 1.
3. **Time Integration**: `ode23tb` on `f_RP_bub_collapse_comp`, with an event at `Rdot = 0` to locate the maximum radius.
4. **Stress Extraction**: the stress state at the maximum radius.

### Internal Functions
- `f_maxwell_growth_iter`: returns the maximum radius for a trial initial wall velocity.
- `f_RP_bub_collapse_comp`: Nondimensional Keller–Miksis ODE for the growth with the stress evolution. For the Maxwell branch the stress relaxes toward the elastic stress `Se` of the (quadratic) neo-Hookean model, the same `Ze/R^3` as stress 3 and 4 in [`f_stress`](f_stress.md): `Sdot = -4*Rdot/(Re*R*De) + dSe/dt - (S - Se)/De`.

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

## Notes
- Against the previous `fminsearch` shooting, over the default collapse parameters with `Req` 0.1–0.5 and `De` 0.01–1, the stress agrees to `3e-8` relative. At `Req = 0.05` no subsonic start reaches `R = 1` (the maximum radius stays at 0.93); the previous version returned the stress at that non-root without notice.
- Integrates with `ode23tb`.
