# `imr.physics.stress_power`

`imr.physics.stress_power` computes the medium heating by the stress, `taugradu = Br*Phi`, when `medtherm = 1`. The medium temperature equation adds it: `Tmdot = advection + diffusion + taugradu`. Both solvers call it.

The heating is the full stress power, `Phi = tau:D`. This is the entropy-elastic (rubber-like) convention: elastic work is released as heat under loading and absorbed under unloading.

With the local strain rate `a = R^2*Rdot/r^3 = (Rdot/R)*yT^-3` and `s = tau_rr - tau_thetatheta`, the stress power is `tau:D = -2*a*s`. Only this normal-stress difference enters; the trace of the stress does not.

## Inputs

| Name | Description |
|---|---|
| `stress` | Stress model (0–5) |
| `Req`, `R`, `Rdot` | Equilibrium radius, radius and wall velocity |
| `Ca`, `Br`, `Re8`, `LAM` | Cauchy, Brinkman and Reynolds numbers, and the solvent fraction |
| `alphax` | Strain-stiffening parameter (stress 2) |
| `iDRe` | Factor of the non-Newtonian stress term (`1/DRe` or `LAM/DRe`, 0 when inactive) |
| `yT2`, `yT3`, `iyT3`, `iyT4`, `iyT6` | Powers of the medium coordinate `y = r/R` |
| `fnu` | Non-Newtonian viscosity factor at the medium points (the `fq` output of `imr.physics.viscosity`); scalar 0 for Newtonian runs |
| `sp` | Local polymer stress difference at the medium points (stress 3–5); ignored otherwise |

## Outputs

| Name | Description |
|---|---|
| `taugradu` | Medium heating `Br*tau:D` at the medium points |

## Notes

- **Where `sp` comes from:**
    - `imr.solve_spectral` evaluates the Zener field or the stress 5 field at the medium points with `ZZT`.
    - `imr.solve_nodal` stores the Zener field at the medium nodes.
- **Stress 5 two-state form:** it carries only integrals of the stress, not a local field. So stress 5 with medium heat transfer uses the stress field in `imr.solve_spectral`, and is rejected in `imr.solve_nodal` (see [`imr.configuration.validate_case`](f_checkers.md)).
- **Far-field point:** at `yT = Inf` the elastic term is undefined. The solvers hold the temperature there fixed (`Tmdot(end) = 0`).
- **Batched states:** `R`, `Rdot` and `sp` may hold one state per column, for the batched right-hand side of [`imr.solve_nodal`](../../method/time-integration.md#batched-right-hand-side).

For the algorithm and assumptions, see the [topic guide](../../theory/stress-heating.md). Numerical measurements are retained as [historical reports](../../validation/convergence.md).
