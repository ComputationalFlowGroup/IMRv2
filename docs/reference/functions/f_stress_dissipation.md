# `imr.physics.stress_power`

`imr.physics.stress_power` computes the medium heating by the stress, `taugradu = Br*Phi`, when `medtherm = 1`. The medium temperature equation adds it: `Tmdot = advection + diffusion + taugradu`. Both solvers call it.

The heating is the full stress power, `Phi = tau:D`. This is the entropy-elastic (rubber-like) convention: elastic work is released as heat under loading and absorbed under unloading.

With the local strain rate `a = R^2*Rdot/r^3 = (Rdot/R)*yT^-3` and `s = tau_rr - tau_thetatheta`, the stress power is `tau:D = -2*a*s`. Only this normal-stress difference enters; the trace of the stress does not.

## Inputs

| Name | Description |
|---|---|
| `stress` | Stress model (0–11) |
| `Req`, `R`, `Rdot` | Equilibrium radius, radius and wall velocity |
| `Ca`, `Br`, `Re8`, `LAM` | Cauchy, Brinkman and Reynolds numbers, and the solvent fraction |
| `alphax` | Strain-stiffening parameter (stress 2, 6); the next series term reuses it for the 2nd-order Fung spring (3, 7) |
| `C2` | Mooney-Rivlin second-modulus ratio (stress 4, 8); unused otherwise |
| `iDRe` | Factor of the non-Newtonian stress term (`1/DRe` or `LAM/DRe`, 0 when inactive) |
| `yT2`, `yT3`, `iyT3`, `iyT4`, `iyT6` | Powers of the medium coordinate `y = r/R` |
| `fnu` | Non-Newtonian viscosity factor at the medium points (the `fq` output of `imr.physics.viscosity`); scalar 0 for Newtonian runs |
| `sp` | Local polymer stress difference at the medium points (stress 5–11); ignored otherwise |

## Outputs

| Name | Description |
|---|---|
| `taugradu` | Medium heating `Br*tau:D` at the medium points |

## Notes

- **Where `sp` comes from:**
    - `imr.solve_spectral` evaluates the Zener field or the stress 9–11 field at the medium points with `ZZT`.
    - `imr.solve_nodal` stores the Zener field at the medium nodes.
- **Stress 9 two-state form:** it carries only integrals of the stress, not a local field. So stress 9 with medium heat transfer uses the stress field in `imr.solve_spectral`, and is rejected in `imr.solve_nodal` (see [`imr.configuration.validate_case`](f_checkers.md)). Stress 10 and 11 always run as that stress field, so this does not arise for them.
- **Stress 3, 4 (Kelvin-Voigt):** the elastic stress difference `se` is the 2nd-order Fung or Mooney-Rivlin field in place of the neo-Hookean or quadratic one; the instantaneous viscosity and non-Newtonian terms are unchanged.
- **Stress 7, 8 (Zener):** carry the corresponding spring's local field exactly as stress 5, 6 do, via `imr.physics.zener_field`.
- **Far-field point:** at `yT = Inf` the elastic term is undefined. The solvers hold the temperature there fixed (`Tmdot(end) = 0`).
- **Batched states:** `R`, `Rdot` and `sp` may hold one state per column, for the batched right-hand side of [`imr.solve_nodal`](../../method/time-integration.md#batched-right-hand-side).

For the algorithm and assumptions, see the [topic guide](../../theory/stress-heating.md). Numerical measurements are retained as [historical reports](../../validation/convergence.md).
