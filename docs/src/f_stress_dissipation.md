# `f_stress_dissipation.m` Documentation

## Overview

`f_stress_dissipation` computes the medium heating by the stress, `taugradu = Br*Phi`, when `medtherm = 1`. The medium temperature equation adds it: `Tmdot = advection + diffusion + taugradu`. Both solvers call it.

The heating is the full stress power, `Phi = tau:D`. This is the entropy-elastic (rubber-like) convention: elastic work is released as heat under loading and absorbed under unloading.

With the local strain rate `a = R^2*Rdot/r^3 = (Rdot/R)*yT^-3` and `s = tau_rr - tau_thetatheta`, the stress power is `tau:D = -2*a*s`. Only this normal-stress difference enters; the trace of the stress does not.

## Heating per model

| `stress` | `taugradu / Br` |
|---|---|
| 0 | `0` |
| 1 | `12*mu_KV*a^2 - 2*a*se`, neo-Hookean `se` |
| 2 | the same with the quadratic `se` |
| 3, 4 | `12*mu_s*a^2 - 2*a*sp`, with the Zener polymer field `sp` of [`f_zener_field`](f_zener_field.md) |
| 5 | `12*mu_s*a^2 - 2*a*(tau_rr - tau_tt)`, from the spectral stress field |

**Viscosities**, with `fnu` the non-Newtonian factor at the medium points:

- `mu_KV = 1/Re8 + fnu*iDRe`: stress 1 and 2, where the whole viscosity acts instantly.
- `mu_s = LAM/Re8 + fnu*iDRe`: stress 3–5, where only the solvent part acts instantly.

`iDRe` is `1/DRe` for stress 1–2 and `LAM/DRe` for stress 3–5, the same factor as in the radial equation (set in `f_call_params`).

**Elastic stress difference.** With the hoop stretch `lam = r/r0`, where `r0^3 = r^3 - R^3 + Req^3`:

- neo-Hookean: `se = (lam^-4 - lam^2)/Ca`;
- quadratic: the same times `1 + alphax*(lam^-4 + 2*lam^2 - 3)`.

In the code, `lam^2 = yT2./x2` and `lam^-4 = x2.^2.*iyT4`, with `x2 = (yT3 - 1 + (Req/R)^3).^(2/3)`.

## Total heating

Integrating over the medium gives `4*pi*R^3*int_1^Inf y^2*Phi dy = -4*pi*R^2*Rdot*S`, with `S` the stress integral of the radial equation. This holds for every model. For stress 3 and 4 it equals `16*pi*LAM/Re8*R*Rdot^2 + 24*pi*R^2*Rdot*intfnu*iDRe - 4*pi*Rdot*Z1/R`.

## Inputs

| Name | Description |
|---|---|
| `stress` | Stress model (0–5) |
| `Req`, `R`, `Rdot` | Equilibrium radius, radius and wall velocity |
| `Ca`, `Br`, `Re8`, `LAM` | Cauchy, Brinkman and Reynolds numbers, and the solvent fraction |
| `alphax` | Strain-stiffening parameter (stress 2) |
| `iDRe` | Factor of the non-Newtonian stress term (`1/DRe` or `LAM/DRe`, 0 when inactive) |
| `yT2`, `yT3`, `iyT3`, `iyT4`, `iyT6` | Powers of the medium coordinate `y = r/R` |
| `fnu` | Non-Newtonian viscosity factor at the medium points (the `fq` output of `f_viscosity`); scalar 0 for Newtonian runs |
| `sp` | Local polymer stress difference at the medium points (stress 3–5); ignored otherwise |

## Outputs

| Name | Description |
|---|---|
| `taugradu` | Medium heating `Br*tau:D` at the medium points |

## Notes

- **Where `sp` comes from:**
    - `f_imr_spectral` evaluates the Zener field or the stress 5 field at the medium points with `ZZT`.
    - `f_imr_nodal` stores the Zener field at the medium nodes.
- **Stress 5 two-state form:** it carries only integrals of the stress, not a local field. So stress 5 with medium heat transfer uses the stress field in `f_imr_spectral`, and is rejected in `f_imr_nodal` (see [`f_checkers`](f_checkers.md)).
- **Far-field point:** at `yT = Inf` the elastic term is undefined. The solvers hold the temperature there fixed (`Tmdot(end) = 0`).
- **Batched states:** `R`, `Rdot` and `sp` may hold one state per column, for the batched right-hand side of [`f_imr_nodal`](f_imr_nodal.md#batched-right-hand-side).

## Verification (GNU Octave)

- **Energy budget.** The total heating from quadrature on the solver grids matches `-4*pi*R^2*Rdot*S` to:
    - `1e-8`–`5e-5` (relative) with Chebyshev and spectral grids;
    - `5e-4` with finite differences at `Mt = 80` (second-order quadrature error).

    This held for stress 1–5, Newtonian and Carreau.
- **Kelvin–Voigt limit.** At the same state, the stress 3 (4) heating tends to the stress 1 (2) heating as `De -> 0`: the relative difference is `4e-2`, `4e-3` and `4e-4` for `De = 1e-2`, `1e-3` and `1e-4`.
