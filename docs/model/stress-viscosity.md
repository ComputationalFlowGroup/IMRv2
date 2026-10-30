# Stress and viscosity

`imr.physics.material_stress` supplies the integrated material stress to the radial equation. `imr.physics.viscosity` adds a shear-dependent contribution when `nu_model` is active.

| `stress` | Constitutive family | History |
|---|---|---|
| 0 | No stress | None |
| 1 | Neo-Hookean Kelvin–Voigt | Closed form |
| 2 | Quadratic neo-Hookean Kelvin–Voigt | Closed form |
| 3 | Linear Maxwell, Jeffreys or Zener with neo-Hookean spring | One stress integral, plus a local field for medium heating |
| 4 | Same family with quadratic spring | One stress integral, plus a local field for medium heating |
| 5 | UCM or Oldroyd-B | Two integrals or a spectral local field |
| 6 | PTT (linear Phan-Thien-Tanner) | Spectral local field only |
| 7 | Giesekus | Spectral local field only |
| 8 | 2nd-order Fung Kelvin–Voigt | Closed form |
| 9 | Linear Maxwell, Jeffreys or Zener with 2nd-order Fung spring | One stress integral, plus a local field for medium heating |
| 10 | Mooney-Rivlin Kelvin–Voigt | Closed form |
| 11 | Linear Maxwell, Jeffreys or Zener with Mooney-Rivlin spring | One stress integral, plus a local field for medium heating |

Actual labels depend on `Ca`, `Re8` and `LAM`; see the canonical [display-label table](../reference/functions/f_display.md#constitutive-model-label). Stress 5–7 force `Ca = Inf`. Stress 3–7, 9 and 11 require `De >= 1e-4`; use the small-De limit named by the checker below that threshold. `eps3` is accepted and displayed, and enters only the Giesekus model (stress 7); it is otherwise unused.

## PTT and Giesekus (stress 6, 7)

Stress 6 and 7 generalize the UCM/Oldroyd-B closure (stress 5) with the nonlinear terms of Warnez and Johnsen (2015) (their Eq. 13; see [references](../reference/references.md)), available only as the spectral stress field of `imr.solve_spectral` (`imr.configuration.validate_case` forces `spectral = 1`; neither model has a closed-form, two-ODE reduction, so `imr.solve_nodal` rejects them). Along a material point with `a = R^2*Rdot/r^3`, the field obeys

```text
D(tau_rr)/Dt = -tau_rr/De*exp(eps2*Re8*De*(tau_rr+2*tau_tt))
    - 4*a*tau_rr - eps3*Re8*tau_rr^2 - 4*etap*a/De
D(tau_tt)/Dt = -tau_tt/De*exp(eps2*Re8*De*(tau_rr+2*tau_tt))
    + 2*a*tau_tt - eps3*Re8*tau_tt^2 + 2*etap*a/De
```

with `etap = (1-LAM)/Re8`. Stress 5 is the `eps2 = 0`, `eps3 = 0` case, recovered exactly. Stress 6 (PTT) fixes their linear closure `eps2 = 1` and takes `eps3 = 0`. Stress 7 (Giesekus) keeps `eps2 = 0` and takes the input `eps3` in `(0, 0.5]` as their mobility factor; `eps3 = 0` makes stress 7 identical to stress 5 (a warning says so). `imr.solvers.build_spectral` (`f_spectral_stress`) implements the field, and `imr.numerics.jacobian_spectral` differentiates both new terms analytically for `jacobian = 2`.

## 2nd-order Fung (stress 8, 9)

Stress 2's quadratic spring is the first-order term of a series expansion of a Fung-type exponential strain energy in the invariant `x = I1 - 3`: its local elastic stress difference is `se = (lam^-4 - lam^2)/Ca*(1 + alphax*x)`, the derivative of `1 + alphax*x` truncated after its first term. Stress 8 (Kelvin-Voigt) and 9 (the same spring in the linear Maxwell/Jeffreys/Zener family of stress 3, 4) add the next term, `se = (lam^-4 - lam^2)/Ca*(1 + alphax*x + (alphax^2/2)*x^2)`, reusing `alphax`: no new input is needed. `alphax = 0` reduces stress 8 to stress 1 and stress 9 to stress 3 exactly, the same reduction stress 2 and 4 already have at `alphax = 0`.

The wall (radial-equation) stress integral `S` was cross-checked against `funIMRsolver.m`'s `fung2` branch ([`InertialMicrocavitationRheometry/IMR_v1`](https://github.com/InertialMicrocavitationRheometry/IMR_v1)), after correcting that source's reversed stretch convention (`R/Req` there, `Req/R` here). The local field `se` above was derived independently, as the next term of the same series stress 2 already truncates, then checked two ways before being trusted: its integral `S = 2*int_1^Inf se/y dy` reproduces the cross-checked wall formula exactly (not just approximately) over a range of `R`, `Req`, `Ca` and `alphax`, and its `Rst`-derivative matches the wall formula's own `Sdot` by construction (`imr.physics.material_stress`). Unlike stress 8's Tier-1 predecessor (a linear spring, since removed), this local field is exact, not a coarser approximation, so the medium heating of stress 8/9 is on the same footing as stress 1-4's.

## Mooney-Rivlin (stress 10, 11)

The two-invariant Mooney-Rivlin strain energy adds one term beyond neo-Hookean's single invariant. Its local elastic stress difference is `se = (1/Ca)*[(lam^-4 - lam^2) - C2*(lam^4 - lam^-2)]`, with `C2` the new, dimensionless second-modulus ratio (`C2 = 0` reduces exactly to neo-Hookean). Stress 10 (Kelvin-Voigt) and 11 (the same spring in the linear Maxwell/Jeffreys/Zener family) use it.

Neither source codebase implements Mooney-Rivlin, so both the local field and its wall (radial-equation) integral are derived here from scratch, not ported. The new `C2` term's integral collapses, under the substitution `v = lam` already used for the existing terms, to the elementary `2*int_1^Inf (lam^4 - lam^-2)/y dy = 2/Rst - Rst^2 - 1`, giving

```text
S = -(5 - 4*Rst - Rst^4)/(2*Ca) - (C2/Ca)*(2/Rst - Rst^2 - 1) - 4/Re8*Rdot/R - 6*intfnu*iDRe
```

Checked the same way as the Fung term: the field's own integral reproduces this closed form exactly, and its `Rst`-derivative matches the closed form's own `Sdot`, over a range of `R`, `Req`, `Ca` and `C2`.

## Non-Newtonian viscosity

`mu = mu8 + (mu0 - mu8)*f(x)`, where `x = lambda*|gammadot|` and `lambda` is the dimensionless `v_lambda` input.

| `nu_model` | Model | `f(x)` | Conditions |
|---|---|---|---|
| 0 | Newtonian | No additional contribution | Default |
| 1 | Carreau | `(1 + x^2)^((n-1)/2)` | `v_nc = n` |
| 2 | Carreau–Yasuda | `(1 + x^a)^((n-1)/a)` | `v_a = a > 0`, `v_nc = n` |
| 3 | Powell–Eyring | `asinh(x)/x` | Continuous value 1 at zero |
| 4 | Modified Powell–Eyring | `log(1 + x)/x^n` | `v_nc = n <= 1`, `v_lambda > 0` |
| 5 | Cross | `1/(1 + x^n)` | `v_nc = n > 0` |
| 6 | Simplified Cross | `1/(1 + x)` | No exponent |
| 7 | Modified Cross | `(1 + x^n)^(-a)` | `v_a > 0`, `v_nc > 0` |

The modified Powell–Eyring model has zero added viscosity at zero shear for `n < 1`, and factor 1 for `n = 1`; do not assume every exponent yields zero-shear viscosity `mu0`. The other laws approach factor 1 as `v_lambda` tends to zero.

The contribution vanishes at `DRe = 0`. For the Maxwell/Zener/UCM family (stress 3–7, 9, 11) it acts on the solvent share and also vanishes at `LAM = 0`. Its common factor is `iDRe = 1/DRe` for the Kelvin-Voigt family (stress 1, 2, 8, 10), `LAM/DRe` for stress 3–7, 9, 11, and zero when inactive. See [integral theory](../theory/viscosity-integrals.md), [heating theory](../theory/stress-heating.md) and [compatibility](../reference/compatibility.md).
