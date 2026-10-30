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

Actual labels depend on `Ca`, `Re8` and `LAM`; see the canonical [display-label table](../reference/functions/f_display.md#constitutive-model-label). Stress 5–7 force `Ca = Inf`. Stress 3–7 require `De >= 1e-4`; use the small-De limit named by the checker below that threshold. `eps3` is accepted and displayed, and enters only the Giesekus model (stress 7); it is otherwise unused.

## PTT and Giesekus (stress 6, 7)

Stress 6 and 7 generalize the UCM/Oldroyd-B closure (stress 5) with the nonlinear terms of Warnez and Johnsen (2015) (their Eq. 13; see [references](../reference/references.md)), available only as the spectral stress field of `imr.solve_spectral` (`imr.configuration.validate_case` forces `spectral = 1`; neither model has a closed-form, two-ODE reduction, so `imr.solve_nodal` rejects them). Along a material point with `a = R^2*Rdot/r^3`, the field obeys

```text
D(tau_rr)/Dt = -tau_rr/De*exp(eps2*Re8*De*(tau_rr+2*tau_tt))
    - 4*a*tau_rr - eps3*Re8*tau_rr^2 - 4*etap*a/De
D(tau_tt)/Dt = -tau_tt/De*exp(eps2*Re8*De*(tau_rr+2*tau_tt))
    + 2*a*tau_tt - eps3*Re8*tau_tt^2 + 2*etap*a/De
```

with `etap = (1-LAM)/Re8`. Stress 5 is the `eps2 = 0`, `eps3 = 0` case, recovered exactly. Stress 6 (PTT) fixes their linear closure `eps2 = 1` and takes `eps3 = 0`. Stress 7 (Giesekus) keeps `eps2 = 0` and takes the input `eps3` in `(0, 0.5]` as their mobility factor; `eps3 = 0` makes stress 7 identical to stress 5 (a warning says so). `imr.solvers.build_spectral` (`f_spectral_stress`) implements the field, and `imr.numerics.jacobian_spectral` differentiates both new terms analytically for `jacobian = 2`.

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

The contribution vanishes at `DRe = 0`. For stress 3–5 it acts on the solvent share and also vanishes at `LAM = 0`. Its common factor is `iDRe = 1/DRe` for stress 1–2, `LAM/DRe` for stress 3–5, and zero when inactive. See [integral theory](../theory/viscosity-integrals.md), [heating theory](../theory/stress-heating.md) and [compatibility](../reference/compatibility.md).
