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

Actual labels depend on `Ca`, `Re8` and `LAM`; see the canonical [display-label table](../reference/functions/f_display.md#constitutive-model-label). Stress 5 forces `Ca = Inf`. Stress 3–5 require `De >= 1e-4`; use the small-De limit named by the checker below that threshold. `eps3` is accepted and displayed but does not enter the implemented stress equations.

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
