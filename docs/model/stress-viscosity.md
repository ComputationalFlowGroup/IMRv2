# Stress and viscosity

`imr.physics.material_stress` supplies the integrated material stress to the radial equation. `imr.physics.viscosity` adds a shear-dependent contribution when `nu_model` is active.

| `stress` | Constitutive family | History |
|---|---|---|
| 0 | No stress | None |
| 1 | Neo-Hookean Kelvin–Voigt | Closed form |
| 2 | Quadratic neo-Hookean Kelvin–Voigt | Closed form |
| 3 | 2nd-order Fung Kelvin–Voigt | Closed form |
| 4 | Mooney-Rivlin Kelvin–Voigt | Closed form |
| 5 | Linear Maxwell, Jeffreys or Zener with neo-Hookean spring | One stress integral, plus a local field for medium heating |
| 6 | Same family with quadratic spring | One stress integral, plus a local field for medium heating |
| 7 | Same family with 2nd-order Fung spring | One stress integral, plus a local field for medium heating |
| 8 | Same family with Mooney-Rivlin spring | One stress integral, plus a local field for medium heating |
| 9 | UCM or Oldroyd-B | Two integrals or a spectral local field |
| 10 | PTT (linear Phan-Thien-Tanner) | Spectral local field only |
| 11 | Giesekus | Spectral local field only |

The numbering groups by family and complexity: 1–4 the Kelvin-Voigt family (no stress history), 5–8 the matching Maxwell/Jeffreys/Zener family (one stress integral), 9 UCM/Oldroyd-B (two integrals), 10–11 PTT/Giesekus (spectral field only, no closed form). Within each spring-bearing family the spring itself goes neo-Hookean, quadratic, 2nd-order Fung, Mooney-Rivlin.

Actual labels depend on `Ca`, `Re8` and `LAM`; see the canonical [display-label table](../reference/functions/describe_case.md#constitutive-model-label). Stress 9–11 force `Ca = Inf`. Stress 5–11 require `De >= 1e-4`; use the small-De limit named by the checker below that threshold. `eps3` is accepted and displayed, and enters only the Giesekus model (stress 11); it is otherwise unused.

## 2nd-order Fung (stress 3, 7)

Stress 2's quadratic spring is the first-order term of a series expansion of a Fung-type exponential strain energy in the invariant $x=I_1-3$: its local elastic stress difference is

$$
s_e=\frac{\lambda^{-4}-\lambda^2}{\mathrm{Ca}}(1+\alpha_x x),
$$

the derivative of $1+\alpha_x x$ truncated after its first term. Stress 3 (Kelvin-Voigt) and 7 (the same spring in the linear Maxwell/Jeffreys/Zener family of stress 5, 6) add the next term,

$$
s_e=\frac{\lambda^{-4}-\lambda^2}{\mathrm{Ca}}\left(1+\alpha_x x+\frac{\alpha_x^2}{2}x^2\right),
$$

reusing $\alpha_x$: no new input is needed. $\alpha_x=0$ reduces stress 3 to stress 1 and stress 7 to stress 5 exactly, the same reduction stress 2 and 6 already have at $\alpha_x=0$.

The wall (radial-equation) stress integral $S$ was cross-checked against `funIMRsolver.m`'s `fung2` branch ([`InertialMicrocavitationRheometry/IMR_v1`](https://github.com/InertialMicrocavitationRheometry/IMR_v1)), after correcting that source's reversed stretch convention ($R/R_{eq}$ there, $R_{eq}/R$ here). The local field $s_e$ above was derived independently, as the next term of the same series stress 2 already truncates, then checked two ways before being trusted: its integral $S=2\int_1^\infty s_e/y\,dy$ reproduces the cross-checked wall formula exactly (not just approximately) over a range of $R$, $R_{eq}$, $\mathrm{Ca}$ and $\alpha_x$, and its $R_{eq}/R$-derivative matches the wall formula's own `Sdot` by construction (`imr.physics.material_stress`). This local field is exact, not an approximation, so the medium heating of stress 3/7 is on the same footing as stress 1, 2, 5, 6's.

## Mooney-Rivlin (stress 4, 8)

The two-invariant Mooney-Rivlin strain energy adds one term beyond neo-Hookean's single invariant. Its local elastic stress difference is

$$
s_e=\frac{1}{\mathrm{Ca}}\Big[(\lambda^{-4}-\lambda^2)-C_2(\lambda^4-\lambda^{-2})\Big],
$$

with $C_2$ the new, dimensionless second-modulus ratio ($C_2=0$ reduces exactly to neo-Hookean). Stress 4 (Kelvin-Voigt) and 8 (the same spring in the linear Maxwell/Jeffreys/Zener family) use it.

Neither source codebase implements Mooney-Rivlin, so both the local field and its wall (radial-equation) integral are derived here from scratch, not ported. The new $C_2$ term's integral collapses, under the substitution $v=\lambda$ already used for the existing terms, to the elementary

$$
2\int_1^\infty \frac{\lambda^4-\lambda^{-2}}{y}\,dy = \frac{2}{R_{st}}-R_{st}^2-1,\qquad R_{st}=\frac{R_{eq}}{R},
$$

giving

$$
S = -\frac{5-4R_{st}-R_{st}^4}{2\,\mathrm{Ca}} - \frac{C_2}{\mathrm{Ca}}\left(\frac{2}{R_{st}}-R_{st}^2-1\right) - \frac{4}{\mathrm{Re}_8}\frac{\dot R}{R} - 6\,\mathtt{intfnu}\,\mathtt{iDRe}.
$$

Checked the same way as the Fung term: the field's own integral reproduces this closed form exactly, and its $R_{st}$-derivative matches the closed form's own `Sdot`, over a range of $R$, $R_{eq}$, $\mathrm{Ca}$ and $C_2$.

## PTT and Giesekus (stress 10, 11)

Stress 10 and 11 generalize the UCM/Oldroyd-B closure (stress 9) with the nonlinear terms of Warnez and Johnsen (2015) (their Eq. 13; see [references](../reference/references.md)), available only as the spectral stress field of `imr.solve_spectral` (`imr.configuration.validate_case` forces `spectral = 1`; neither model has a closed-form, two-ODE reduction, so `imr.solve_nodal` rejects them). Along a material point with $a=R^2\dot R/r^3$, the field obeys

$$
\frac{D\tau_{rr}}{Dt} = -\frac{\tau_{rr}}{\mathrm{De}}\exp\!\big(\epsilon_2\mathrm{Re}_8\mathrm{De}(\tau_{rr}+2\tau_{\theta\theta})\big) - 4a\tau_{rr} - \epsilon_3\mathrm{Re}_8\tau_{rr}^2 - \frac{4\eta_p a}{\mathrm{De}},
$$

$$
\frac{D\tau_{\theta\theta}}{Dt} = -\frac{\tau_{\theta\theta}}{\mathrm{De}}\exp\!\big(\epsilon_2\mathrm{Re}_8\mathrm{De}(\tau_{rr}+2\tau_{\theta\theta})\big) + 2a\tau_{\theta\theta} - \epsilon_3\mathrm{Re}_8\tau_{\theta\theta}^2 + \frac{2\eta_p a}{\mathrm{De}},
$$

with $\eta_p=(1-\mathrm{LAM})/\mathrm{Re}_8$. Stress 9 is the $\epsilon_2=0$, $\epsilon_3=0$ case, recovered exactly. Stress 10 (PTT) fixes their linear closure $\epsilon_2=1$ and takes $\epsilon_3=0$. Stress 11 (Giesekus) keeps $\epsilon_2=0$ and takes the input `eps3` in $(0,0.5]$ as their mobility factor; `eps3 = 0` makes stress 11 identical to stress 9 (a warning says so). `imr.solvers.build_spectral` (`f_spectral_stress`) implements the field, and `imr.numerics.jacobian_spectral` differentiates both new terms analytically for `jacobian = 2`.

## Non-Newtonian viscosity

$$
\mu = \mu_8 + (\mu_0-\mu_8)f(x),
$$

where $x=\lambda|\dot\gamma|$ and $\lambda$ is the dimensionless `v_lambda` input.

| `nu_model` | Model | $f(x)$ | Conditions |
|---|---|---|---|
| 0 | Newtonian | No additional contribution | Default |
| 1 | Carreau | $(1+x^2)^{(n-1)/2}$ | `v_nc = n` |
| 2 | Carreau–Yasuda | $(1+x^a)^{(n-1)/a}$ | `v_a = a > 0`, `v_nc = n` |
| 3 | Powell–Eyring | $\sinh^{-1}(x)/x$ | Continuous value 1 at zero |
| 4 | Modified Powell–Eyring | $\ln(1+x)/x^n$ | `v_nc = n <= 1`, `v_lambda > 0` |
| 5 | Cross | $1/(1+x^n)$ | `v_nc = n > 0` |
| 6 | Simplified Cross | $1/(1+x)$ | No exponent |
| 7 | Modified Cross | $(1+x^n)^{-a}$ | `v_a > 0`, `v_nc > 0` |

The modified Powell–Eyring model has zero added viscosity at zero shear for `n < 1`, and factor 1 for `n = 1`; do not assume every exponent yields zero-shear viscosity `mu0`. The other laws approach factor 1 as `v_lambda` tends to zero.

The contribution vanishes at `DRe = 0`. For the Maxwell/Zener/UCM family (stress 5–11) it acts on the solvent share and also vanishes at `LAM = 0`. Its common factor is `iDRe = 1/DRe` for the Kelvin-Voigt family (stress 1–4), `LAM/DRe` for stress 5–11, and zero when inactive. See [integral theory](../theory/viscosity-integrals.md), [heating theory](../theory/stress-heating.md) and [compatibility](../reference/compatibility.md).
