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

## Neo-Hookean and quadratic (stress 1, 2, 5, 6)

With the hoop stretch $\lambda=r/r_0$, where $r_0^3=r^3-R^3+R_{eq}^3$, and $x=\lambda^{-4}+2\lambda^2-3$ (the invariant term $I_1-3$), the neo-Hookean local elastic stress difference is

$$
s_e=\frac{\lambda^{-4}-\lambda^2}{\mathrm{Ca}}
$$

for stress 1 (Kelvin-Voigt) and 5 (the same spring in the linear Maxwell/Jeffreys/Zener family). The quadratic spring (stress 2, 6) adds the strain-stiffening factor $1+\alpha_x x$, the first-order term of the series $\exp(\alpha_x x)$:

$$
s_e=\frac{\lambda^{-4}-\lambda^2}{\mathrm{Ca}}(1+\alpha_x x).
$$

## 2nd-order Fung (stress 3, 7)

Stress 3 (Kelvin-Voigt) and 7 (the same spring in the linear Maxwell/Jeffreys/Zener family of stress 5, 6) add the next term of the series stress 2 and 6 already truncate after their first order,

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

In the code, `lam^2 = yT2./x2` and `lam^-4 = x2.^2.*iyT4`, with `x2 = (yT3 - 1 + (Req/R)^3).^(2/3)`; `lam^4 = (yT2.*ix2).^2` and `lam^-2 = x2./yT2` for the Mooney-Rivlin term.

## PTT and Giesekus (stress 10, 11)

Stress 10 and 11 generalize the UCM/Oldroyd-B closure (stress 9) with the nonlinear terms of Warnez and Johnsen (2015) (their Eq. 13; see [references](../reference/references.md)), available only as the spectral stress field of `imr.solve_spectral` (`imr.configuration.validate_case` forces `spectral = 1`; neither model has a closed-form, two-ODE reduction, so `imr.solve_nodal` rejects them). Along a material point with $a=R^2\dot R/r^3$, the field obeys

$$
\frac{D\tau_{rr}}{Dt} = -\frac{\tau_{rr}}{\mathrm{De}}\exp\!\big(\epsilon_2\mathrm{Re}_8\mathrm{De}(\tau_{rr}+2\tau_{\theta\theta})\big) - 4a\tau_{rr} - \epsilon_3\mathrm{Re}_8\tau_{rr}^2 - \frac{4\eta_p a}{\mathrm{De}},
$$

$$
\frac{D\tau_{\theta\theta}}{Dt} = -\frac{\tau_{\theta\theta}}{\mathrm{De}}\exp\!\big(\epsilon_2\mathrm{Re}_8\mathrm{De}(\tau_{rr}+2\tau_{\theta\theta})\big) + 2a\tau_{\theta\theta} - \epsilon_3\mathrm{Re}_8\tau_{\theta\theta}^2 + \frac{2\eta_p a}{\mathrm{De}},
$$

with $\eta_p=(1-\mathrm{LAM})/\mathrm{Re}_8$. Stress 9 is the $\epsilon_2=0$, $\epsilon_3=0$ case, recovered exactly. Stress 10 (PTT) fixes their linear closure $\epsilon_2=1$ and takes $\epsilon_3=0$. Stress 11 (Giesekus) keeps $\epsilon_2=0$ and takes the input `eps3` in $(0,0.5]$ as their mobility factor; `eps3 = 0` makes stress 11 identical to stress 9 (a warning says so). `imr.solvers.build_spectral` (`f_spectral_stress`) implements the field, and `imr.numerics.jacobian_spectral` differentiates both new terms analytically for `jacobian = 2`.

## Medium heating

`imr.physics.stress_power` computes the medium heating by the stress, $\mathtt{taugradu}=\mathrm{Br}\,\Phi$, when `medtherm = 1`. The medium temperature equation adds it: $\dot T_m = \text{advection}+\text{diffusion}+\mathtt{taugradu}$. Both solvers call it. The heating is the full stress power, $\Phi=\boldsymbol\tau:\mathbf D$ — the entropy-elastic (rubber-like) convention: elastic work is released as heat under loading and absorbed under unloading.

With the local strain rate $a=\dfrac{R^2\dot R}{r^3}=\dfrac{\dot R}{R}y^{-3}$ ($y=r/R$) and $s=\tau_{rr}-\tau_{\theta\theta}$, the stress power is $\boldsymbol\tau:\mathbf D=-2as$. Only this normal-stress difference enters; the trace of the stress does not.

| `stress` | $\mathtt{taugradu}/\mathrm{Br}$ |
|---|---|
| 0 | $0$ |
| 1 | $12\mu_{KV}a^2-2as_e$, neo-Hookean $s_e$ |
| 2 | the same with the quadratic $s_e$ |
| 3 | the same with the 2nd-order Fung $s_e$ |
| 4 | the same with the Mooney-Rivlin $s_e$ |
| 5, 6, 7, 8 | $12\mu_s a^2-2as_p$, with the Zener polymer field $s_p$ of [`imr.physics.zener_field`](../reference/functions/zener_field.md) |
| 9, 10, 11 | $12\mu_s a^2-2a(\tau_{rr}-\tau_{\theta\theta})$, from the spectral stress field (UCM/Oldroyd-B, PTT or Giesekus) |

Each $s_e$ above is the spring formula of the matching stress 1–4 model, stated in the sections above. **Viscosities**, with $f_\nu$ the non-Newtonian factor at the medium points:

- $\mu_{KV}=1/\mathrm{Re}_8+f_\nu\,\mathtt{iDRe}$: stress 1–4, where the whole viscosity acts instantly.
- $\mu_s=\mathrm{LAM}/\mathrm{Re}_8+f_\nu\,\mathtt{iDRe}$: stress 5–11, where only the solvent part acts instantly.

`iDRe` is $1/\mathtt{DRe}$ for stress 1–4 and $\mathrm{LAM}/\mathtt{DRe}$ for stress 5–11, the same factor as in the radial equation (set in `imr.configuration.prepare_case`).

### Zener local field

For stress 5–8, integrating the heating over the medium gives

$$
4\pi R^3\int_1^\infty y^2\Phi\,dy = -4\pi R^2\dot R\,S,
$$

with $S$ the stress integral of the radial equation. This holds for every model. For stress 5, 6, 7 and 8 it equals

$$
16\pi\frac{\mathrm{LAM}}{\mathrm{Re}_8}R\dot R^2+24\pi R^2\dot R\,\mathtt{intfnu}\,\mathtt{iDRe}-\frac{4\pi\dot R Z_1}{R}.
$$

On the grid $y=r/R$, the local polymer stress difference $s_p$ obeys

$$
\left.\frac{\partial s_p}{\partial t}\right|_y = -\dot y\,\frac{\partial s_p}{\partial y} - \frac{s_p-s_{KV}}{\mathrm{De}} - 3\frac{\dot R}{R}\left(1-y^{-3}\right)s_p,
$$

$$
\dot y = \frac{\dot R}{R}\left(y^{-2}-y\right),
$$

$$
s_{KV} = s_e - 6\eta_p a,\qquad a=\frac{\dot R}{R}y^{-3},\qquad \eta_p=\frac{1-\mathrm{LAM}}{\mathrm{Re}_8}.
$$

- $s_e$ is the neo-Hookean elastic stress difference for stress 5, the quadratic one for stress 6, the 2nd-order Fung one for stress 7, or the Mooney-Rivlin one for stress 8.
- The caller computes the advection $-\dot y\,\partial s_p/\partial y$ on its own grid.
- **Boundary conditions:** none at the wall, where $\dot y=0$. $s_p$ is 0 at infinity.
- **Initial condition:** $s_p=0$ from rest.

**Why this law.** It is the local law whose integral reproduces the solver's stress state exactly:

$$
2\int_1^\infty \frac{s_p}{y}\,dy = \frac{Z_1}{R^3},
$$

where $Z_1$ obeys

$$
\mathrm{De}\,\dot Z_1 + Z_1 = Z_e - 4\eta_p R^2\dot R
$$

with $Z_e = R^3 S_e$. Written for the weighted field $w=y^{-3}s_p$, the law reads $\dfrac{Dw}{Dt} = -\dfrac{w-y^{-3}s_{KV}}{\mathrm{De}}$ along a material point. From rest its solution is

$$
y^{-3}s_p(t) = \int_0^t K(t-t')\,y'^{-3}s_{KV}(t')\,dt',\qquad K(t)=\frac{e^{-t/\mathrm{De}}}{\mathrm{De}},
$$

where primes denote the same material point at time $t'$. Incompressibility gives $r^3-r'^3=R^3-R'^3$, so $dr/r = r'^2\,dr'/r^3$. Hence

$$
2\int_R^\infty \frac{s_p}{r}\,dr = \frac{1}{R^3}\int_0^t K\,R'^3\left[2\int_{R'}^\infty \frac{s_{KV}}{r'}\,dr'\right]dt' = \frac{1}{R^3}\int_0^t K\,R'^3\,S_{KV}(t')\,dt',
$$

with $S_{KV}=S_e-4\eta_p\dot R/R$, the stress integral of $s_{KV}$. This is $Z_1/R^3$, the solution of the $Z_1$ equation from rest.

The $Z_1$ equation fixes only this weighted integral. The law above is the simplest local field consistent with it, not the only one. It sets the radial distribution of the heating, but neither $R(t)$ nor the total heating. So: $Z_1$ still drives $R(t)$, unchanged; the field only distributes the heating in $r$; the total heating is fixed by $Z_1$.

**In the solvers.**

- `imr.solve_spectral`: $s_p$ is expanded in the stress 9 basis $T_n(z)-1$ with `Nv` modes, stored after $Z_1$, and evaluated at the medium points with `ZZT`. After a run, the solver prints the largest mismatch of $2\,\mathtt{cdd}'e$ against $Z_1/R^3$ as a check.
- `imr.solve_nodal`: $s_p$ is stored at the medium nodes `yT` after $Z_1$: the first node is the wall, and the last, at infinity, is held at 0. The advection reuses the $\dot R$ part of the medium temperature advection, `xi1_sq./(Lt*R).*(Rdot./yT2.*om_yT3/2).*(dmatrix_Tm*sp)`, which equals $-\dot y\,\partial s_p/\partial y$.
- Collapse (nodal only): the history of the field before $t=0$ is not computed. As a stopgap, the initial field takes the shape of the elastic stress difference of the initial radius, scaled so that the integral identity holds for the initial $Z_1$ from `imr.physics.initial_stress`. Without elasticity ($\mathrm{Ca}=\infty$), it takes the viscous shape $y^{-3}$ instead.

### Spectral stress field

For stress 9, 10 and 11, `imr.solve_spectral` expands $\tau_{rr}$ and $\tau_{\theta\theta}$ in $T_n(z)-1$, $n=1,\ldots,N_v$, with

$$
z = 1 - \frac{2}{1 + (r/R - 1)/L_v}
$$

and collocation points $z_j=\cos(\pi j/N_v)$. If $c$ and $d$ are their coefficients, the stress integral is

$$
S = 2\,\mathtt{cdd}'(c-d) - 4\frac{\mathrm{LAM}}{\mathrm{Re}_8}\frac{\dot R}{R} - 6\,\mathtt{intfnu}\,\mathtt{iDRe},
$$

the same for all three models: only the local evolution of $c$ and $d$ differs (the PTT and Giesekus terms added to the UCM/Oldroyd-B closure above, from [Warnez and Johnsen (2015)](../reference/references.md)).

The weights are

$$
\mathtt{cdd}(n) = 2L_v\int_{-1}^{1} \frac{T_n(x)-1}{(L_v(1+x)+1-x)(1-x)}\,dx.
$$

One Gauss–Legendre rule with $N_v+64$ points in $\theta=\arccos(x)$ evaluates them, where the integrands are smooth. The local field supplies the polymer part of medium stress power. See [reported comparisons](../verification/convergence.md#reported-spectral-stress-comparisons) for historical quadrature, timing and resolution measurements.

Stress 5, 6, 7 and 8 keep their integral $Z_1$ and use the same basis for the local field when medium heat is on. For stress 9 without that field, only $Z_1$ and $Z_2$ are available; the [compatibility rules](../reference/compatibility.md) explain why these cannot supply local heating. Stress 10 and 11 always run with the field (`spectral = 1` is forced), so this restriction does not apply to them.

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

The contribution vanishes at `DRe = 0`. For the Maxwell/Zener/UCM family (stress 5–11) it acts on the solvent share and also vanishes at `LAM = 0`. Its common factor is `iDRe = 1/DRe` for the Kelvin-Voigt family (stress 1–4), `LAM/DRe` for stress 5–11, and zero when inactive.

`imr.physics.viscosity` reduces this radial non-Newtonian stress to two integrals. One has a closed form; the other uses a fixed graded Gauss–Legendre rule. With $\dot\gamma=\dot\gamma_R(R/r)^3$, $\dot\gamma_R=-2\dot R/R$, and the substitution $u=R/r$, the integrals over $r$ in $[R,\infty)$ become integrals over $u$ in $[0,1]$, with $x_R=\lambda|\dot\gamma_R|$:

$$
I_1=\int_R^\infty \frac{\dot\gamma}{r}f\,dr=\dot\gamma_R\int_0^1 u^2 f(x_R u^3)\,du,
$$

$$
I_2=\int_R^\infty \frac{f+xf'}{r^4}\,dr
=R^{-3}\int_0^1 u^2\,(f+xf')(x_R u^3)\,du
=\frac{f(x_R)}{3R^3}.
$$

$I_2$ is exact: with $x=x_R u^3$, $u^2\,du=\dfrac{dx}{3x_R}$, and $f+xf'=\dfrac{d(xf)}{dx}$, so the integral is

$$
\frac{1}{3R^3}\int_0^{x_R} \frac{d(xf)}{dx}\,dx = \frac{x_R f(x_R)}{3x_R R^3} = \frac{f(x_R)}{3R^3}.
$$

$I_1$ uses a fixed Gauss–Legendre rule: 10 points on each of the 15 panels

$$
[0,2^{-14}],\ [2^{-14},2^{-13}],\ \ldots,\ [1/2,1],
$$

150 nodes, computed once. The panels shrink toward $u=0$ because at high shear rates the integrand changes on the scale $u\sim x_R^{-1/3}$. See [reported comparisons](../verification/convergence.md) and [compatibility](../reference/compatibility.md).
