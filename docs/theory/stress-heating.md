# Stress power and local fields

`imr.physics.stress_power` uses the full stress power, including elastic work, for medium heating.
The following derivation states that convention and the local Zener field chosen to reproduce the stress integral.

`imr.physics.stress_power` computes the medium heating by the stress, $\mathtt{taugradu}=\mathrm{Br}\,\Phi$, when `medtherm = 1`. The medium temperature equation adds it: $\dot T_m = \text{advection}+\text{diffusion}+\mathtt{taugradu}$. Both solvers call it.

The heating is the full stress power,

$$
\Phi = \boldsymbol\tau:\mathbf D.
$$

This is the entropy-elastic (rubber-like) convention: elastic work is released as heat under loading and absorbed under unloading.

With the local strain rate

$$
a=\frac{R^2\dot R}{r^3}=\frac{\dot R}{R}\,y^{-3},\qquad y=\frac{r}{R},
$$

and $s=\tau_{rr}-\tau_{\theta\theta}$, the stress power is

$$
\boldsymbol\tau:\mathbf D = -2as.
$$

Only this normal-stress difference enters; the trace of the stress does not.

## Heating per model

| `stress` | $\mathtt{taugradu}/\mathrm{Br}$ |
|---|---|
| 0 | $0$ |
| 1 | $12\mu_{KV}a^2-2as_e$, neo-Hookean $s_e$ |
| 2 | the same with the quadratic $s_e$ |
| 3 | the same with the 2nd-order Fung $s_e$ |
| 4 | the same with the Mooney-Rivlin $s_e$ |
| 5, 6, 7, 8 | $12\mu_s a^2-2as_p$, with the Zener polymer field $s_p$ of [`imr.physics.zener_field`](../reference/functions/f_zener_field.md) |
| 9, 10, 11 | $12\mu_s a^2-2a(\tau_{rr}-\tau_{\theta\theta})$, from the spectral stress field (UCM/Oldroyd-B, PTT or Giesekus) |

**Viscosities**, with $f_\nu$ the non-Newtonian factor at the medium points:

- $\mu_{KV}=1/\mathrm{Re}_8+f_\nu\,\mathtt{iDRe}$: stress 1–4, where the whole viscosity acts instantly.
- $\mu_s=\mathrm{LAM}/\mathrm{Re}_8+f_\nu\,\mathtt{iDRe}$: stress 5–11, where only the solvent part acts instantly.

`iDRe` is $1/\mathtt{DRe}$ for stress 1–4 and $\mathrm{LAM}/\mathtt{DRe}$ for stress 5–11, the same factor as in the radial equation (set in `imr.configuration.prepare_case`).

**Elastic stress difference.** With the hoop stretch $\lambda=r/r_0$, where $r_0^3=r^3-R^3+R_{eq}^3$, and $x=\lambda^{-4}+2\lambda^2-3$ (the invariant term $I_1-3$):

- neo-Hookean (1, 5): $s_e=\dfrac{\lambda^{-4}-\lambda^2}{\mathrm{Ca}}$;
- quadratic (2, 6): the same times $1+\alpha_x x$, the first-order term of the series $\exp(\alpha_x x)$;
- 2nd-order Fung (3, 7): the same times $1+\alpha_x x+\dfrac{\alpha_x^2}{2}x^2$, the next term of the same series -- unlike stress 3's own wall (radial-equation) formula, which was cross-checked against an independent codebase, this local field is derived here (the natural next term of the series) and verified only by its own consistency: it reproduces the independently-sourced wall formula exactly under $S = 2\int_1^\infty \frac{s_e}{y}\,dy$, and its $R_{eq}/R$-derivative matches the wall formula's own `Sdot` (see `imr.physics.material_stress`);
- Mooney-Rivlin (4, 8): $s_e=\dfrac{1}{\mathrm{Ca}}\Big[(\lambda^{-4}-\lambda^2)-C_2(\lambda^4-\lambda^{-2})\Big]$, the neo-Hookean term plus one new term scaled by the second-modulus ratio $C_2$; this field and its wall integral are both derived here from the two-invariant Mooney-Rivlin strain energy, with no external source, and cross-checked the same way (integral and derivative consistency).

In the code, `lam^2 = yT2./x2` and `lam^-4 = x2.^2.*iyT4`, with `x2 = (yT3 - 1 + (Req/R)^3).^(2/3)`; `lam^4 = (yT2.*ix2).^2` and `lam^-2 = x2./yT2` for the Mooney-Rivlin term.

## Total heating

Integrating over the medium gives

$$
4\pi R^3\int_1^\infty y^2\Phi\,dy = -4\pi R^2\dot R\,S,
$$

with $S$ the stress integral of the radial equation. This holds for every model. For stress 5, 6, 7 and 8 it equals

$$
16\pi\frac{\mathrm{LAM}}{\mathrm{Re}_8}R\dot R^2+24\pi R^2\dot R\,\mathtt{intfnu}\,\mathtt{iDRe}-\frac{4\pi\dot R Z_1}{R}.
$$

## Equation

On the grid $y=r/R$, the field obeys

$$
\left.\frac{\partial s_p}{\partial t}\right|_y = -\dot y\,\frac{\partial s_p}{\partial y} - \frac{s_p-s_{KV}}{\mathrm{De}} - 3\frac{\dot R}{R}\left(1-y^{-3}\right)s_p,
$$

$$
\dot y = \frac{\dot R}{R}\left(y^{-2}-y\right),
$$

$$
s_{KV} = s_e - 6\eta_p a,\qquad a=\frac{\dot R}{R}y^{-3},\qquad \eta_p=\frac{1-\mathrm{LAM}}{\mathrm{Re}_8}.
$$

- $s_e$ is the neo-Hookean elastic stress difference for stress 5, the quadratic one for stress 6, the 2nd-order Fung one for stress 7, or the Mooney-Rivlin one for stress 8 (see [`imr.physics.stress_power`](../reference/functions/f_stress_dissipation.md)).
- The caller computes the advection $-\dot y\,\partial s_p/\partial y$ on its own grid.
- **Boundary conditions:** none at the wall, where $\dot y=0$. $s_p$ is 0 at infinity.
- **Initial condition:** $s_p=0$ from rest.

## Why this law

It is the local law whose integral reproduces the solver's stress state exactly:

$$
2\int_1^\infty \frac{s_p}{y}\,dy = \frac{Z_1}{R^3},
$$

where $Z_1$ obeys

$$
\mathrm{De}\,\dot Z_1 + Z_1 = Z_e - 4\eta_p R^2\dot R
$$

with $Z_e = R^3 S_e$.

**Proof.** Written for the weighted field $w=y^{-3}s_p$, the law reads $\dfrac{Dw}{Dt} = -\dfrac{w-y^{-3}s_{KV}}{\mathrm{De}}$ along a material point. From rest its solution is

$$
y^{-3}s_p(t) = \int_0^t K(t-t')\,y'^{-3}s_{KV}(t')\,dt',\qquad K(t)=\frac{e^{-t/\mathrm{De}}}{\mathrm{De}},
$$

where primes denote the same material point at time $t'$.

Incompressibility gives $r^3-r'^3=R^3-R'^3$, so $dr/r = r'^2\,dr'/r^3$. Hence

$$
2\int_R^\infty \frac{s_p}{r}\,dr = \frac{1}{R^3}\int_0^t K\,R'^3\left[2\int_{R'}^\infty \frac{s_{KV}}{r'}\,dr'\right]dt' = \frac{1}{R^3}\int_0^t K\,R'^3\,S_{KV}(t')\,dt',
$$

with $S_{KV}=S_e-4\eta_p\dot R/R$, the stress integral of $s_{KV}$. This is $Z_1/R^3$, the solution of the $Z_1$ equation from rest.

The $Z_1$ equation fixes only this weighted integral. The law above is the simplest local field consistent with it, not the only one. It sets the radial distribution of the heating, but neither $R(t)$ nor the total heating.

So:

- $Z_1$ still drives $R(t)$, unchanged.
- The field only distributes the heating in $r$.
- The total heating is fixed by $Z_1$.

## In the solvers

- **`imr.solve_spectral`.** $s_p$ is expanded in the stress 9 basis $T_n(z)-1$ with `Nv` modes, stored after $Z_1$, and evaluated at the medium points with `ZZT`. After a run, the solver prints the largest mismatch of $2\,\mathtt{cdd}'e$ against $Z_1/R^3$ as a check.
- **`imr.solve_nodal`.** $s_p$ is stored at the medium nodes `yT` after $Z_1$: the first node is the wall, and the last, at infinity, is held at 0. The advection reuses the $\dot R$ part of the medium temperature advection, `xi1_sq./(Lt*R).*(Rdot./yT2.*om_yT3/2).*(dmatrix_Tm*sp)`, which equals $-\dot y\,\partial s_p/\partial y$.
- **Collapse (nodal only).** The history of the field before $t=0$ is not computed. As a stopgap, the initial field takes the shape of the elastic stress difference of the initial radius, scaled so that the integral identity holds for the initial $Z_1$ from `imr.physics.initial_stress`. Without elasticity ($\mathrm{Ca}=\infty$), it takes the viscous shape $y^{-3}$ instead.

## Spectral stress field

For stress 9, 10 and 11, `imr.solve_spectral` expands $\tau_{rr}$ and $\tau_{\theta\theta}$ in $T_n(z)-1$, $n=1,\ldots,N_v$, with

$$
z = 1 - \frac{2}{1 + (r/R - 1)/L_v}
$$

and collocation points $z_j=\cos(\pi j/N_v)$. If $c$ and $d$ are their coefficients, the stress integral is

$$
S = 2\,\mathtt{cdd}'(c-d) - 4\frac{\mathrm{LAM}}{\mathrm{Re}_8}\frac{\dot R}{R} - 6\,\mathtt{intfnu}\,\mathtt{iDRe},
$$

the same for all three models: only the local evolution of $c$ and $d$ differs (see [stress and viscosity](../model/stress-viscosity.md#ptt-and-giesekus-stress-10-11) for the PTT and Giesekus terms added to the UCM/Oldroyd-B closure, from [Warnez and Johnsen (2015)](../reference/references.md)).

The weights are

$$
\mathtt{cdd}(n) = 2L_v\int_{-1}^{1} \frac{T_n(x)-1}{(L_v(1+x)+1-x)(1-x)}\,dx.
$$

One Gauss–Legendre rule with $N_v+64$ points in $\theta=\arccos(x)$ evaluates them, where the integrands are smooth. The local field supplies the polymer part of medium stress power. See [reported comparisons](../validation/convergence.md#reported-spectral-stress-comparisons) for historical quadrature, timing and resolution measurements.

Stress 5, 6, 7 and 8 keep their integral $Z_1$ and use the same basis for the local field when medium heat is on. For stress 9 without that field, only $Z_1$ and $Z_2$ are available; the [compatibility rules](../reference/compatibility.md) explain why these cannot supply local heating. Stress 10 and 11 always run with the field (`spectral = 1` is forced), so this restriction does not apply to them.
