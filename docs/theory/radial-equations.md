# Radial equations and equation of state

`src/+imr/+physics/radial_acceleration.m` computes wall acceleration. All quantities below are dimensionless: radius is scaled by $R_0$, pressure by $P_\infty$, density by $\rho_\infty$, velocity by $U_c=\sqrt{P_\infty/\rho_\infty}$ time by $R_0/U_c$, and liquid enthalpy by $U_c^2$. A dot denotes a derivative with respect to that dimensionless time. See [parameters](../model/parameters.md) for SI conversion.

## Pressure forms

| Symbol | Meaning | Code |
|---|---|---|
| $R$, $U$ | Radius and wall velocity | `R`, `Rdot` |
| $P$, $p_w$, $p_f$ | Total bubble, liquid-wall and absolute far-field pressures | `P`, wall closure, `1 + Pf8` |
| $f$, $\dot f$ | Far-field perturbation and time derivative | `Pf8`, `Pf8dot` |
| $S$, $\dot S_0$ | Signed stress integral and acceleration-independent rate | `S`, `Sdot` |
| $W$, $A_S$ | Surface-tension factor and acceleration-dependent stress coefficient | `iWe`, expression below |
| $C$, $\rho$, $h$, $H$ | Sound speed, liquid density, reference enthalpy and wall/far-field difference | Radial closure |

Radius must remain positive and the acceleration denominator must be nonsingular. Pressure-form Keller–Miksis keeps the ambient sound speed fixed and retains first-order wall-Mach effects. Enthalpy models additionally require admissible pressures and positive density/sound speed on their barotrope. See [compatibility](../reference/compatibility.md) for thermal restrictions.

Let $U=\dot R$, let $P$ be total bubble pressure (gas plus enabled vapor), and let $f(t)$ be the imposed pressure above ambient (`Pf8`). Then the far-field absolute pressure is $p_f=1+f$. Define the liquid pressure at the wall and the driving difference by

$$
p_w=P-\frac{W}{R}+S,\qquad Q=p_w-p_f,\qquad W=\mathrm{We}^{-1}.
$$

Here $W/R$ is the Laplace-pressure drop and $S$ is the integrated material-stress contribution, including its sign convention. In an inviscid liquid with no surface tension or forcing, $Q=P-1$. Rayleigh–Plesset (`radial = 1`) is

$$
\ddot R=\frac{Q-\tfrac32 U^2}{R}.
$$

The stress derivative contains terms proportional to acceleration. The code's `Sdot` is the remaining term, denoted $\dot S_0$ here; the coefficient of $-\ddot R$ is collected in $A_S$:

$$
\frac{dS}{dt}=\dot S_0-A_S\ddot R,\qquad
A_S=\frac{\mathtt{JdotA}}{R}
+6\,\mathtt{ddintfnu}\,\mathtt{iDRe}.
$$

`JdotA` accounts for instantaneous Newtonian/solvent stress, while `ddintfnu` is the acceleration coefficient of the non-Newtonian viscosity integral's derivative. Their definitions and the solvent share in `iDRe` are described in [stress and viscosity](../model/stress-viscosity.md).

Pressure-form Keller–Miksis (`radial = 2`) uses the fixed ambient sound speed $C=\mathtt{cstar}$:

$$
\ddot R=
\frac{
 (1+U/C)Q+
 \dfrac{R}{C}\left(\dot P+\dfrac{WU}{R^2}+\dot S_0-\dot f\right)
 -\dfrac32\left(1-\dfrac{U}{3C}\right)U^2
}{
 (1-U/C)R+\dfrac{R A_S}{C}
}.
$$

The terms containing $1/C$ represent finite sound-speed effects; sending $C$ to infinity recovers Rayleigh–Plesset. Moving the acceleration-dependent stress term to the denominator avoids using an unknown acceleration in the numerator. [Fluid preparation](../api/f_dim2nondim.md#fluid-profiles) determines ambient density and sound speed once; this branch does not evaluate a wall EoS during integration.

## Enthalpy forms

Radial 3–6 use the **same barotrope at the wall and far field**. For its dimensionless density $\rho(p)$ define

$$
h(p)=\int_1^p\frac{dq}{\rho(q)},\qquad
H=h(p_w)-h(p_f),\qquad h'(p)=\frac{1}{\rho(p)}.
$$

This pressure-integrated liquid enthalpy is distinct from the temperature-only caloric slice returned by a material profile. Its time derivative, after separating acceleration-dependent stress, is

$$
\dot H_0=
\frac{\dot P+WU/R^2+\dot S_0}{\rho(p_w)}
-\frac{\dot f}{\rho(p_f)}.
$$

The resulting acceleration is

$$
\ddot R=
\frac{(1+U/C)H+(R/C)\dot H_0
-\tfrac32(1-U/(3C))U^2}
{(1-U/C)R+R A_S/(C\rho(p_w))}.
$$

Keller–Miksis enthalpy (`radial = 3, 5`) uses `cstar`; Gilmore (`radial = 4, 6`) uses the local wall sound speed. Matching wall and far-field pressures and derivatives gives zero acceleration at rest. The previous implementation subtracted $f$ and $\dot f$ directly, approximating the far-field density by one. The corrected code evaluates both enthalpy terms with the selected EoS; driven enthalpy histories can change. Unforced calculations retain the same reference state.

### Tait relation

For $B=\mathtt{gama}$, $n=\mathtt{nstate}$, $a=1+B$ and $m=(n-1)/n$, radial 3–4 use

$$
\rho(p)=\left(\frac{p+B}{a}\right)^{1/n},\qquad
h(p)=\frac{a}{m}\left[\left(\frac{p+B}{a}\right)^m-1\right],\qquad
C(p)^2=\frac{n(p+B)}{\rho(p)}.
$$

The pressure shift must satisfy $p+B>0$. At $n=1$, the finite limit is $h(p)=a\ln((p+B)/a)$. The code uses this logarithmic limit exactly and evaluates the other exponents with `log1p` and `expm1` to preserve small pressure changes and exponents near one. In the source, `sam = a`, `no = m`, `hB = h(p_w)`, `hF = h(p_f)` and `hH = 1/rho(p_w)`.

### Existing Hugoniot-based relation

For radial 5–6, write $s=\mathtt{hugoniot\_s}$, $b=(n-1)/2=\mathtt{nog}$ and $\mu=\rho-1$. The implemented relation is

$$
\frac{p}{C_*^2}=g(\mu)
=\frac{\mu(1+b\mu)}{(1-s\mu)^2},\qquad C_*=\mathtt{cstar}.
$$

The code solves the quadratic for the small physical root, which tends to $p/C_*^2$ at low pressure. Integration by parts gives

$$
h(p)=\left[\frac{q}{1+\mu(q)}\right]_1^p
+C_*^2\int_{\mu(1)}^{\mu(p)}\frac{g(\mu)}{(1+\mu)^2}\,d\mu.
$$

Partial fractions at $\mu=-1$ and $\mu=1/s$ give the closed form in the nested `f_mie_gruneisen_eos_scalar`; `log1p` differences reduce cancellation near the reference state. These branches implement a pressure-based liquid model, rather than a complete arbitrary-temperature thermodynamic EoS. They are distinct from the Murnaghan property family.

See [compatibility](../reference/compatibility.md) for thermal restrictions, [numerical checks](../validation/numerics.md) for forcing regressions, and [reported comparisons](../validation/convergence.md#reported-radial-enthalpy-comparison) for historical unforced measurements.
