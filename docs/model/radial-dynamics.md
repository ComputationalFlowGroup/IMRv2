# Radial dynamics

`src/+imr/+physics/radial_acceleration.m` evaluates the wall acceleration for six choices of `radial`. The selector changes the liquid compressibility model and how pressure enters the radial equation.

| `radial` | Equation | Liquid description |
|---|---|---|
| 1 | Rayleigh–Plesset | Incompressible |
| 2 | Keller–Miksis, pressure form | Constant far-field sound speed |
| 3 | Keller–Miksis, enthalpy form | Tait EoS |
| 4 | Gilmore | Tait EoS, local sound speed |
| 5 | Keller–Miksis, enthalpy form | Mie–Grüneisen EoS |
| 6 | Gilmore | Mie–Grüneisen EoS, local sound speed |

An equation of state (EoS) relates pressure to the liquid density and acoustic quantities. The enthalpy forms use a pressure-only relation at the far-field temperature; they cannot be combined with medium heat transfer in this implementation. Radial 5–6 are further limited to the polytropic bubble. See [compatibility](../reference/compatibility.md) for complete requirements.

`cstar`, `gama`, `nstate` and `hugoniot_s` are defined in [parameters](parameters.md#acoustic-options). Choosing radial 5 or 6 does not supply a validated material parameter set. Material profiles and radial EoS coefficients have distinct roles; see [fluid profiles](fluids.md).

All quantities below are dimensionless: radius is scaled by $R_0$, pressure by $P_\infty$, density by $\rho_\infty$, velocity by $U_c=\sqrt{P_\infty/\rho_\infty}$, time by $R_0/U_c$, and liquid enthalpy by $U_c^2$. A dot denotes a derivative with respect to that dimensionless time. See [parameters](parameters.md) for SI conversion. The implementation notes for the Tait and Mie–Grüneisen EoS evaluation live in the [equation-of-state method](../method/equation-of-state.md).

## Pressure forms

| Symbol | Meaning | Code |
|---|---|---|
| *R*, *U* | Radius and wall velocity | `R`, `Rdot` |
| *P*, *p_w*, *p_f* | Total bubble, liquid-wall and absolute far-field pressures | `P`, wall closure, `1 + Pf8` |
| *f*, *ḟ* | Far-field perturbation and time derivative | `Pf8`, `Pf8dot` |
| *S*, *Ṡ₀* | Signed stress integral and acceleration-independent rate | `S`, `Sdot` |
| *W*, *A_S* | Surface-tension factor and acceleration-dependent stress coefficient | `iWe`, expression below |
| *C*, *ρ*, *h*, *H* | Sound speed, liquid density, reference enthalpy and wall/far-field difference | Radial closure |

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

`JdotA` accounts for instantaneous Newtonian/solvent stress, while `ddintfnu` is the acceleration coefficient of the non-Newtonian viscosity integral's derivative. Their definitions and the solvent share in `iDRe` are described in [stress and viscosity](stress-viscosity.md).

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

The terms containing $1/C$ represent finite sound-speed effects; sending $C$ to infinity recovers Rayleigh–Plesset. Moving the acceleration-dependent stress term to the denominator avoids using an unknown acceleration in the numerator. [Fluid preparation](../reference/nondimensionalize.md#fluid-profiles) determines ambient density and sound speed once; this branch does not evaluate a wall EoS during integration.

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

The Tait relation (radial 3–4) and the existing Hugoniot-based Mie–Grüneisen relation (radial 5–6) that supply $\rho(p)$, $h(p)$ and $C(p)$ above are stated in the [equation-of-state method](../method/equation-of-state.md), along with their numerically robust evaluation.

See [compatibility](../reference/compatibility.md) for thermal restrictions, [numerical checks](../verification/numerics.md) for forcing regressions, and [reported comparisons](../verification/convergence.md#reported-radial-enthalpy-comparison) for historical unforced measurements.
