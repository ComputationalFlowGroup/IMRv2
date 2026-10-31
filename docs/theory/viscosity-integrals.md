# Viscosity integrals

`imr.physics.viscosity` reduces the radial non-Newtonian stress to two integrals. One has a closed form; the other uses a fixed graded Gauss–Legendre rule.

With $\dot\gamma=\dot\gamma_R(R/r)^3$, $\dot\gamma_R=-2\dot R/R$, and the substitution $u=R/r$, the integrals over $r$ in $[R,\infty)$ become integrals over $u$ in $[0,1]$, with $x_R=\lambda|\dot\gamma_R|$:

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

150 nodes, computed once. The panels shrink toward $u=0$ because at high shear rates the integrand changes on the scale $u\sim x_R^{-1/3}$. See [reported comparisons](../validation/convergence.md).
