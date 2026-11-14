# Equation-of-state evaluation

The enthalpy radial forms (`radial = 3–6`, see [radial dynamics](../model/radial-dynamics.md#enthalpy-forms)) need the barotrope $\rho(p)$, its enthalpy integral $h(p)$ and sound speed $C(p)$ evaluated robustly, including at small pressure changes and near reference exponents. This page states the two implemented relations and the numerical care their evaluation takes.

### Tait relation

For $B=\mathtt{gama}$, $n=\mathtt{nstate}$, $a=1+B$ and $m=(n-1)/n$, radial 3–4 use

$$
\rho(p)=\left(\frac{p+B}{a}\right)^{1/n},\qquad
h(p)=\frac{a}{m}\left[\left(\frac{p+B}{a}\right)^m-1\right],\qquad
C(p)^2=\frac{n(p+B)}{\rho(p)}.
$$

The pressure shift must satisfy $p+B>0$. At $n=1$, the finite limit is $h(p)=a\ln((p+B)/a)$. The code uses this logarithmic limit exactly. It evaluates the other exponents with `log1p` and `expm1`. This preserves small pressure changes and exponents near one. In the source, `sam = a`, `no = m`, `hB = h(p_w)`, `hF = h(p_f)` and `hH = 1/rho(p_w)`.

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

Partial fractions at $\mu=-1$ and $\mu=1/s$ give the closed form in the nested `f_mie_gruneisen_eos_scalar`. `log1p` differences reduce cancellation near the reference state. These branches implement a pressure-based liquid model, rather than a complete arbitrary-temperature thermodynamic EoS. They are distinct from the Murnaghan property family.

See [compatibility](../reference/compatibility.md) for thermal restrictions, [numerical checks](../verification/numerics.md) for forcing regressions, and [reported comparisons](../verification/convergence.md#reported-radial-enthalpy-comparison) for historical unforced measurements.
