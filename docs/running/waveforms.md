# Pressure waveforms

`imr.physics.far_field_pressure(t,{[om,ee,tw,dt,mn,wave_type],poly,dpoly})` returns a pressure **perturbation** $f$ and its time derivative. Absolute far-field pressure is $p_f=1+f$, in units of ambient pressure. Time is in units of $t_0=R_0/U_c$; derivatives are with respect to that dimensionless time. See the [radial equations](../theory/radial-equations.md#pressure-forms) for how $f$ and $\dot f$ drive acceleration.

| Symbol | Input | Meaning |
|---|---|---|
| *A* | `ee` | Signed dimensionless amplitude |
| *ω* | `om` | Angular frequency times *t₀* |
| *w* | `tw` | Gaussian width or step ending time |
| *d* | `dt` | Gaussian/histotripsy center time |
| *m* | `mn` | Histotripsy exponent |

**Sign convention:** Positive `ee` always drives compressive waves (positive far-field pressure rise); negative `ee` drives tensile (rarefaction) waves. Physically: tension pulls material apart (negative stress), and compression pushes it together (positive stress). This convention is consistent across all waveform types. `imr.nondimensionalize` scales SI `pa`, `omega`, `tw` and `dt` before preparation. Scalar and array times preserve their input shape.

## Analytic choices

### Constant pressure: type 0

$$
f(t)=A,\qquad \dot f(t)=0.
$$

The legacy code calls this an impulse, but it is a constant pressure offset over the simulation, not a Dirac impulse. Both factories bypass the waveform call for this case.

### Gaussian pulse: type 1

For positive width $w$, set $z=(t-d)/w$:

$$
f(t)=A e^{-z^2},\qquad
\dot f(t)=-\frac{2A(t-d)}{w^2}e^{-z^2}.
$$

Positive amplitude gives a compression pulse (pressure rise) centered at $d$; its derivative changes sign there. Negative amplitude gives a rarefaction (tensile) pulse. The pulse has infinite mathematical support and decays smoothly. The [Gaussian suite](../validation/gaussian.md) checks return toward equilibrium after forcing.

### Histotripsy pulse: type 2

For $\omega>0$, $m>0$ and $\phi=\omega(t-d)/2$, inside $|t-d|<\pi/\omega$:

$$
f(t)=A\cos^{2m}\phi,\qquad
\dot f(t)=-A m\omega\cos^{2m-1}\phi\sin\phi.
$$

The function and returned derivative are zero outside the support and at its endpoints. At $m=0$, pressure instead equals $A$ on the **closed** support interval, including both endpoints, and zero elsewhere; the returned derivative is zero. This convention avoids evaluating `0*Inf` and is tested in the [numerical suite](../validation/numerics.md).

Endpoint derivatives are a pointwise numerical convention. A rectangular pulse has distributional impulses at its jumps; these are not represented. For $0<m<1/2$ the interior derivative has an unbounded one-sided endpoint limit, while $m=1/2$ gives a derivative jump. Use $m>1/2$ for a continuously differentiable pulse. The half-angle form avoids cancellation near endpoints for fractional exponents.

### Finite step: type 3

For $w\geq0$:

$$
f(t)=\begin{cases}A,&t<w,\\0,&t\geq w,\end{cases}
\qquad \dot f(t)=0.
$$

Positive amplitude $A$ (positive `ee`) gives a constant compressive pressure rise before the step at $t=w$; negative amplitude gives a constant tensile pressure drop (rarefaction). The derivative omits the distributional impulse at $w$. With `tw = 0` the perturbation is zero for every allowed nonnegative simulation time. For a discontinuous drive, choose the requested time span and step resolution deliberately; a zero returned derivative does not make the waveform smooth.

## Measured waveforms

The three measured driving-pressure waveforms (`wave_type = -1, -2, -3`) are from:

E. Hersey, M. Rodriguez and E. Johnsen, "Dynamics of an oscillating microbubble in a blood-like Carreau fluid," *J. Acoust. Soc. Am.* **153**(3), 1836–1845 (2023). <https://doi.org/10.1121/10.0017342>

Please cite this paper if you use these waveforms.

## Files

| File | Contents |
|---|---|
| `patterson_waveform_data.mat` | The three waveforms, stored as dimensionless time and pressure columns: mostly linear (ML), moderately nonlinear (MN) and highly nonlinear (HN) |
| `r_patterson_data.m` | Fits a cubic spline to each waveform and its derivative, plots them, and saves them as `ml.mat`, `mn.mat` and `hn.mat`. It uses `fnder` from the Curve Fitting Toolbox. |

Preparation loads the spline fits from `src/`, through `imr.configuration.prepare_case`:

| `wave_type` | File | Waveform |
|---|---|---|
| -1 | `src/data/d_hn.mat` | Highly nonlinear |
| -2 | `src/data/d_mn.mat` | Moderately nonlinear |
| -3 | `src/data/d_ml.mat` | Mostly linear |

Each file holds the spline `poly` of the dimensionless pressure and its derivative `dpoly`, over dimensionless times up to 15. [`imr.physics.far_field_pressure`](../reference/functions/f_pinfinity.md) evaluates them.

Measured choices use the stored spline directly: `ee`, `dt` and `om` do not rescale or shift it. Outside its stored breakpoints, `ppval` extrapolates; preparation warns if the requested time interval extends outside the data. This is not additional measured evidence beyond that interval.
