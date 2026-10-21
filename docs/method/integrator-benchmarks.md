# Integrator comparison

`s_compare_integrators.m` compares the stiff time integrators of the nodal solver `f_imr_nodal` on six cases, from a polytropic Rayleigh–Plesset bubble to mass transfer with Chebyshev differentiation:

| Integrator | `method` | Notes |
|------------|----------|-------|
| `ode23tb` | 23 (default) | MATLAB code |
| `ode15s` | 15 | MATLAB code |
| CVODES | 3 or `'cvodes'` | Compiled SUNDIALS code through the `ode` object (`Solver = 'cvodesstiff'`), MATLAB R2024a or later |

Each integrator runs with two Jacobians:

- the analytic Jacobian (`jacobian = 2`, the default)
- the finite-difference Jacobian from one batched right-hand-side call (`jacobian = 3`)

## Running it

In MATLAB, from this folder:

```matlab
s_compare_integrators
```

The script adds `../../src` to the path itself.

For each case and integrator it prints:

- the error: the largest `|R - R_ref|` over 2001 output times, divided by `max |R_ref|`
- the wall time: the fastest of `nrep = 3` runs, after one warm-up run
- the speedup over `ode23tb` with the analytic Jacobian

The reference `R_ref` is `ode15s` at tolerance `1e-12` (with `ode23tb` as the fallback). The compared runs use `1e-8`, the default. A summary gives each integrator's geometric-mean speedup and worst error. The results go to `integrators.mat` and `integrators.csv`.

Without the `ode` object, the CVODES runs are skipped with a message.

## What to expect

CVODES is compiled, but it still calls the model's right-hand side and Jacobian, which are MATLAB code. It saves the integrator's own work: step control, error estimates and the Newton iterations.

The batched Jacobian costs about 3 right-hand-side calls. The analytic one costs 25–40. See the batched right-hand side section of the `f_imr_nodal` documentation for the measurements.

The error column shows whether a faster combination is also as accurate. Differences at the level of the tolerance (`1e-8` relative per step, a few `1e-7` to `1e-6` accumulated over a run) are expected.
