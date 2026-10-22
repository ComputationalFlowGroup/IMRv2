# Convergence and reported comparisons

This page preserves numerical measurements reported in the previous function documentation.
They were not rerun for this documentation change. Many lack a bundled reproduction driver or captured log; treat their digits and timings as historical reports, not current acceptance criteria.
The executable assertions live in the [test suite pages](tests.md).

## Reported measurements: `imr.solve_nodal`

On thermal cases with `Req = 0.5` and `t` in [0, 3]:

- `scheme = 0` gave bit-identical `R(t)` to `f_imr_fd` over 98 cases: radial 1–4, stress 0–5, polytropic, bubble, bubble and medium, and full mass-transfer models, plus a non-Newtonian and a Gaussian-waveform case. This held before the bracketed mass-transfer wall closure; with it, the cases without mass transfer stay bit-identical and the mass-transfer cases differ at the integrator tolerance.
- `scheme = 1` and `imr.solve_spectral` converge toward each other with resolution over 32 cases (radial 1–4, stress 0/1/3/5, with and without medium heat transfer). The largest relative deviation in `R(t)` over the cases is `2.8e-3`, `4.8e-4` and `2.6e-4` at `Nt = Mt = 12, 20, 30`; the median is `1.4e-4`, `1.0e-5` and `3.4e-6`.
- Against Richardson-extrapolated finite differences (`Nt = Mt = 200, 400`) for Keller–Miksis with neo-Hookean Kelvin–Voigt stress, `scheme = 1` at `Nt = Mt = 25` is within `1e-4` with bubble heat transfer only, where finite differences at 200 points are `2.4e-3` away. With medium heat transfer it reaches `2e-4` at 40 points, the same as `imr.solve_spectral`.
- With mass transfer (Keller–Miksis, neo-Hookean Kelvin–Voigt, bubble and medium heat transfer), `scheme = 1` converges with resolution: against its own solution at 80 points the largest deviation is `1.2e-2`, `3.6e-3`, `2.6e-4` and `5.5e-5` at 20, 30, 40 and 60 points for `Req = 0.5`, and `4.9e-2`, `2.0e-2`, `2.5e-3` and `1.1e-4` for `Req = 0.4`. At 40 points it is within `1e-3` of Richardson-extrapolated finite differences (`Nt = Mt = 200, 400`), the uncertainty of that reference.
- The analytic Jacobian matches a central-difference Jacobian to a scaled error of `3e-9` for `scheme = 0` and `4e-10` for `scheme = 1`, with and without mass transfer.

## Reported measurements: `imr.physics.initial_stress`

- Against the previous `fminsearch` shooting, over the default collapse parameters with `Req` 0.1–0.5 and `De` 0.01–1, the stress agrees to `3e-8` relative. At `Req = 0.05` no subsonic start reaches `R = 1` (the maximum radius stays at 0.93); the previous version returned the stress at that non-root without notice.
- Integrates with `ode23tb`.

## Reported measurements: `imr.numerics.jacobian_nodal`

Against a central-difference Jacobian of the right-hand side over radial 1–4, stress 0/1/3/5, viscosity models 0/1/5, with and without medium heat and mass transfer, the maximum scaled difference is below `3e-9` for finite differences and `4e-10` for Chebyshev.

## Reported measurements: `imr.numerics.jacobian_spectral`

Against a central-difference Jacobian over radial 1–4, stress 0/1/3/5, Newtonian and Carreau viscosity, polytropic, cold-medium and full thermal models, and the spectral stress representation (102 configurations), the maximum scaled difference is `3e-7`.

## Reported measurements: `imr.physics.stress_power`

- **Energy budget.** The total heating from quadrature on the solver grids matches `-4*pi*R^2*Rdot*S` to:
    - `1e-8`–`5e-5` (relative) with Chebyshev and spectral grids;
    - `5e-4` with finite differences at `Mt = 80` (second-order quadrature error).

    This held for stress 1–5, Newtonian and Carreau.
- **Kelvin–Voigt limit.** At the same state, the stress 3 (4) heating tends to the stress 1 (2) heating as `De -> 0`: the relative difference is `4e-2`, `4e-3` and `4e-4` for `De = 1e-2`, `1e-3` and `1e-4`.

## Reported measurements: `imr.physics.zener_field`

The largest error in `2*int sp/y dy - Z1/R^3` over a run, relative to `max |Z1/R^3|` (Keller–Miksis, `Req = 0.5`):

| Mode | Stress 3 | Stress 4 |
|---|---|---|
| Nodal, finite differences, `Mt = 40, 80, 160` | `5.8e-4`, `1.4e-4`, `3.5e-5` | `1.4e-3`, `3.5e-4`, `8.9e-5` |
| Nodal, Chebyshev, `Mt = 12, 20, 30` | `2.0e-5`, `4.9e-7`, `7.7e-8` | `9.1e-5`, `7.4e-6`, `4.7e-7` |
| Spectral, `Nv = 8, 16, 32` | `1.9e-3`, `3.6e-5`, `5.9e-7` | `7.6e-3`, `6.2e-4`, `1.5e-5` |

That is second-order convergence for finite differences and spectral convergence otherwise.

## Reported measurements: `imr.numerics.differentiation`

`scheme = 0` reproduces the solver's operator bitwise for `nodes` in {10, 100, 200}, both orders and both grids. The Chebyshev interior operators differentiate `cos(pi*y/2)` and `y^2` to 1e-11 or better for 8–32 nodes and annihilate constants to 1e-13; the exterior operators converge spectrally (`exp(xi)`: 6e-6 at 8 nodes, 4e-14 at 16 nodes).

## Reported cross-mode regression agreement

In the cross-mode report, Chebyshev and the spectral solver agree with finite differences to about 2% on the cases without mass transfer. The mass-transfer cases differ by up to 11% between the two nodal schemes. At `Nt = Mt = 20` (Chebyshev) and 70 (finite differences), neither nodal scheme is resolved for mass transfer: the convergence study in the `imr.solve_nodal` documentation needs about 40 Chebyshev points. The suite checks that results are reproduced, not that they are accurate, so these grids are adequate.

## Reported batched-Jacobian measurements


- Every column of a batch matches a single-column call to `4e-13` (relative), over 372 model combinations: radial 1–6, stress 0–5, four viscosity models, the thermal and mass-transfer configurations, and both schemes.
- The batched Jacobian matches the analytic one to a scaled error below `1e-5` in every thermal case.
- Over the 74 regression cases, `jacobian = 3` reproduces the `jacobian = 2` solution to within the integrator tolerance (`R(t)` differs by at most `9e-7`, mostly `1e-8` or less).
- `benchmarks/integrators/s_compare_integrators.m` compares `ode23tb`, `ode15s` and CVODES with both Jacobians against a tight-tolerance reference, in MATLAB.

## Reported viscosity quadrature comparison

Against an adaptive `integral` at `1e-13`, the rule is accurate to `1.3e-11` relative for all seven models, exponents 0.3–0.8, and `x_R` from `1e-4` to `1e9`.

Compared with the previous version, which called `integral` twice at `1e-8`, the outputs agree to `4e-11`, with one fixed graded Gauss rule in place of two adaptive integrations. A run with a non-Newtonian model takes about half the time.

## Legacy scaling data

`benchmarks/scaling` retains old `.mat` outputs and the scripts that generated related studies. The current generators explicitly say the stored data used radial 3, while the current code uses radial 2 to allow medium heat transfer. Do not describe those datasets as current-model validation or silently overwrite them.

The Slurm templates load MATLAB R2019a, below the documented supported version, and depend on site-specific `module`, `sbatch` and `matlab-threaded` commands. The plotting scripts also use a finite-resolution history named `Rexact`; that name does not make it an analytic solution. Regenerate a complete, consistent study and record its revision before interpreting new convergence numbers.

## Reported collapse setup comparison

The earlier parser page reported that the recomputed collapse equilibrium radius reproduces the input `Req` over 0.05–0.9 and `Pv` values 0.03 and 0.3. No reproduction driver or run log accompanies that statement here. The initial-pressure and equilibrium consistency should be checked for the chosen case.

## Reported radial enthalpy comparison

- Mie–Grüneisen enthalpy: with `s = hugoniot_s` and `A = p/Cstar^2 = g(mu) = mu*(1 + nog*mu)/(1 - s*mu)^2`, `h(P) = int_1^P dp/(1 + mu)` is integrated by parts, `h = [p/(1 + mu)] + Cstar^2*int g/(1 + mu)^2 dmu`, and the partial fractions of `g/(1 + mu)^2` (poles at `mu = -1` and `mu = 1/s`) give it in closed form. The differences are written so that they hold for `s = 0` and near `P = 1`. Against an adaptive `integral` at `1e-12` the closed form agrees to `2.4e-13` relative to `|P - 1|` over 396 states (`Cstar` 5–1000, `nog` 0.5–2.5, `s` 0–2, `P` from -0.5 to `1e4`). It replaces an adaptive `integral` on every right-hand-side call, and radial 5 and 6 runs take about half the time, with `R(t)` unchanged to `1e-13`.

## Reported wall-closure comparisons

**Warm start.** Consecutive calls come from nearby states, so each call first runs Newton with the analytic derivative from the previous root. The result is accepted once a Newton step falls below `1e-10` with every iterate inside the bracket, within 6 steps. The first call of a run, and any call where Newton leaves the bracket, meets a non-finite value or does not converge, uses the bracketed solve below instead. Over 360 trajectory states the warm start was accepted in 16382 of 16392 calls, after 2.75 Newton steps on average, and its roots agree with the reference to `5e-16` in `theta_w`. The residual has exactly one sign change in `(0, 1)` at every one of those states (checked on 20001 points), so the warm start and the bracketed solve find the same root. The wall value now depends on the previous call only through round-off: in the regression cases with mass transfer, `R(t)` moves by at most `9e-7` (median `3e-8`), at the level of the integrator tolerance.

**Solve.** The bracket is fixed at `kv_w` in `[1e-13, 1 - 1e-13]`; no initial guess is used. The residual is evaluated on 64 equally spaced points of the bracket, including both ends, in one vectorized call. If the ends have the same sign, the solver stops with an error that reports the time and state, since that indicates a physical or indexing problem rather than a solver problem. Otherwise the grid interval with the first sign change becomes the new bracket. Two such rounds narrow it by a factor of about `63^2 = 4e3`, about as much as 12 bisection steps. Three Newton steps with the analytic derivative `dr/d(kv_w)` follow, each kept inside the final interval. The bracketed solve has no adaptive termination, so it always costs the same: 2 vectorized evaluations of 64 points, 3 Newton evaluations and 1 final evaluation. Against 3 rounds of 32 points, this layout gives the same roots (round-off, 360 trajectory states) at 1.0 to 1.4 times less cost per call.

Over 360 states from mild and violent collapses, with `kv_w` between `0.0015` and `0.9992`, the root agrees with a tight solve of the `imr.physics.wall_mass_balance` balance to `4e-16` in `theta_w`. It needs a fixed number of residual evaluations, where the earlier closure took 20 scalar bisection steps with central-difference Newton.

## Reported spectral stress comparisons

The former spectral solver page reported that Gauss–Legendre weights with `Nv + 64` points agree with adaptive integration in `theta` to `1e-12` for `Lv` 0.3–100 and `Nv` up to 150. It reported setup time of 0.05 s instead of 0.6 s for `Nv = 150` when replacing `Nv` adaptive integrations in `x`.

For mild oscillations, the same page reported exponential convergence with `Nv` and agreement of the spectral stress field with the two-state form to about `1e-11` from `Nv = 16`: Keller–Miksis, `Req = 0.7` and 0.9, `Lv = 0.5` and 3. Strong collapses need more modes to resolve the wall stress layer.

## Reported Kirchhoff transform correction

The former nodal page reported that replacing an inverse assuming `beta = 1 - alpha` with the actual conductivity coefficients changed `R(t)` by up to `2.4e-4` over its regression cases. For the default gas coefficients, `alpha_g + beta_g = 1.0024`; this history is retained as a reported result, not a rerun.

## Reported local Jacobian checks

The former nodal method notes reported agreement of Zener-field rows with central differences to `1e-9`. The spectral method notes reported `5e-10` for stress 1, 3, 4 and 5 in one set of checks. These narrower reports differ in scope from the broader Jacobian comparisons above and are not interchangeable acceptance thresholds.
