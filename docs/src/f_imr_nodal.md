# `f_imr_nodal.m` Documentation

## Overview
`f_imr_nodal` solves the PDEs for thermal transport (inside the bubble and in the surrounding medium) and vapor mass transfer, coupled to the radial (Rayleigh–Plesset family) bubble dynamics and viscoelastic stress models. The fields are carried as node values. The `scheme` input selects the spatial discretization, and the differentiation matrices from [`f_diff_mat`](f_diff_mat.md) are the only difference between the two:

| `scheme` | Discretization | Available models |
|---|---|---|
| 0 (default) | Second-order finite differences, uniform bubble grid | All, including mass transfer and collapse |
| 1 | Chebyshev collocation | All, including mass transfer and collapse |

It uses a non-dimensionalized formulation and integrates the semi-discrete system in time with MATLAB's ODE solvers through `f_odesolve`. The solver was `f_imr_fd` before the Chebyshev scheme was added. Without mass transfer, `scheme = 0` reproduces `f_imr_fd` bit for bit; with mass transfer, the bracketed wall closure moves `R(t)` by integrator-tolerance amounts.

---

## Major Sections of the Code

### 1. Initialization
- Calls `f_call_params('nodal', varargin{:})`, which reads the inputs and checks them with `f_checkers`.
- Unpacks solver settings, physical properties, initial conditions, and waveform settings.

### 2. Grids and Differentiation Matrices
- `f_diff_mat` returns the node positions `y = r/R` in the bubble and `xi` in the medium, with the first- and second-derivative matrices for the chosen `scheme`. The medium coordinate is mapped to the radius with the stretching length `Lt`.
- Finite differences: uniform grids, second-order central differences in the interior, a second-order one-sided row at the bubble wall, and the symmetry condition at the bubble center.
- Chebyshev: Chebyshev points, with the bubble operators built for even functions so that the symmetry condition at the center holds exactly.
- The bubble-wall heat flux balance, with or without mass transfer, uses the wall rows of the first-derivative matrices: three-point one-sided stencils for finite differences, and every node of each field for Chebyshev.

### 3. Time Integration
- Assembles the state vector `[R; Rdot; P; theta; Tm; kv; Z]` and integrates it with `f_odesolve` (`ode15s`, `ode23tb` or CVODES).
- The bubble temperature is carried as the Kirchhoff variable `theta = alpha*(T^2 - 1)/2 + beta*(T - 1)` of the conductivity `K = alpha*T + beta`, and recovered as `T = (sqrt((alpha + beta)^2 + 2*alpha*theta) - beta)/alpha`, where the square root is `K`. `alpha` and `beta` are the mixture values of the gas and vapor coefficients: of `kv0` without mass transfer, of the local `kv` with it. Both are used as given, so the transform, its inverse and the conductivity in the diffusion term agree. Before, the inverse assumed `beta = 1 - alpha` while the diffusion term used the actual `beta` (`alpha_g + beta_g = 1.0024` in the default case), which changed `R(t)` by up to `2.4e-4` over the regression cases.
- The bubble-wall temperature is set by the heat flux balance: in closed form for thermal-only runs (`f_wall_theta_exact`, `fzero` only as a fallback) and, with mass transfer, by a solve for the wall vapor mass fraction: Newton from the previous root, with a fixed-cost bracketed solve as the fallback (see [`f_wall_bc`](f_wall_bc.md)).
- The stiff integrators receive the analytic Jacobian of `f_jacobian_nodal` for `bubtherm = 1` (with or without medium heat and mass transfer), or the sparsity pattern from `f_jpattern` otherwise (dense for Chebyshev); see the `jacobian` option.
- With `jacobian = 3` they receive a finite-difference Jacobian from one batched right-hand-side call (see [Batched right-hand side](#batched-right-hand-side)). It works for every model.

### 4. Output Assembly
- Returns `[t, R, Rdot, P, T, Tm, kv]` through `varargout`, and optionally the struct of the run's inputs as an 8th output. All outputs are dimensionless.
- The wall entries of the state are not integrated (their time derivatives are zero), so the wall values of `T`, `Tm` and `kv` are recovered at every output time from the same wall closure as the right-hand side. `T(:,end)` and `Tm(:,1)` are both the wall temperature.

---

## Inputs

Inputs are dimensionless name/value pairs (see `default_case`) or a struct from [`f_imr_inputs`](f_imr_inputs.md), routed through `f_call_params` and checked by `f_checkers`, then unpacked locally.

| Variable          | Description |
|-------------------|-------------|
| `radial`           | Radial dynamics model (1–6; 5 and 6 are Mie–Grüneisen and polytropic only) |
| `bubtherm`         | Enable bubble thermal model (0 or 1) |
| `medtherm`         | Enable medium thermal model (0 or 1) |
| `stress`           | Stress model selection (0–5) |
| `eps3`             | Second-order strain term, in [0, 0.5]; currently unused by the stress models |
| `masstrans`        | Enable mass transfer (0/1, requires `bubtherm = 1` and `vapor = 1`) |
| `scheme`           | Differentiation matrices: 0 finite difference (default), 1 Chebyshev |
| `method`           | Stiff time integrator: 23 or `'ode23tb'` (default), 15 or `'ode15s'`, 3 or `'cvodes'` (CVODES through the `ode` object, MATLAB R2024a or later) |
| `reltol`, `abstol` | Integrator tolerances (default `1e-8`) |
| `jacobian`         | Jacobian for the stiff integrators: 0 dense numerical, 1 sparsity pattern, 2 analytic (default), 3 batched finite differences. The analytic Jacobian needs `bubtherm = 1`; otherwise 2 falls back to the sparsity pattern. |
| `spectral`         | Spectral stress representation; `f_checkers` rejects `spectral = 1` with stress 5 in this solver |
| `divisions`, `Nv`, `Nt`, `Mt` | Minimum number of time steps, and grid sizes |
| `Lt`               | Grid stretching length in the medium |
| `Rzero`, `Uzero`   | Initial radius and velocity (non-dimensional) |
| `Pb_star`, `Pv_star` | Non-dimensional pressures |
| `Req`              | Equilibrium radius |
| `alphax`           | Strain-stiffening parameter |
| `Szero`            | Initial stress state |
| `tspan`            | Time vector (non-dimensional) |
| `progdisplay`      | Show progress during run (0/1); prints on every right-hand-side call |
| `Cstar`, `GAMa`, `kappa`, `nstate`, `hugoniot_s` | Acoustic properties |
| `om`, `ee`, `tw`, `dt`, `mn`, `wave_type` | Acoustic waveform parameters |
| `wave_poly`, `wave_dpoly` | Waveform spline fits (for `wave_type < 0`) |
| `We`, `Re8`, `DRe` | Dimensionless fluid dynamics numbers |
| `v_a`, `v_nc`, `Ca`, `LAM`, `De`, `JdotA` | Viscoelastic model parameters |
| `nu_model`, `v_lambda_star`, `iDRe` | Viscosity model settings (`nu_model` 0–7) |

---

## Outputs

| Output | Description |
|--------|-------------|
| `t`    | Time vector |
| `R`    | Bubble radius history |
| `Rdot` | Wall velocity history |
| `P`    | Bubble pressure history |
| `T`    | Bubble temperature field (if `bubtherm = 1`, otherwise empty) |
| `Tm`   | Medium temperature field (if `medtherm = 1`, otherwise the far-field value) |
| `kv`   | Vapor concentration field (if `masstrans = 1`, otherwise empty) |
| `inputs` | Optional 8th output: the struct of every input of the run, as from [`f_imr_inputs`](f_imr_inputs.md) with the given changes |

---

## Notes

- The finite difference discretization is second-order accurate; the Chebyshev discretization converges spectrally.
- **Medium heating with stress 3 and 4:** the solver carries the local polymer stress `sp` at the medium nodes after `Z1` ([`f_zener_field`](f_zener_field.md)). The state is `[R; Rdot; P; theta; Tm; kv; Z1; sp]`.
- **Stress 5 with medium heat transfer:** not available in this solver. Its two-state form has no local stress for the heating, so use `f_imr_spectral`.
- Collapse (`collapse = 1`) requires bubble and medium heat transfer, vapor and mass transfer, and with stress 5 also `stress0`.
- See [`compatibility`](../compatibility.md) for which features run together.
- Mesh resolution must be tuned carefully for high Weissenberg numbers or sharp waveforms.
- Cost per right-hand-side call: the constant factors of the medium equation and of the bubble temperature map are computed once, with the same arithmetic as before, so results are bit-identical; with `wave_type = 0` the far-field pressure is the constant `ee` and `f_pinfinity` is not called. In GNU Octave this made a right-hand-side call 11–31% cheaper and most runs 8–17% faster.
- A nested function shares every variable name that the solver function also uses, so the nested functions use their own loop variable names (the output loop is `it_out`).

## Batched right-hand side

The right-hand side `SVBDODE` takes a matrix `X` with one state per column and returns one column of derivatives per state. Every operation is elementwise across the columns, in the same order as for a single state. The helpers it calls work the same way: `f_stress`, `f_radial_eq`, `f_viscosity`, `f_stress_dissipation`, and the wall closures. A single column therefore gives the same result as before bit for bit. The integrators call it with one column.

`jacobian = 3` uses the batch. It evaluates the state and its `n` perturbations in one call, `F = SVBDODE(t, [X, X + diag(del)])`, and forms `J = (F(:,2:end) - F(:,1))./del.'`. The step is `sqrt(eps)*max(|X|, 1)`, rounded so that the perturbed state is exactly `X + del`. Entries outside the `f_jpattern` pattern are set to zero. `J` is returned dense: it is computed dense, dense LU is fast at these sizes, and a sparse `J` whose number of nonzeros changes between calls breaks Octave's `ode15s`.

Why it pays: an interpreted operation costs about the same for one number as for a few hundred. A batch of `n + 1` states costs about 3 right-hand-side calls, not `n + 1`. The analytic Jacobian is built from 55 sparse dual-number operations and costs 25–40 right-hand-side calls.

Measured in GNU Octave at the default grids (`Nt = Mt = 25`):

| Case | States | Right-hand side | Analytic `J` | Batched FD `J` |
|------|--------|-----------------|--------------|----------------|
| Bubble heat transfer | 28 | 0.12 ms | 4.6 ms | 0.30 ms |
| Bubble and medium heat transfer | 53 | 0.25 ms | 6.7 ms | 0.65 ms |
| Kelvin–Voigt, Zener stress | 54 | 0.27 ms | 6.9 ms | 0.68 ms |
| Carreau viscosity | 53 | 0.36 ms | 8.5 ms | 1.2 ms |
| Mass transfer | 78 | 0.44 ms | 15.0 ms | 1.4 ms |
| Mass transfer, Chebyshev | 78 | 0.53 ms | 15.8 ms | 1.5 ms |
| Medium heat transfer, `Nt = Mt = 40` | 83 | 0.24 ms | 6.8 ms | 1.0 ms |

Checks:

- Every column of a batch matches a single-column call to `4e-13` (relative), over 372 model combinations: radial 1–6, stress 0–5, four viscosity models, the thermal and mass-transfer configurations, and both schemes.
- The batched Jacobian matches the analytic one to a scaled error below `1e-5` in every thermal case.
- Over the 74 regression cases, `jacobian = 3` reproduces the `jacobian = 2` solution to within the integrator tolerance (`R(t)` differs by at most `9e-7`, mostly `1e-8` or less). The total run time went from 43 s to 28 s in Octave, where `ode15s` is itself compiled. MATLAB's `ode15s` and `ode23tb` are interpreted, so the gain there is expected to differ.
- `benchmarks/integrators/s_compare_integrators.m` compares `ode23tb`, `ode15s` and CVODES with both Jacobians against a tight-tolerance reference, in MATLAB.

## Verification of the Chebyshev scheme

Measured in GNU Octave (`ode15s` standing in for `ode23tb`), on thermal cases with `Req = 0.5` and `t` in [0, 3]:

- `scheme = 0` gave bit-identical `R(t)` to `f_imr_fd` over 98 cases: radial 1–4, stress 0–5, polytropic, bubble, bubble and medium, and full mass-transfer models, plus a non-Newtonian and a Gaussian-waveform case. This held before the bracketed mass-transfer wall closure; with it, the cases without mass transfer stay bit-identical and the mass-transfer cases differ at the integrator tolerance.
- `scheme = 1` and `f_imr_spectral` converge toward each other with resolution over 32 cases (radial 1–4, stress 0/1/3/5, with and without medium heat transfer). The largest relative deviation in `R(t)` over the cases is `2.8e-3`, `4.8e-4` and `2.6e-4` at `Nt = Mt = 12, 20, 30`; the median is `1.4e-4`, `1.0e-5` and `3.4e-6`.
- Against Richardson-extrapolated finite differences (`Nt = Mt = 200, 400`) for Keller–Miksis with neo-Hookean Kelvin–Voigt stress, `scheme = 1` at `Nt = Mt = 25` is within `1e-4` with bubble heat transfer only, where finite differences at 200 points are `2.4e-3` away. With medium heat transfer it reaches `2e-4` at 40 points, the same as `f_imr_spectral`.
- With mass transfer (Keller–Miksis, neo-Hookean Kelvin–Voigt, bubble and medium heat transfer), `scheme = 1` converges with resolution: against its own solution at 80 points the largest deviation is `1.2e-2`, `3.6e-3`, `2.6e-4` and `5.5e-5` at 20, 30, 40 and 60 points for `Req = 0.5`, and `4.9e-2`, `2.0e-2`, `2.5e-3` and `1.1e-4` for `Req = 0.4`. At 40 points it is within `1e-3` of Richardson-extrapolated finite differences (`Nt = Mt = 200, 400`), the uncertainty of that reference.
- The analytic Jacobian matches a central-difference Jacobian to a scaled error of `3e-9` for `scheme = 0` and `4e-10` for `scheme = 1`, with and without mass transfer.
