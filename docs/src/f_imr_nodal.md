# `f_imr_nodal.m` Documentation

## Overview
`f_imr_nodal` solves the PDEs for thermal transport (inside the bubble and in the surrounding medium) and vapor mass transfer, coupled to the radial (Rayleigh–Plesset family) bubble dynamics and viscoelastic stress models. The fields are carried as node values. The `scheme` input selects the spatial discretization, and the differentiation matrices from [`f_diff_mat`](f_diff_mat.md) are the only difference between the two:

| `scheme` | Discretization | Available models |
|---|---|---|
| 0 (default) | Second-order finite differences, uniform bubble grid | All, including mass transfer and collapse |
| 1 | Chebyshev collocation | All, including mass transfer and collapse |

It uses a non-dimensionalized formulation and integrates the semi-discrete system in time with MATLAB's ODE solvers through `f_odesolve`. The solver was `f_imr_fd` before the Chebyshev scheme was added; `scheme = 0` reproduces `f_imr_fd` bit for bit.

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
- Assembles the state vector `[R; Rdot; P; theta; Tm; kv; Z]` and integrates it with `f_odesolve` (`ode15s`, `ode23tb`, or `ode45`).
- The bubble-wall temperature is set by the heat flux balance: in closed form for thermal-only runs (`f_wall_theta_exact`, `fzero` only as a fallback) and, with mass transfer, by a warm-started Newton iteration on the transcendental balance of `f_wall_mass_bc` (with `fzero` as a fallback).
- The stiff integrators receive the analytic Jacobian of `f_jacobian_nodal` for `bubtherm = 1` (with or without medium heat and mass transfer), or the sparsity pattern from `f_jpattern` otherwise (dense for Chebyshev); see the `jacobian` option.

### 4. Output Assembly
- Returns `[t, R, Rdot, P, T, Tm, kv]` through `varargout`. All outputs are dimensionless.

---

## Inputs

Inputs are dimensionless name/value pairs (see `default_case`), routed through `f_call_params` and checked by `f_checkers`, then unpacked locally.

| Variable          | Description |
|-------------------|-------------|
| `radial`           | Radial dynamics model (1–6; 5 and 6 are Mie–Grüneisen and polytropic only) |
| `bubtherm`         | Enable bubble thermal model (0 or 1) |
| `medtherm`         | Enable medium thermal model (0 or 1) |
| `stress`           | Stress model selection (0–5) |
| `eps3`             | Second-order strain term, in [0, 0.5] |
| `masstrans`        | Enable mass transfer (0/1, requires `bubtherm = 1` and `vapor = 1`) |
| `scheme`           | Differentiation matrices: 0 finite difference (default), 1 Chebyshev |
| `method`           | Time integrator (15 = `ode15s`, 23 = `ode23tb`, 45 = `ode45`) |
| `reltol`, `abstol` | Integrator tolerances (default `1e-8`) |
| `jacobian`         | Jacobian for the stiff integrators: 0 dense numerical, 1 sparsity pattern, 2 analytic (default). The analytic Jacobian needs `bubtherm = 1`; otherwise 2 falls back to the sparsity pattern. |
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

---

## Notes

- The finite difference discretization is second-order accurate; the Chebyshev discretization converges spectrally.
- Medium heating for `stress >= 3` currently uses the Kelvin–Voigt closed form for the dissipation; `f_checkers` warns (`IMR:dissipation`) when this applies.
- Collapse (`collapse = 1`) requires bubble and medium heat transfer, vapor and mass transfer.
- See [`compatibility`](../compatibility.md) for which features run together.
- Mesh resolution must be tuned carefully for high Weissenberg numbers or sharp waveforms.

## Verification of the Chebyshev scheme

Measured in GNU Octave (`ode15s` standing in for `ode23tb`), on thermal cases with `Req = 0.5` and `t` in [0, 3]:
- `scheme = 0` gives bit-identical `R(t)` to `f_imr_fd` over 98 cases: radial 1–4, stress 0–5, polytropic, bubble, bubble and medium, and full mass-transfer models, plus a non-Newtonian and a Gaussian-waveform case.
- `scheme = 1` and `f_imr_spectral` converge toward each other with resolution over 32 cases (radial 1–4, stress 0/1/3/5, with and without medium heat transfer). The largest relative deviation in `R(t)` over the cases is `2.8e-3`, `4.8e-4` and `2.6e-4` at `Nt = Mt = 12, 20, 30`; the median is `1.4e-4`, `1.0e-5` and `3.4e-6`.
- Against Richardson-extrapolated finite differences (`Nt = Mt = 200, 400`) for Keller–Miksis with neo-Hookean Kelvin–Voigt stress, `scheme = 1` at `Nt = Mt = 25` is within `1e-4` with bubble heat transfer only, where finite differences at 200 points are `2.4e-3` away. With medium heat transfer it reaches `2e-4` at 40 points, the same as `f_imr_spectral`.
- With mass transfer (Keller–Miksis, neo-Hookean Kelvin–Voigt, bubble and medium heat transfer), `scheme = 1` converges with resolution: against its own solution at 80 points the largest deviation is `1.2e-2`, `3.6e-3`, `2.6e-4` and `5.5e-5` at 20, 30, 40 and 60 points for `Req = 0.5`, and `4.9e-2`, `2.0e-2`, `2.5e-3` and `1.1e-4` for `Req = 0.4`. At 40 points it is within `1e-3` of Richardson-extrapolated finite differences (`Nt = Mt = 200, 400`), the uncertainty of that reference.
- The analytic Jacobian matches a central-difference Jacobian to a scaled error of `3e-9` for `scheme = 0` and `4e-10` for `scheme = 1`, with and without mass transfer.
