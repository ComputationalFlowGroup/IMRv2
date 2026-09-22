# `f_imr_fd.m` Documentation

## Overview
`f_imr_fd` is a second-order accurate finite difference solver of the PDEs for thermal transport (inside the bubble and in the surrounding medium) and vapor mass transfer, coupled to the radial (Rayleigh–Plesset family) bubble dynamics and viscoelastic stress models.
It uses a non-dimensionalized formulation and integrates the semi-discrete system in time with MATLAB's ODE solvers through `f_odesolve`.

---

## Major Sections of the Code

### 1. Initialization
- Calls `f_call_params(varargin{:})` to set up all problem parameters.
- Unpacks solver settings, physical properties, initial conditions, and waveform settings.

### 2. Mesh and Finite Difference Matrix Setup
- Builds a uniform grid inside the bubble and a stretched grid (length scale `Lt`) in the medium.
- Builds sparse first- and second-derivative matrices (`f_finite_diff_mat`, via `spdiags`): second-order central differences in the interior, a second-order one-sided closure at the bubble wall, and the symmetry condition at the bubble center.

### 3. Time Integration
- Assembles the state vector `[R; Rdot; P; theta; Tm; kv; Z]` and integrates it with `f_odesolve` (`ode15s`, `ode23tb`, or `ode45`).
- The bubble-wall temperature is set by the heat flux balance: in closed form for thermal-only runs and by a warm-started Newton iteration (with `fzero` fallback) when mass transfer is active.
- The stiff integrators receive an analytic Jacobian (`bubtherm = 1`, `masstrans = 0`) or the sparsity pattern from `f_jpattern` otherwise; see the `jacobian` option.

### 4. Output Assembly
- Returns `[t, R, Rdot, P, T, Tm, kv]` through `varargout`, in dimensional units if `dimout = 1`.

---

## Inputs

Inputs are passed as name/value pairs, routed through `f_call_params`, and then unpacked locally.

| Variable          | Description |
|-------------------|-------------|
| `radial`           | Radial dynamics model (1–6) |
| `bubtherm`         | Enable bubble thermal model (0 or 1) |
| `medtherm`         | Enable medium thermal model (0 or 1) |
| `stress`           | Stress model selection (0–5) |
| `eps3`             | Include second-order strain term (0/1) |
| `masstrans`        | Enable mass transfer (0/1, requires `bubtherm = 1`) |
| `method`           | Time integrator (15 = `ode15s`, 23 = `ode23tb`, 45 = `ode45`) |
| `reltol`, `abstol` | Integrator tolerances (default `1e-8`) |
| `jacobian`         | Jacobian for the stiff integrators (0 = dense numerical, 1 = sparsity pattern, 2 = analytic where implemented, default) |
| `spectral`         | Enable spectral stress representation (should be 0 here) |
| `divisions`, `Nv`, `Nt`, `Mt` | Maximum step control and grid sizes |
| `Lt`               | Grid stretching length in the medium |
| `Rzero`, `Uzero`   | Initial radius and velocity (non-dimensional) |
| `Pb_star`, `Pv_star`, `P8` | Non-dimensional pressures |
| `T8`               | Far-field temperature |
| `Req`              | Equilibrium radius |
| `alphax`           | Strain-stiffening parameter |
| `Szero`            | Initial stress state |
| `tspan`            | Time vector (non-dimensional) |
| `dimensionalout`   | Output in physical units (0/1) |
| `progdisplay`      | Show progress during run (0/1); prints on every right-hand-side call |
| `Cstar`, `GAMa`, `kappa`, `nstate`, `hugoniot_s` | Acoustic properties |
| `om`, `ee`, `tw`, `dt`, `mn`, `wave_type` | Acoustic waveform parameters |
| `wave_poly`, `wave_dpoly` | Waveform spline fits (for `wave_type < 0`) |
| `We`, `Re8`, `DRe` | Dimensionless fluid dynamics numbers |
| `v_a`, `v_nc`, `Ca`, `LAM`, `De`, `JdotA` | Viscoelastic model parameters |
| `nu_model`, `v_lambda_star`, `zeNO`, `iDRe` | Viscosity model settings (`nu_model` 0–2) |

---

## Outputs

| Output | Description |
|--------|-------------|
| `t`    | Time vector |
| `R`    | Bubble radius history |
| `Rdot` | Wall velocity history |
| `P`    | Bubble pressure history |
| `T`    | Bubble temperature field (if `bubtherm = 1`) |
| `Tm`   | Medium temperature field (if `medtherm = 1`, otherwise the far-field value) |
| `kv`   | Vapor concentration field (if `masstrans = 1`) |

---

## Notes

- The spatial discretization is second-order accurate.
- Medium heating for `stress >= 3` currently uses the Kelvin–Voigt closed form for the dissipation; `f_call_params` warns (`IMR:dissipation`) when this applies.
- Mesh resolution must be tuned carefully for high Weissenberg numbers or sharp waveforms.
