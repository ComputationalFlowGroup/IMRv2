# `f_imr_spectral.m` Documentation

## Overview
`f_imr_spectral` is a Chebyshev spectral collocation solver for the bubble dynamics (Rayleigh–Plesset family) with optional heat transfer in the bubble and the surrounding medium and viscoelastic stress models, including a spectral representation of the UCM / Oldroyd-B stress field. All inputs and outputs are dimensionless (see `default_case`). Vapor mass transfer is not available in this solver (`f_checkers` rejects `masstrans = 1`); use `f_imr_nodal`. Because `f_checkers` requires mass transfer for collapse, collapse is not available either. See [`compatibility`](../compatibility.md).

---

## Major Sections of the Code

### 1. Initialization
- Calls `f_call_params('spectral', varargin{:})`, which reads and checks the inputs.

### 2. Collocation grids and operators
- Even Chebyshev expansion of the bubble temperature and a Chebyshev expansion of the medium temperature on the stretched coordinate (length `Lt`).
- The combined operators from Chebyshev coefficients to collocation values of the first and second derivatives (`gDA`, `gDDA`, `mDA`, `mDDA`) are formed once.
- The thermal boundary condition system has a single time-dependent row; it is factored once and updated with a Sherman–Morrison rank-one correction on every call.
- The constant factors of the medium temperature equation and of the spectral stress transport are formed once, with the same arithmetic as in the right-hand side, so results are bit-identical. With `wave_type = 0` the far-field pressure is the constant `ee` and `f_pinfinity` is not called.

### 3. Stress models
- `stress = 0–4`: closed-form or single-state models through `f_stress`.
- `stress = 3, 4` with `medtherm = 1`: `Z1`, then the `Nv` Chebyshev coefficients of the local polymer stress field ([`f_zener_field`](f_zener_field.md)), in the same basis as the stress 5 field, for the medium heating. After a run the solver prints the largest mismatch of the field's stress integral against `Z1/R^3`.
- `stress = 5`, `spectral = 0`: UCM / Oldroyd-B through the two stress states `Z1`, `Z2`. With `medtherm = 1`, `f_checkers` switches to the stress field (`spectral = 1`) with a warning, since the heating needs the local stress.
- `stress = 5`, `spectral = 1`: UCM / Oldroyd-B with the polymer stresses `tau_rr` and `tau_tt` expanded in `T_n(z) - 1`, `n = 1..Nv`, on `z = 1 - 2/(1 + (r/R - 1)/Lv)` (collocation points `z_j = cos(pi*j/Nv)`). The stress integral is `S = 2*cdd'*(c - d) - 4*LAM/Re8*Rdot/R - 6*intfnu*iDRe`. The weights `cdd(n) = 2*Lv*int_{-1}^{1} (T_n(x) - 1)/((Lv*(1 + x) + 1 - x)*(1 - x)) dx` come from one Gauss–Legendre rule with `Nv + 64` points in `theta = acos(x)`, where the integrands are smooth; it matches an adaptive integral in `theta` to `1e-12` for `Lv` 0.3–100 and `Nv` up to 150, and replaces `Nv` adaptive integrals in `x` (0.6 s against 0.05 s of setup at `Nv = 150`). The stress field also provides the medium heating `tau:grad(u)`. For a mild oscillation the representation converges exponentially in `Nv` and agrees with the two-state form to about `1e-11` from `Nv = 16` (Keller–Miksis, `Req = 0.7` and `0.9`, `Lv = 0.5` and `3`, measured in GNU Octave); strong collapses produce a steep stress layer at the wall that needs many modes, where the two-state form is the cheaper choice for `R(t)`.

### 4. Time integration
- The bubble temperature is carried as the Kirchhoff variable `theta = alpha_g*(T^2 - 1)/2 + beta_g*(T - 1)` of the gas conductivity `K = alpha_g*T + beta_g`, recovered as `T = (sqrt((alpha_g + beta_g)^2 + 2*alpha_g*theta) - beta_g)/alpha_g`. With medium heat transfer the wall condition scales by the wall conductivity `alpha_g*T + beta_g`. Both coefficients are used as given; before, the inverse assumed `beta_g = 1 - alpha_g`.
- `f_odesolve` with the analytic Jacobian of `f_jacobian_spectral` when `jacobian = 2` (default).

### 5. Output
- Returns `[t, R, Rdot, P, T, Tm]` (dimensionless); the seventh output is empty. The optional 8th output is the struct of the run's inputs (see [`f_imr_inputs`](f_imr_inputs.md)).

---

## Inputs

Dimensionless name/value pairs, see `default_case`, or a struct from [`f_imr_inputs`](f_imr_inputs.md). Solver-specific notes:

| Variable | Description |
|----------|-------------|
| `radial` | Radial dynamics model (1–6; 5 and 6 are Mie–Grüneisen and polytropic only) |
| `stress` | Stress model (0–5) |
| `spectral` | Spectral stress representation for `stress = 5` (0/1) |
| `nv`, `lv` | Number of stress modes and stretching length of the stress coordinate |
| `nt`, `mt`, `lt` | Thermal expansion orders and medium stretching length |
| `jacobian` | 2 analytic Jacobian (default); 0 and 1 give a dense numerical Jacobian; 3 (batched finite differences) is nodal only and raises an error here |
