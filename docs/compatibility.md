# Model compatibility

These grids show which features can run together in each solver: `f_imr_nodal`, with either differentiation scheme, and `f_imr_spectral`. Every feature appears on both axes, and only the upper triangle is filled in.

Legend:
- **✓** the two features run together.
- **✗** they cannot run together. `f_checkers` rejects the run with an `INPUT ERROR`.
- **–** the diagonal, where a feature is paired with itself.
- **·** alternative settings of the same option, such as two radial models, which cannot be selected at once.

A ✗ on the diagonal means the solver does not have that feature at all, so its whole row is ✗.

The radial models fall into three groups with the same compatibility:

| Radial | Model | Liquid | Group |
|---|---|---|---|
| 1 | Rayleigh–Plesset | Incompressible | 1–2 |
| 2 | Keller–Miksis, pressure form | Constant sound speed | 1–2 |
| 3 | Keller–Miksis, enthalpy form | Tait equation of state | 3–4 |
| 4 | Gilmore | Tait equation of state | 3–4 |
| 5 | Keller–Miksis, enthalpy form | Mie–Grüneisen equation of state | 5–6 |
| 6 | Gilmore | Mie–Grüneisen equation of state | 5–6 |

The enthalpy forms (radial 3–6) describe the liquid with an equation of state in pressure only, at the far-field temperature. Medium heat transfer would change the liquid temperature, so it is not available with them, and neither is collapse, which requires it. Bubble heat transfer and mass transfer remain available with radial 3–4: without medium heat transfer the wall stays at the far-field temperature, and the wall vapor fraction is the saturated value at that temperature and the bubble pressure. The liquid then acts as an infinite heat reservoir, which is consistent with the equation of state but neglects the evaporative cooling and condensation heating of the wall. Radial 5–6 remain limited to the polytropic bubble.

### `f_imr_nodal`, either scheme (finite difference or Chebyshev)

| | Feature | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---:|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| 1 | Radial 1–2 | – | · | · | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ |
| 2 | Radial 3–4 |  | – | · | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✗ | ✓ |
| 3 | Radial 5–6 |  |  | – | ✓ | ✗ | ✗ | ✗ | ✓ | ✓ | ✓ | ✗ | ✓ | ✗ | ✓ |
| 4 | Vapor pressure |  |  |  | – | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ |
| 5 | Bubble heat transfer |  |  |  |  | – | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ |
| 6 | Medium heat transfer |  |  |  |  |  | – | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ |
| 7 | Mass transfer |  |  |  |  |  |  | – | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ |
| 8 | Kelvin–Voigt stress (1, 2) |  |  |  |  |  |  |  | – | · | · | · | ✓ | ✓ | ✓ |
| 9 | Maxwell / Jeffreys / Zener (3, 4) |  |  |  |  |  |  |  |  | – | · | · | ✓ | ✓ | ✓ |
| 10 | UCM / Oldroyd-B, two states (5) |  |  |  |  |  |  |  |  |  | – | · | ✓ | ✓ | ✓ |
| 11 | UCM / Oldroyd-B, stress field (5, `spectral = 1`) |  |  |  |  |  |  |  |  |  |  | ✗ | ✗ | ✗ | ✗ |
| 12 | Non-Newtonian viscosity |  |  |  |  |  |  |  |  |  |  |  | – | ✓ | ✓ |
| 13 | Collapse |  |  |  |  |  |  |  |  |  |  |  |  | – | ✓ |
| 14 | Pressure waveform |  |  |  |  |  |  |  |  |  |  |  |  |  | – |

### `f_imr_spectral`

| | Feature | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---:|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| 1 | Radial 1–2 | – | · | · | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ |
| 2 | Radial 3–4 |  | – | · | ✓ | ✓ | ✗ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ |
| 3 | Radial 5–6 |  |  | – | ✓ | ✗ | ✗ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ |
| 4 | Vapor pressure |  |  |  | – | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ |
| 5 | Bubble heat transfer |  |  |  |  | – | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ |
| 6 | Medium heat transfer |  |  |  |  |  | – | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ |
| 7 | Mass transfer |  |  |  |  |  |  | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ |
| 8 | Kelvin–Voigt stress (1, 2) |  |  |  |  |  |  |  | – | · | · | · | ✓ | ✗ | ✓ |
| 9 | Maxwell / Jeffreys / Zener (3, 4) |  |  |  |  |  |  |  |  | – | · | · | ✓ | ✗ | ✓ |
| 10 | UCM / Oldroyd-B, two states (5) |  |  |  |  |  |  |  |  |  | – | · | ✓ | ✗ | ✓ |
| 11 | UCM / Oldroyd-B, stress field (5, `spectral = 1`) |  |  |  |  |  |  |  |  |  |  | – | ✓ | ✗ | ✓ |
| 12 | Non-Newtonian viscosity |  |  |  |  |  |  |  |  |  |  |  | – | ✗ | ✓ |
| 13 | Collapse |  |  |  |  |  |  |  |  |  |  |  |  | ✗ | ✗ |
| 14 | Pressure waveform |  |  |  |  |  |  |  |  |  |  |  |  |  | – |

The finite difference and Chebyshev schemes of `f_imr_nodal` were run separately and gave the same grid. One run is worth knowing about: Rayleigh–Plesset (radial 1) with collapse on 12 Chebyshev points compressed the bubble to `R = 0.005` and `P = 3.6e6`, above the ceiling `C1_pv` of the saturation curve. The wall closure then stops with its sign-change error, since no wall vapor fraction satisfies the saturation condition there. Rayleigh–Plesset has no acoustic damping, so strong collapses are best run with radial 2.

## Requirements

Some features switch on others. A ✓ means the pair runs with these requirements met:
- **Medium heat transfer** requires bubble heat transfer.
- **Mass transfer** requires vapor pressure and bubble heat transfer.
- **Collapse** requires bubble and medium heat transfer, vapor pressure and mass transfer, so only radial 1–2.
- **Medium heat transfer** is not available with the enthalpy forms, radial 3–6.

Other rules:
- **Stress 3, 4 and 5** require `De > 0`.
- **Stress 5** sets `Ca = Inf`, so it has no elastic contribution.
- **Medium heating with stress 3 to 5** uses the Kelvin–Voigt closed form, which is inconsistent for history-dependent stress, and `f_checkers` warns. The stress field option (feature 10) computes the heating from the stress field itself.

## Why the spectral solver cannot collapse

Collapse only sets the initial state in `f_call_params`: the equilibrium radius, a zero wall velocity and, for stress 3 and 4, the initial stress. No solver refers to it afterwards, and the initial-state calculation uses the vapor and gas properties but nothing from mass transfer.

`f_imr_spectral` is blocked only because `f_checkers` requires mass transfer for collapse, and it has no mass transfer. That requirement dates from the original code. With it relaxed to bubble and medium heat transfer plus vapor pressure, `f_imr_nodal` and `f_imr_spectral` ran a collapse with stress 0, 1, 3 and 5.

## How this was checked

Every cell comes from running the solvers with GNU Octave, `f_imr_nodal` in both schemes. Each run pairs the two features on a Keller–Miksis base case, radial 2, and each ran in its own Octave process. Radial 1–6 were each paired with every non-radial feature in all three solver modes, and each gave the same result as the rest of its group, apart from the Rayleigh–Plesset collapse run described above. A ✓ means the run finished with a finite `R(t)`.

A separate sweep ran both solvers over every radial model, thermal configuration and stress model: 504 runs, plus 36 more for the remaining thermal combinations. Every non-Newtonian viscosity model, driving waveform, integrator and Jacobian option was also run on a thermal case in both solvers.

Octave needs two stand-ins for MATLAB behavior:
- **No `ode23tb`.** `ode15s` with a dense Jacobian stood in for it.
- **A different `fzero`.** Octave's `fzero` throws an error where MATLAB's returns `NaN`, which breaks the collapse initial state. A stand-in with MATLAB's interval search was used.

One Octave-only problem remains: `f_imr_nodal` with `method = 15` crashes Octave's `ode15s` on the sparse Jacobian. MATLAB is not affected.
