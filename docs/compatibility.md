# Model compatibility

These grids show which features can run together in each solver. Every feature appears on both axes, and only the upper triangle is filled in.

Legend:
- **✓** the two features run together.
- **✗** they cannot run together. `f_checkers` rejects the run with an `INPUT ERROR`.
- **–** the diagonal, where a feature is paired with itself.
- **·** alternative settings of the same option, such as two radial models, which cannot be selected at once.

A ✗ on the diagonal means the solver does not have that feature at all, so its whole row is ✗.

### `f_imr_fd`

| | Feature | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---:|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| 1 | Radial 1–4 (RP, K–M, Gilmore; Tait) | – | · | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ |
| 2 | Radial 5–6 (Mie–Grüneisen) |  | – | ✓ | ✗ | ✗ | ✗ | ✓ | ✓ | ✓ | ✗ | ✓ | ✗ | ✓ |
| 3 | Vapor pressure |  |  | – | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ |
| 4 | Bubble heat transfer |  |  |  | – | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ |
| 5 | Medium heat transfer |  |  |  |  | – | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ |
| 6 | Mass transfer |  |  |  |  |  | – | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ |
| 7 | Kelvin–Voigt stress (1, 2) |  |  |  |  |  |  | – | · | · | · | ✓ | ✓ | ✓ |
| 8 | Maxwell / Jeffreys / Zener (3, 4) |  |  |  |  |  |  |  | – | · | · | ✓ | ✓ | ✓ |
| 9 | UCM / Oldroyd-B, two states (5) |  |  |  |  |  |  |  |  | – | · | ✓ | ✓ | ✓ |
| 10 | UCM / Oldroyd-B, stress field (5, `spectral = 1`) |  |  |  |  |  |  |  |  |  | ✗ | ✗ | ✗ | ✗ |
| 11 | Non-Newtonian viscosity |  |  |  |  |  |  |  |  |  |  | – | ✓ | ✓ |
| 12 | Collapse |  |  |  |  |  |  |  |  |  |  |  | – | ✓ |
| 13 | Pressure waveform |  |  |  |  |  |  |  |  |  |  |  |  | – |

### `f_imr_spectral`

| | Feature | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---:|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| 1 | Radial 1–4 (RP, K–M, Gilmore; Tait) | – | · | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ |
| 2 | Radial 5–6 (Mie–Grüneisen) |  | – | ✓ | ✗ | ✗ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ |
| 3 | Vapor pressure |  |  | – | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ |
| 4 | Bubble heat transfer |  |  |  | – | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ |
| 5 | Medium heat transfer |  |  |  |  | – | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ |
| 6 | Mass transfer |  |  |  |  |  | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ |
| 7 | Kelvin–Voigt stress (1, 2) |  |  |  |  |  |  | – | · | · | · | ✓ | ✗ | ✓ |
| 8 | Maxwell / Jeffreys / Zener (3, 4) |  |  |  |  |  |  |  | – | · | · | ✓ | ✗ | ✓ |
| 9 | UCM / Oldroyd-B, two states (5) |  |  |  |  |  |  |  |  | – | · | ✓ | ✗ | ✓ |
| 10 | UCM / Oldroyd-B, stress field (5, `spectral = 1`) |  |  |  |  |  |  |  |  |  | – | ✓ | ✗ | ✓ |
| 11 | Non-Newtonian viscosity |  |  |  |  |  |  |  |  |  |  | – | ✗ | ✓ |
| 12 | Collapse |  |  |  |  |  |  |  |  |  |  |  | ✗ | ✗ |
| 13 | Pressure waveform |  |  |  |  |  |  |  |  |  |  |  |  | – |

## Requirements

Some features switch on others. A ✓ means the pair runs with these requirements met:
- **Medium heat transfer** requires bubble heat transfer.
- **Mass transfer** requires vapor pressure and bubble heat transfer.
- **Collapse** requires bubble and medium heat transfer, vapor pressure and mass transfer.

Other rules:
- **Stress 3, 4 and 5** require `De > 0`.
- **Stress 5** sets `Ca = Inf`, so it has no elastic contribution.
- **Medium heating with stress 3 to 5** uses the Kelvin–Voigt closed form, which is inconsistent for history-dependent stress, and `f_checkers` warns. The stress field option (feature 10) computes the heating from the stress field itself.

## Why the spectral solver cannot collapse

Collapse only sets the initial state in `f_call_params`: the equilibrium radius, a zero wall velocity and, for stress 3 and 4, the initial stress. Neither solver refers to it afterwards, and the initial-state calculation uses the vapor and gas properties but nothing from mass transfer.

The spectral solver is blocked only because `f_checkers` requires mass transfer for collapse, and the spectral solver has no mass transfer. That requirement dates from the original code. With it relaxed to bubble and medium heat transfer plus vapor pressure, both solvers ran a collapse with stress 0, 1, 3 and 5.

## Planned

`f_imr_nodal` is a planned solver with node values and either the finite-difference or the Chebyshev operator from `f_diff_mat`. Only `f_diff_mat` exists so far.

## How this was checked

Every cell comes from running both solvers with GNU Octave on `main` as of the merge of PR #4. Each run pairs the two features on a Keller–Miksis base case, or radial 5 for the Mie–Grüneisen models, and each ran in its own Octave process. A ✓ means the run finished with a finite `R(t)`.

A separate sweep ran both solvers over every radial model, thermal configuration and stress model: 504 runs, plus 36 more for the remaining thermal combinations. Every non-Newtonian viscosity model, driving waveform, integrator and Jacobian option was also run on a thermal case in both solvers.

Octave needs two stand-ins for MATLAB behavior:
- **No `ode23tb`.** `ode15s` with a dense Jacobian stood in for it.
- **A different `fzero`.** Octave's `fzero` throws an error where MATLAB's returns `NaN`, which breaks the collapse initial state. A stand-in with MATLAB's interval search was used.

One Octave-only problem remains: `f_imr_fd` with `method = 15` crashes Octave's `ode15s` on the sparse Jacobian. MATLAB is not affected.
