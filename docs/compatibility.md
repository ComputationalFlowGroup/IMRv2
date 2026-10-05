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
| 6 | Medium heat transfer |  |  |  |  |  | – | ✓ | ✓ | ✓ | ✗ | ✗ | ✓ | ✓ | ✓ |
| 7 | Mass transfer |  |  |  |  |  |  | – | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ |
| 8 | Kelvin–Voigt stress (1, 2) |  |  |  |  |  |  |  | – | · | · | · | ✓ | ✓ | ✓ |
| 9 | Maxwell / Jeffreys / Zener (3, 4) |  |  |  |  |  |  |  |  | – | · | · | ✓ | ✓ | ✓ |
| 10 | UCM / Oldroyd-B, two states (5) |  |  |  |  |  |  |  |  |  | – | · | ✓ | ✗ | ✓ |
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

- **Non-Newtonian viscosity** requires a stress model with a viscous term, stress 1–5; `f_checkers` rejects it with stress 0. It enters every stress model the same way, `-6*intfnu*iDRe` in the stress integral, scaled by `1/DRe` for stress 1–2 and `LAM/DRe` for stress 3–5, and the medium heating uses the local non-Newtonian viscosity.
- **Stress 3, 4 and 5** require `De >= 1e-4`. Below that the relaxation is far faster than the bubble motion, the model is its `De -> 0` limit, and `f_checkers` asks for the simpler model instead: stress 1 for stress 3, stress 2 for stress 4, and stress 1 with `Ca = Inf` for stress 5. With a non-Newtonian viscosity, the simpler model takes `DRe/LAM` in place of `DRe`, since stress 3–5 apply it to the solvent share `LAM`.
- **Collapse with stress 3 or 4** starts from the stress of a growth from `Req` to `R = 1` ([`f_init_stress`](src/f_init_stress.md)), including the quadratic term for stress 4. If no subsonic growth reaches `R = 1` (very small `Req`), the run stops and asks for `stress0`.
- **Collapse with stress 0–2** starts from rest at the maximum radius. These models carry no stress history, so there is no initial stress to set, and the run says so.
- **Collapse with stress 5** is rejected unless `stress0` is given: the stress at the maximum radius is not computed for UCM and Oldroyd-B, and starting from zero stress would be silently wrong.
- **Stress 5** sets `Ca = Inf`, so it has no elastic contribution.
- **Medium heating with stress 3 to 5** comes from the local polymer stress, not the Kelvin–Voigt form (see [`f_stress_dissipation`](src/f_stress_dissipation.md)).
    - **Stress 3 and 4:** carry that stress as a field next to `Z1` ([`f_zener_field`](src/f_zener_field.md)).
    - **Stress 5:** the two-state form carries only integrals of the stress. With medium heat transfer, `f_imr_spectral` switches to the stress field (feature 11) with a warning, and `f_imr_nodal` rejects the run: its UCM / Oldroyd-B with medium heat transfer cell is ✗.

## Why the spectral solver cannot collapse

Collapse only sets the initial state in `f_call_params`: the equilibrium radius, a zero wall velocity and, for stress 3 and 4, the initial stress (stress 5 needs `stress0`). No solver refers to it afterwards, and the initial-state calculation uses the vapor and gas properties but nothing from mass transfer.

`f_imr_spectral` is blocked only because `f_checkers` requires mass transfer for collapse, and it has no mass transfer. That requirement dates from the original code. With it relaxed to bubble and medium heat transfer plus vapor pressure, `f_imr_nodal` and `f_imr_spectral` ran a collapse with stress 0, 1, 3 and 5.

## How this was checked

Each cell pairs two features on a Keller–Miksis base case, radial 2, in both solvers and both schemes of `f_imr_nodal`. Each radial model (1–6) was paired with every other feature, and every member of a feature group was run: stress 1 and 2, and stress 3 and 4 (with `alphax = 0.3` for the quadratic models). Non-Newtonian runs use Kelvin–Voigt stress unless the other feature sets a stress model. A ✓ means every run of the cell finished with a finite `R(t)`, and every ✗ is an `INPUT ERROR` from `f_checkers`. The one other failure is the Rayleigh–Plesset collapse run described above.

The test suites keep the matrix honest: [`f_test_cases`](toolchain/f_test_cases.md) builds every suite from the combinations `f_checkers` accepts, and [`s_test_suite_permutations`](toolchain/s_test_suite_permutations.md) runs every accepted combination of the solver options briefly and fails on any error that is not an `INPUT ERROR`.
