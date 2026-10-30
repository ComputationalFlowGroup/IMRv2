# Model compatibility

These grids show which features can run together in each solver: `imr.solve_nodal`, with either differentiation scheme, and `imr.solve_spectral`. Every feature appears on both axes, and only the upper triangle is filled in.

Legend:

- **✓** the two features run together.
- **✗** they cannot run together. `imr.configuration.validate_case` rejects the run with an `INPUT ERROR`.
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

### `imr.solve_nodal`, either scheme (finite difference or Chebyshev)

| | Feature | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---:|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| 1 | Radial 1–2 | – | · | · | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✗ |
| 2 | Radial 3–4 |  | – | · | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✗ | ✓ | ✗ |
| 3 | Radial 5–6 |  |  | – | ✓ | ✗ | ✗ | ✗ | ✓ | ✓ | ✓ | ✗ | ✓ | ✗ | ✓ | ✗ |
| 4 | Vapor pressure |  |  |  | – | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✗ |
| 5 | Bubble heat transfer |  |  |  |  | – | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✗ |
| 6 | Medium heat transfer |  |  |  |  |  | – | ✓ | ✓ | ✓ | ✗ | ✗ | ✓ | ✓ | ✓ | ✗ |
| 7 | Mass transfer |  |  |  |  |  |  | – | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✗ |
| 8 | Kelvin–Voigt stress (1, 2, 8, 10) |  |  |  |  |  |  |  | – | · | · | · | ✓ | ✓ | ✓ | · |
| 9 | Maxwell / Jeffreys / Zener (3, 4, 9, 11) |  |  |  |  |  |  |  |  | – | · | · | ✓ | ✓ | ✓ | · |
| 10 | UCM / Oldroyd-B, two states (5) |  |  |  |  |  |  |  |  |  | – | · | ✓ | ✗ | ✓ | · |
| 11 | UCM / Oldroyd-B, stress field (5, `spectral = 1`) |  |  |  |  |  |  |  |  |  |  | ✗ | ✗ | ✗ | ✗ | · |
| 12 | Non-Newtonian viscosity |  |  |  |  |  |  |  |  |  |  |  | – | ✓ | ✓ | ✗ |
| 13 | Collapse |  |  |  |  |  |  |  |  |  |  |  |  | – | ✓ | ✗ |
| 14 | Pressure waveform |  |  |  |  |  |  |  |  |  |  |  |  |  | – | ✗ |
| 15 | PTT / Giesekus, stress field (6, 7, `spectral = 1` forced) |  |  |  |  |  |  |  |  |  |  |  |  |  |  | ✗ |

### `imr.solve_spectral`

| | Feature | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---:|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| 1 | Radial 1–2 | – | · | · | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ |
| 2 | Radial 3–4 |  | – | · | ✓ | ✓ | ✗ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ |
| 3 | Radial 5–6 |  |  | – | ✓ | ✗ | ✗ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ |
| 4 | Vapor pressure |  |  |  | – | ✓ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ |
| 5 | Bubble heat transfer |  |  |  |  | – | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ |
| 6 | Medium heat transfer |  |  |  |  |  | – | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ |
| 7 | Mass transfer |  |  |  |  |  |  | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ |
| 8 | Kelvin–Voigt stress (1, 2, 8, 10) |  |  |  |  |  |  |  | – | · | · | · | ✓ | ✗ | ✓ | · |
| 9 | Maxwell / Jeffreys / Zener (3, 4, 9, 11) |  |  |  |  |  |  |  |  | – | · | · | ✓ | ✗ | ✓ | · |
| 10 | UCM / Oldroyd-B, two states (5) |  |  |  |  |  |  |  |  |  | – | · | ✓ | ✗ | ✓ | · |
| 11 | UCM / Oldroyd-B, stress field (5, `spectral = 1`) |  |  |  |  |  |  |  |  |  |  | – | ✓ | ✗ | ✓ | · |
| 12 | Non-Newtonian viscosity |  |  |  |  |  |  |  |  |  |  |  | – | ✗ | ✓ | ✓ |
| 13 | Collapse |  |  |  |  |  |  |  |  |  |  |  |  | ✗ | ✗ | ✗ |
| 14 | Pressure waveform |  |  |  |  |  |  |  |  |  |  |  |  |  | – | ✓ |
| 15 | PTT / Giesekus, stress field (6, 7, `spectral = 1` forced) |  |  |  |  |  |  |  |  |  |  |  |  |  |  | – |

The finite difference and Chebyshev schemes of `imr.solve_nodal` were run separately and gave the same grid. One run is worth knowing about: Rayleigh–Plesset (radial 1) with collapse on 12 Chebyshev points compressed the bubble to `R = 0.005` and `P = 3.6e6`, above the ceiling `C1_pv` of the saturation curve. The wall closure then stops with its sign-change error, since no wall vapor fraction satisfies the saturation condition there. Rayleigh–Plesset has no acoustic damping, so strong collapses are best run with radial 2.

## Requirements

Some features switch on others. A ✓ means the pair runs with these requirements met:

- **Medium heat transfer** requires bubble heat transfer.
- **Mass transfer** requires vapor pressure and bubble heat transfer.
- **Collapse** requires bubble and medium heat transfer, vapor pressure and mass transfer, so only radial 1–2.
- **Medium heat transfer** is not available with the enthalpy forms, radial 3–6.

Other rules:

- **Non-Newtonian viscosity** requires a stress model with a viscous term, stress 1–11; `imr.configuration.validate_case` rejects it with stress 0. It enters every stress model the same way, `-6*intfnu*iDRe` in the stress integral, scaled by `1/DRe` for the Kelvin-Voigt family (stress 1, 2, 8, 10) and `LAM/DRe` for stress 3–7, 9, 11, and the medium heating uses the local non-Newtonian viscosity.
- **Stress 3 through 7, 9 and 11** require `De >= 1e-4`. Below that the relaxation is far faster than the bubble motion, the model is its `De -> 0` limit, and `imr.configuration.validate_case` asks for the simpler model instead: stress 1 for stress 3, stress 2 for stress 4, stress 8 for stress 9, stress 10 for stress 11, and stress 1 with `Ca = Inf` for stress 5–7. With a non-Newtonian viscosity, the simpler model takes `DRe/LAM` in place of `DRe`, since stress 3–7, 9, 11 apply it to the solvent share `LAM`.
- **Collapse with stress 3 or 4** starts from the stress of a growth from `Req` to `R = 1` ([`imr.physics.initial_stress`](functions/f_init_stress.md)), including the quadratic term for stress 4. If no subsonic growth reaches `R = 1` (very small `Req`), the run stops and asks for `stress0`.
- **Collapse with stress 0–2, 8 or 10** starts from rest at the maximum radius. These models carry no stress history, so there is no initial stress to set, and the run says so.
- **Collapse with stress 5, 6, 7, 9 or 11** is rejected unless `stress0` is given: the stress at the maximum radius is not computed for UCM/Oldroyd-B, PTT, Giesekus, or the 2nd-order Fung or Mooney-Rivlin Zener springs, and starting from zero stress would be silently wrong. `imr.physics.initial_stress` is not extended to the new springs' Zener variants (9, 11), unlike their Kelvin-Voigt variants (8, 10), which need no such extension. Collapse is unsupported in `imr.solve_spectral` regardless (feature 13's diagonal is ✗), so this only matters if that restriction is ever relaxed.
- **Stress 5, 6 and 7** set `Ca = Inf`, so they have no elastic contribution.
- **Stress 6 (PTT) and 7 (Giesekus)** generalize stress 5 with the nonlinear terms of [Warnez and Johnsen (2015)](references.md); see [stress and viscosity](../model/stress-viscosity.md#ptt-and-giesekus-stress-6-7). Neither has a closed-form, two-ODE reduction, so `imr.configuration.validate_case` rejects them in `imr.solve_nodal` outright and forces `spectral = 1` in `imr.solve_spectral` (with a warning if it was 0).
- **Stress 8 and 9** add the next term of the series stress 2 and 4 already truncate after their first order, reusing `alphax`; see [2nd-order Fung](../model/stress-viscosity.md#2nd-order-fung-stress-8-9). **Stress 10 and 11** add a second elastic modulus `C2` to the neo-Hookean spring; see [Mooney-Rivlin](../model/stress-viscosity.md#mooney-rivlin-stress-10-11). `alphax = 0` (stress 8, 9) or `C2 = 0` (stress 10, 11) reduces exactly to the corresponding neo-Hookean model (1, 3).
- **Medium heating with stress 3 to 7, 9 and 11** comes from the local polymer stress, not the Kelvin–Voigt form (see [`imr.physics.stress_power`](functions/f_stress_dissipation.md)).
    - **Stress 3, 4, 9 and 11:** carry that stress as a field next to `Z1` ([`imr.physics.zener_field`](functions/f_zener_field.md)).
    - **Stress 5:** the two-state form carries only integrals of the stress. With medium heat transfer, `imr.solve_spectral` switches to the stress field (feature 11) with a warning, and `imr.solve_nodal` rejects the run: its UCM / Oldroyd-B with medium heat transfer cell is ✗.
    - **Stress 6 and 7:** always run as the spectral stress field (feature 15), so medium heat transfer uses the local polymer stress the same way feature 11 does.

## Historical investigation: spectral collapse

The following paragraph preserves a historical experiment with relaxed checks, not a supported feature or a run performed for this documentation change. Current `imr.configuration.validate_case` still rejects spectral collapse.

Collapse only sets the initial state in `imr.configuration.prepare_case`: the equilibrium radius, a zero wall velocity and, for stress 3 and 4, the initial stress (stress 5, 9 and 11 need `stress0`). No solver refers to it afterwards, and the initial-state calculation uses the vapor and gas properties but nothing from mass transfer.

`imr.solve_spectral` is blocked only because `imr.configuration.validate_case` requires mass transfer for collapse, and it has no mass transfer. That requirement dates from the original code. With it relaxed to bubble and medium heat transfer plus vapor pressure, `imr.solve_nodal` and `imr.solve_spectral` ran a collapse with stress 0, 1, 3 and 5.

## Historical matrix checks

These are previously reported checks, not a fresh execution of the matrix. Static ✓ marks describe intended compatibility subject to every listed prerequisite, not guaranteed success for arbitrary inputs.

Each cell pairs two features on a Keller–Miksis base case, radial 2, in both solvers and both schemes of `imr.solve_nodal`. Each radial model (1–6) was paired with every other feature, and every member of a feature group was run: stress 1 and 2, and stress 3 and 4 (with `alphax = 0.3` for the quadratic models). Non-Newtonian runs use Kelvin–Voigt stress unless the other feature sets a stress model. A ✓ means every run of the cell finished with a finite `R(t)`, and every ✗ is an `INPUT ERROR` from `imr.configuration.validate_case`. The one other failure is the Rayleigh–Plesset collapse run described above.

The test suites keep the matrix honest: [`f_test_cases`](../validation/tests.md) builds every suite from the combinations `imr.configuration.validate_case` accepts, and [`s_test_suite_permutations`](../validation/permutations.md) runs its enumerated, accepted option combinations briefly and fails on any error that is not an `INPUT ERROR`.
