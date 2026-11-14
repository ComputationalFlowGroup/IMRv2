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

The enthalpy forms (radial 3–6) describe the liquid with an equation of state in pressure only, at the far-field temperature. Medium heat transfer would change the liquid temperature, so these forms do not support it. Collapse needs medium heat transfer, so these forms do not support collapse either. Bubble heat transfer and mass transfer remain available with radial 3–4. Without medium heat transfer, the wall stays at the far-field temperature. The wall vapor fraction is then the saturated value at that temperature and the bubble pressure. The liquid acts as an infinite heat reservoir. This is consistent with the equation of state, but it neglects the evaporative cooling and condensation heating of the wall. Radial 5–6 remain limited to the polytropic bubble.

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
| 8 | Kelvin–Voigt stress (1, 2, 3, 4) |  |  |  |  |  |  |  | – | · | · | · | ✓ | ✓ | ✓ | · |
| 9 | Maxwell / Jeffreys / Zener (5, 6, 7, 8) |  |  |  |  |  |  |  |  | – | · | · | ✓ | ✓ | ✓ | · |
| 10 | UCM / Oldroyd-B, two states (9) |  |  |  |  |  |  |  |  |  | – | · | ✓ | ✗ | ✓ | · |
| 11 | UCM / Oldroyd-B, stress field (9, `spectral = 1`) |  |  |  |  |  |  |  |  |  |  | ✗ | ✗ | ✗ | ✗ | · |
| 12 | Non-Newtonian viscosity |  |  |  |  |  |  |  |  |  |  |  | – | ✓ | ✓ | ✗ |
| 13 | Collapse |  |  |  |  |  |  |  |  |  |  |  |  | – | ✓ | ✗ |
| 14 | Pressure waveform |  |  |  |  |  |  |  |  |  |  |  |  |  | – | ✗ |
| 15 | PTT / Giesekus, stress field (10, 11, `spectral = 1` forced) |  |  |  |  |  |  |  |  |  |  |  |  |  |  | ✗ |

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
| 8 | Kelvin–Voigt stress (1, 2, 3, 4) |  |  |  |  |  |  |  | – | · | · | · | ✓ | ✗ | ✓ | · |
| 9 | Maxwell / Jeffreys / Zener (5, 6, 7, 8) |  |  |  |  |  |  |  |  | – | · | · | ✓ | ✗ | ✓ | · |
| 10 | UCM / Oldroyd-B, two states (9) |  |  |  |  |  |  |  |  |  | – | · | ✓ | ✗ | ✓ | · |
| 11 | UCM / Oldroyd-B, stress field (9, `spectral = 1`) |  |  |  |  |  |  |  |  |  |  | – | ✓ | ✗ | ✓ | · |
| 12 | Non-Newtonian viscosity |  |  |  |  |  |  |  |  |  |  |  | – | ✗ | ✓ | ✓ |
| 13 | Collapse |  |  |  |  |  |  |  |  |  |  |  |  | ✗ | ✗ | ✗ |
| 14 | Pressure waveform |  |  |  |  |  |  |  |  |  |  |  |  |  | – | ✓ |
| 15 | PTT / Giesekus, stress field (10, 11, `spectral = 1` forced) |  |  |  |  |  |  |  |  |  |  |  |  |  |  | – |

A separate run of the finite difference and Chebyshev schemes of `imr.solve_nodal` gave the same grid. One run is worth knowing about. Rayleigh–Plesset (radial 1) with collapse on 12 Chebyshev points compressed the bubble to `R = 0.005` and `P = 3.6e6`. This is above the ceiling `C1_pv` of the saturation curve. The wall closure then stops with its sign-change error, since no wall vapor fraction satisfies the saturation condition there. Rayleigh–Plesset has no acoustic damping. Run strong collapses with radial 2 instead.

## Requirements

Some features switch on others. A ✓ means the pair runs with these requirements met:

- **Medium heat transfer** requires bubble heat transfer.
- **Mass transfer** requires vapor pressure and bubble heat transfer.
- **Collapse** requires bubble and medium heat transfer, vapor pressure and mass transfer, so only radial 1–2.
- **Medium heat transfer** is not available with the enthalpy forms, radial 3–6.

Other rules:

- **Non-Newtonian viscosity** requires a stress model with a viscous term, stress 1–11. `imr.configuration.validate_case` rejects it with stress 0. It enters every stress model the same way: `-6*intfnu*iDRe` in the stress integral. This term is scaled by `1/DRe` for the Kelvin-Voigt family (stress 1–4) and by `LAM/DRe` for stress 5–11. The medium heating uses the local non-Newtonian viscosity.
- **Stress 5 through 11** require `De >= 1e-4`. Below that threshold, the relaxation is far faster than the bubble motion. The model is then its `De -> 0` limit. `imr.configuration.validate_case` asks for the simpler model instead: stress 1 for stress 5, stress 2 for stress 6, stress 3 for stress 7, stress 4 for stress 8, and stress 1 with `Ca = Inf` for stress 9–11. With a non-Newtonian viscosity, the simpler model takes `DRe/LAM` in place of `DRe`, because stress 5–11 apply it to the solvent share `LAM`.
- **Collapse with stress 5 or 6** starts from the stress of a growth from `Req` to `R = 1` ([`imr.physics.initial_stress`](functions/initial_stress.md)). For stress 6, this growth includes the quadratic term. If no subsonic growth reaches `R = 1` (a very small `Req`), the run stops and asks for `stress0`.
- **Collapse with stress 0–4** starts from rest at the maximum radius. These models carry no stress history. There is no initial stress to set, and the run says so.
- **Collapse with stress 7 through 11** is rejected unless you give `stress0`. The run does not compute the stress at the maximum radius for the 2nd-order Fung or Mooney-Rivlin Zener springs, UCM/Oldroyd-B, PTT or Giesekus. Starting from zero stress would be silently wrong for these models. `imr.physics.initial_stress` does not extend to the new springs' Zener variants (7, 8); their Kelvin-Voigt variants (3, 4) need no such extension. `imr.solve_spectral` does not support collapse at all (feature 13's diagonal is ✗), so this rule only matters if that restriction is ever relaxed.
- **Stress 9, 10 and 11** set `Ca = Inf`, so they have no elastic contribution.
- **Stress 10 (PTT) and 11 (Giesekus)** generalize stress 9 with the nonlinear terms of [Warnez and Johnsen (2015)](references.md); see [stress and viscosity](../model/stress-viscosity.md#ptt-and-giesekus-stress-10-11). Neither has a closed-form, two-ODE reduction. `imr.configuration.validate_case` rejects them outright in `imr.solve_nodal`. In `imr.solve_spectral` it forces `spectral = 1`, with a warning if `spectral` was 0.
- **Stress 3 and 7** add the next term of the series that stress 2 and 6 already truncate after their first order. They reuse `alphax`; see [2nd-order Fung](../model/stress-viscosity.md#2nd-order-fung-stress-3-7). **Stress 4 and 8** add a second elastic modulus `C2` to the neo-Hookean spring; see [Mooney-Rivlin](../model/stress-viscosity.md#mooney-rivlin-stress-4-8). `alphax = 0` (stress 3, 7) or `C2 = 0` (stress 4, 8) reduces exactly to the corresponding neo-Hookean model (1, 5).
- **Medium heating with stress 5 to 11** comes from the local polymer stress, not the Kelvin–Voigt form (see [`imr.physics.stress_power`](functions/stress_power.md)).
    - **Stress 5, 6, 7 and 8:** carry that stress as a field next to `Z1` ([`imr.physics.zener_field`](functions/zener_field.md)).
    - **Stress 9:** the two-state form carries only integrals of the stress. With medium heat transfer, `imr.solve_spectral` switches to the stress field (feature 11) with a warning, and `imr.solve_nodal` rejects the run: its UCM / Oldroyd-B with medium heat transfer cell is ✗.
    - **Stress 10 and 11:** always run as the spectral stress field (feature 15), so medium heat transfer uses the local polymer stress the same way feature 11 does.

## Historical investigation: spectral collapse

The following paragraph preserves a historical experiment with relaxed checks. It is not a supported feature. It is not a run performed for this documentation change. Current `imr.configuration.validate_case` still rejects spectral collapse.

Collapse only sets the initial state in `imr.configuration.prepare_case`: the equilibrium radius, a zero wall velocity and, for stress 5 and 6, the initial stress (stress 7, 8 and 9 need `stress0`). No solver refers to it afterwards. The initial-state calculation uses the vapor and gas properties, but nothing from mass transfer.

`imr.solve_spectral` is blocked for one reason only: `imr.configuration.validate_case` requires mass transfer for collapse, and `imr.solve_spectral` has no mass transfer. That requirement dates from the original code. The historical experiment relaxed it to bubble and medium heat transfer plus vapor pressure. With that relaxed requirement, `imr.solve_nodal` and `imr.solve_spectral` ran a collapse with the no-stress, neo-Hookean Kelvin-Voigt, neo-Hookean Zener and UCM/Oldroyd-B models (stress 0, 1, 5 and 9 in the current numbering).

## Historical matrix checks

These are previously reported checks, not a fresh execution of the matrix. A static ✓ mark describes intended compatibility, subject to every listed prerequisite. A ✓ mark is not a guarantee of success for arbitrary inputs.

Each cell pairs two features on a Keller–Miksis base case, radial 2, in both solvers and both schemes of `imr.solve_nodal`. Each radial model (1–6) was paired with every other feature. Every member of a feature group was run: the neo-Hookean and quadratic neo-Hookean Kelvin-Voigt and Zener models (stress 1, 2, 5 and 6 in the current numbering, with `alphax = 0.3` for the quadratic ones). Non-Newtonian runs use Kelvin–Voigt stress unless the other feature sets a stress model. A ✓ means every run of the cell finished with a finite `R(t)`. Every ✗ is an `INPUT ERROR` from `imr.configuration.validate_case`. The one other failure is the Rayleigh–Plesset collapse run described above.

The test suites keep the matrix honest: [`f_test_cases`](../verification/tests.md) builds every suite from the combinations `imr.configuration.validate_case` accepts, and [`s_test_suite_permutations`](../verification/permutations.md) runs its enumerated, accepted option combinations briefly and fails on any error that is not an `INPUT ERROR`.
