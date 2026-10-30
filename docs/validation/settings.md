# Settings and initialization audit

`tests/suites/s_test_suite_settings.m` declares cases that must be accepted or rejected. Unlike the generated permutation suite, it does not filter cases through the checker first. A regression that rejects a supported boundary therefore fails this suite.

```bash
cd tests
matlab -batch "imr_test_setup; s_test_suite_settings"
```

## Corrected failure paths

| Trigger | Previous behavior | Correction |
|---|---|---|
| Either public solver after the namespace move | Unqualified `prepare_case` could not resolve | Call `imr.configuration.prepare_case` |
| Nonscalar, complex or nonfinite scalar parameters | Some arithmetic ran before validation; `v_a = Inf` could pass | Validate primitive values first and use floating-point arithmetic |
| Integer-valued numeric inputs | MATLAB integer divisions could round derived properties | Normalize effective settings, time and initial stress to `double` |
| `p0star = 0`, `vapor = 0`, `bubtherm = 1` | Thermal equations divided by zero pressure | Require positive initial total pressure for thermal bubbles |
| `p0star = 0`, `masstrans = 1` | Pure-vapor limit divided by a zero gas fraction | Require positive, representable noncondensible gas fraction |
| Frozen thermal mixture with vapor and `rv = 0` or `ra = 0` | Mixture fractions were invalid even without mass diffusion | Check gas constants whenever an active vapor mixture uses them |
| `stress0 = NaN`, complex values or a matrix | Stress-history length alone could pass | Require a finite, real vector with the appropriate length |
| Very small conductivity slope | Cancellation in temperature inversion and wall reconstruction | Rationalized inverse and a wall quadratic in `Tw - 1` |
| Spectral analytic Jacobian with medium heat and stress 0–2 | Tried to expand an empty polymer basis | Form the polymer contribution only when that field exists |
| Disabled excess viscosity (`DRe = 0` or zero solvent share) with large law parameters | An unused viscosity calculation could overflow before multiplication by zero | Skip the inactive contribution in RHS and closure calculations |
| Sonic initial denominator or pressure outside the liquid EoS domain | Could enter integration with an infinite or complex derivative | Check the assembled initial derivative once before integration |
| Collapse initialization without a positive stable equilibrium | Scalar-start root finding could probe negative vapor mass or retry without a bound | Bracket the positive equilibrium radius and reject cases without a stable root |
| Finite radius with nonfinite pressure, velocity or temperature | Shared test runner inspected only radius | Inspect all seven physical outputs for real, finite values |

`eps3` enters only the Giesekus model (`stress = 7`); for every other stress it remains reserved, with no implemented equation term. A nonzero value raises `IMR:unusedParameter` unless `stress = 7`, and `stress = 7` with `eps3 = 0` (identical to stress 5) raises `IMR:parameters`, so both are visible during setup.

## Boundaries preserved

The suite explicitly checks all three discretizations, compatible radial/thermal families, supported Jacobian options, minimum active thermal grids, inactive grids of size zero, `kappa = 1`, and zero vapor pressure. An inactive vapor component does not enter the mixture fraction. A frozen pure-vapor mixture without mass transfer remains accepted when its conductivity is positive, including constant vapor conductivity. The gas wall-closure implementation still requires `alpha_g > 0` for thermal models.

`stress = 6` (PTT) and `7` (Giesekus) are rejected in `imr.solve_nodal` and forced to the spectral stress field in `imr.solve_spectral`. Giesekus at `eps3 = 0` is checked to give the exact same right-hand side, Jacobian and initial state as `stress = 5` (UCM/Oldroyd-B), not only a close one.

All seven non-Newtonian laws are checked at zero wall shear with each viscous stress family. Separate cases ensure a disabled viscosity law cannot contaminate the solution through overflow. Known temperatures are mapped forward and inverted to check the conductivity limit independently.

## Scope and execution status

The assembled-state check performs one RHS evaluation during setup, outside the time integration loop. Scalar compatibility checks alone cannot determine every liquid-domain or wall-boundary failure. This check reports `IMR:initialState` for nonfinite or complex initial derivatives; programming exceptions propagate normally.

An admissible initial state does not guarantee completion for arbitrary forcing. A later collapse can approach zero radius, cross a sonic denominator or leave the EoS/temperature domain. Grid convergence, tolerances and physical validity still require case-specific assessment. The [permutation suite](permutations.md) checks final-time completion for its selected cases; this audit does not establish exhaustive coverage of continuous parameter ranges.

MATLAB is unavailable in this cloud environment. The MATLAB assertions and CI job are prepared but have not been executed here. Static syntax/lint checks and independent arithmetic checks do not substitute for running the suites in MATLAB.
