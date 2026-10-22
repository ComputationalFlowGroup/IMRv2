# Collapse initial conditions

`collapse = 1` initializes the bubble at maximum radius with zero wall velocity. `imr.configuration.prepare_case` constructs that state before time integration; it is not a separate radial equation.

Use [compatibility](../reference/compatibility.md) for the exact prerequisites. The current checker requires bubble heat, medium heat, vapor pressure and mass transfer, so supported collapse uses the nodal solver with radial 1 or 2.

| Stress | Initial history | User action |
|---|---|---|
| 0–2 | No stress history state | Start from rest at maximum radius |
| 3–4 | `imr.physics.initial_stress` shoots a subsonic growth trajectory unless `stress0` is supplied | Provide `stress0` if no admissible growth reaches maximum radius |
| 5 | No automatic maximum-radius stress calculation | Supply the two-state `stress0`; medium heat transfer remains incompatible in the nodal solver |

The stress-5 requirement is necessary but does not remove its medium-heating incompatibility. A full supported collapse configuration with stress 5 is therefore not currently available through the checker; do not infer support from the initial-stress rule alone.

For stress 4, the initial stress calculation includes `alphax`; stress 3 uses zero for that correction. The approximate standalone helper is not called by the solver. See [initial-stress method](../method/initial-stress.md) and [troubleshooting](troubleshooting.md).

The recomputed equilibrium radius and zero initial wall velocity are setup quantities. The existing historical report that the recomputed radius reproduces the input over selected cases is preserved in [reported comparisons](../validation/convergence.md#reported-collapse-setup-comparison), not presented as a new run.
