# MATLAB API

Use the public functions below from MATLAB with `src/` on the path. All solver inputs and outputs are dimensionless; convert physical properties before calling a solver.

| Function | Purpose |
|---|---|
| [`imr.solve`](solve.md) | Common entry point; choose nodal or spectral with `solver` |
| [`imr.solve_nodal`](solve_nodal.md) | Nodal solver with finite-difference or Chebyshev operators |
| [`imr.solve_spectral`](solve_spectral.md) | Spectral solver with optional local stress field |
| [`imr.default_inputs`](default_inputs.md) | Default input struct |
| [`imr.fluid_properties`](fluid_properties.md) | Shared liquid, vapor and saturation evaluation |
| [Fluid profiles](../model/fluids.md) | MATLAB material data files, independent of simulation choices |
| [`imr.nondimensionalize`](nondimensionalize.md) | SI-to-dimensionless conversion and reference scales |
| [`imr.sweep`](sweep.md) | Full-grid parameter study |

Start with the [complete worked call](../getting-started/quick-start.md). Read [configuration](../running/configuration.md) for precedence and [outputs](../running/output.md) for saving and replaying a run. The [function call map](function-call-map.md) links internal helpers; they are not additional solver options.

IMR has no plugin registry. [Extending IMR](extending.md) describes direct changes to defaults, parsing, checks and equations.
