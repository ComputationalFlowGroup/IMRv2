# MATLAB API

Use the public functions below from MATLAB with `src/` on the path. All solver inputs and outputs are dimensionless; convert physical properties before calling a solver.

| Function | Purpose |
|---|---|
| [`imr.solve_nodal`](f_imr_nodal.md) | Nodal solver with finite-difference or Chebyshev operators |
| [`imr.solve_spectral`](f_imr_spectral.md) | Spectral solver with optional local stress field |
| [`imr.default_inputs`](f_imr_inputs.md) | Default input struct |
| [`imr.fluid_properties`](f_fluid_properties.md) | Shared liquid, vapor and saturation evaluation |
| [Fluid profiles](../model/fluids.md) | MATLAB material data files, independent of simulation choices |
| [`imr.nondimensionalize`](f_dim2nondim.md) | SI-to-dimensionless conversion and reference scales |
| [`imr.sweep`](f_imr_sweep.md) | Full-grid parameter study |

Start with the [complete worked call](../getting-started/quick-start.md). Read [configuration](../running/configuration.md) for precedence and [outputs](../running/output.md) for saving and replaying a run. The [code map](../reference/code-map.md) links internal helpers; they are not additional solver options.

IMR has no plugin registry. [Extending IMR](extending.md) describes direct changes to defaults, parsing, checks and equations.
