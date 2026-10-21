# MATLAB API

Use the public functions below from MATLAB with `src/` on the path. All solver inputs and outputs are dimensionless; convert physical properties before calling a solver.

| Function | Purpose |
|---|---|
| [`f_imr_nodal`](f_imr_nodal.md) | Nodal solver with finite-difference or Chebyshev operators |
| [`f_imr_spectral`](f_imr_spectral.md) | Spectral solver with optional local stress field |
| [`f_imr_inputs`](f_imr_inputs.md) | Default input struct |
| [Mercury properties](f_mercury_eos.md) | Liquid/vapor fit coefficients and SI evaluation in base MATLAB |
| [`f_mercury_inputs`](f_mercury_inputs.md) | Ambient mercury preparation for pressure-form KM |
| [`f_dim2nondim`](f_dim2nondim.md) | SI-to-dimensionless conversion and reference scales |
| [`f_imr_sweep`](f_imr_sweep.md) | Full-grid parameter study |

Start with the [complete worked call](../getting-started/quick-start.md). Read [configuration](../running/configuration.md) for precedence and [outputs](../running/output.md) for saving and replaying a run. The [code map](../reference/code-map.md) links internal helpers; they are not additional solver options.

IMR has no plugin registry. [Extending IMR](extending.md) describes direct changes to defaults, parsing, checks and equations.
