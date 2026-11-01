# Examples

Run `imr_setup; run_water_case` from the checkout root for the small, namespaced water example. See [namespaces and migration](../reference/namespaces.md#run-a-case).

For configuring a run (struct vs. name/value pairs) see [configuration](../running/configuration.md); for choosing nodal vs. spectral or finite differences vs. Chebyshev see [choose a solver](../running/solvers.md); for sweeping inputs and saving outputs see [parameter sweeps](../running/sweeps.md) and [outputs](../running/output.md#save-and-replay-a-run). Each fact and recipe lives on its own workflow page.

## Bundled script inventory

| File | What it actually does | Limits |
|---|---|---|
| `examples/s_run_fdmodel.m` | Polytropic nodal Keller–Miksis case with vapor pressure and viscous stress | Does not enable heat or mass transfer |
| `examples/s_run_full_model.m` | Gaussian-driven polytropic Kelvin–Voigt case | Despite the name, bubble heat, medium heat and mass transfer are off |
| `examples/f_multirun.m` | Five large material studies | Expensive parallel study, described in [sweeps](../running/sweeps.md) |
| `examples/matlab-submitjob.sh` | Submits the large study on a cluster | Site-specific `sbatch` and `matlab-threaded` commands |
| `examples/f_contour_field.m` | Legacy stress-contour postprocessing fragment | Needs existing workspace variables and output directories; not standalone |

Run path setup once per MATLAB session, from the checkout root. The setup adds the public source and examples without recursive path expansion:

```matlab
% Run the bundled Gaussian-driven example.
imr_setup;
s_run_full_model;
```

The contour fragment reads `R`, `U`, `f_r`, `max_tau_r`, `lm` and other workspace values. It is retained for provenance, not offered as a complete recipe. The scaling benchmarks also contain legacy data; see the [repository map](../reference/repository.md).
