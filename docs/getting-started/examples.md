# Examples

Run `imr_setup; run_water_case` from the checkout root for the small, namespaced water example. See [namespaces and migration](../reference/namespaces.md#run-a-case).

These recipes use the public MATLAB functions from the checkout root. Their interfaces were checked against source; execution and numerical results remain pending local MATLAB validation.

## Configure a run as a struct

```matlab
% Start from dimensionless defaults and change two model choices.
imr_setup;
p = imr.default_inputs();
p.radial = 2;
p.stress = 1;
[t, R, ~, ~, ~, ~, ~, used] = imr.solve_nodal(p);
```

## Switch spatial discretization

```matlab
% Run the nodal solver with Chebyshev differentiation.
imr_setup;
[t, R] = imr.solve_nodal('radial',2,'stress',1,'bubtherm',1,'scheme',1);
```

## Sweep inputs and save outputs

Use the complete recipes in [parameter sweeps](../running/sweeps.md) and [outputs](../running/output.md#save-and-replay-a-run). Each fact and recipe lives on its own workflow page.

## Bundled script inventory

| File | What it actually does | Limits |
|---|---|---|
| `examples/s_run_fdmodel.m` | Polytropic nodal Keller–Miksis case with vapor pressure and viscous stress | Does not enable heat or mass transfer |
| `examples/s_run_full_model.m` | Gaussian-driven polytropic Kelvin–Voigt case | Despite the name, bubble heat, medium heat and mass transfer are off |
| `examples/s_run_mercury.m` | Pressure-form KM with fitted Hg ambient density/acoustics and equilibrium vapor | Explicit inviscid, zero-surface-tension demonstration, not a measured transport parameter set |
| `examples/f_multirun.m` | Five large material studies | Expensive parallel study, described in [sweeps](../running/sweeps.md) |
| `examples/matlab-submitjob.sh` | Submits the large study on a cluster | Site-specific `sbatch` and `matlab-threaded` commands |
| `examples/f_contour_field.m` | Legacy stress-contour postprocessing fragment | Needs existing workspace variables and output directories; not standalone |

Run path setup once per MATLAB session, from the checkout root. The setup adds the public source and examples without recursive path expansion:

```matlab
% Run the bundled Gaussian-driven example.
imr_setup;
s_run_full_model;
```

The contour fragment reads `R`, `U`, `f_r`, `max_tau_r`, `lm` and other workspace values. It is retained for provenance, not offered as a complete recipe. The scaling benchmarks also contain legacy data; see [reported comparisons](../validation/convergence.md#legacy-scaling-data).
