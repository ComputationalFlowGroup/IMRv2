# Namespaces, organization and migration

Public MATLAB functions live in the `imr` namespace. Add only the parent `src` directory through `imr_setup`; MATLAB resolves the `+` package hierarchy automatically. There is no `private` folder and no package loader or registry.

```text
imr_setup.m
src/
  +imr/
    solve.m
    solve_nodal.m
    solve_spectral.m
    default_inputs.m
    nondimensionalize.m
    fluid_properties.m
    sweep.m
    approximate_initial_stress.m
    +configuration/   defaults, schema, preparation, validation
    +solvers/         shared runner, builders, state layout, display
    +physics/         radial, stress, thermal and vapor closures
    +numerics/        integration, operators and Jacobians
    +fluids/          material profiles
  data/              measured pressure splines
  f_imr_nodal.m      compatibility adapters, no numerical implementation
  ...
examples/
  run_water_case.m
```

Each solver implementation exists once, in `imr.solvers.build_nodal` or `imr.solvers.build_spectral`. The public functions delegate to `imr.solvers.run`, which shares case preparation, integration and reconstruction. Tightly coupled wall helpers remain local functions in their builder. Package functions are callable for developer tests; a namespace organizes names rather than restricting access.

## Run a case

From the checkout root in MATLAB:

```matlab
imr_setup
[args,ref] = imr.nondimensionalize('fluid',@imr.fluids.water, ...
    'r0',50e-6,'req',50e-6,'u0',0.1,'g',0, ...
    'radial',2,'stress',1,'bubtherm',0,'medtherm',0, ...
    'vapor',0,'masstrans',0,'tfin',20e-6);
[t,R,U,P] = imr.solve(args{:},'solver','nodal');
plot(t*ref.t0*1e6,R*ref.R0*1e6)
xlabel('Time (microseconds)')
ylabel('Radius (micrometers)')
```

The bundled `run_water_case` script performs this call. From a terminal at the checkout root:

```bash
matlab -batch "imr_setup; run_water_case"
```

The small initial velocity starts oscillation about the equilibrium radius. Water properties use the documented legacy reference set; this example does not introduce a new calibrated water EoS. Change the `solver` value to `'spectral'` for the spectral solver. Both existing solver wrappers remain supported. See the [Keller–Miksis example](../getting-started/keller-miksis.md) for SI scaling and the [API](../api/index.md) for output shapes.

## Compatibility and naming

| Existing public name | Preferred name |
|---|---|
| `f_imr_nodal` | `imr.solve_nodal` |
| `f_imr_spectral` | `imr.solve_spectral` |
| `f_imr_inputs` | `imr.default_inputs` |
| `f_dim2nondim` | `imr.nondimensionalize` |
| `f_imr_sweep` | `imr.sweep` |
| `f_fluid_properties` | `imr.fluid_properties` |
| `f_approx_init_stress` | `imr.approximate_initial_stress` |
| `default_case` | `imr.configuration.defaults` |
| `f_call_params` | `imr.configuration.legacy_options` |
| `fluids.water`, `fluids.mercury` | `imr.fluids.water`, `imr.fluids.mercury` |

These old public names remain small forwarding adapters with unchanged inputs/output order and no repeated implementation. Internal helpers have moved to their responsibility package; code that called an internal helper directly must use its qualified name from the [code map](code-map.md). Do not `import imr.*` in library code; explicit names avoid ambiguity and make dependencies visible.

Use lowercase `snake_case`, descriptive quantities or verbs, and a primary function matching the filename. For example, `imr.physics.radial_acceleration` is defined as `radial_acceleration` in `src/+imr/+physics/radial_acceleration.m`. Do not put dots in its function declaration or add `+imr` to the search path itself.

## Runtime cost and validation

The organization adds no path setup, file reads or validation inside numerical callbacks. Existing kernels retain their equations and call counts; their names become qualified. The shared runner executes once per simulation. Legacy adapters are not used by the namespaced RHS or Jacobians.

MATLAB name resolution and dispatch performance can vary by release. Source consistency and static checks cannot establish that runtime is unchanged. The [layout benchmark](../method/namespace-benchmark.md) compares all seven original numerical outputs and the original input fields, validates the newly added solver label separately, then times four representative workloads against the pre-namespace checkout. It alternates timing order, isolates source paths and flags a configured slowdown threshold. MATLAB was unavailable in the cloud environment; this timing gate remains for local execution before merging.

`s_test_suite_namespace` checks public namespace calls, legacy adapters and measured waveform loading. The existing numerical suites now call the namespaced implementations directly. No trajectory baselines are regenerated. Markdown and LaTeX source checks remain separate from MATLAB numerical validation.
