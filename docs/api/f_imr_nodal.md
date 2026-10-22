# `f_imr_nodal`

`f_imr_nodal` integrates radial bubble dynamics with optional heat transfer, vapor mass transfer and viscoelastic stress. Fields are carried at nodes; `scheme` chooses the differentiation matrices.

```matlab
[t, R, Rdot, P, T, Tm, kv, inputs] = f_imr_nodal(varargin);
```

`varargin` accepts dimensionless name/value pairs or a struct from [`f_imr_inputs`](f_imr_inputs.md), optionally followed by overriding pairs. See [configuration](../running/configuration.md) and the canonical [parameter table](../model/parameters.md).

| Choice | Behavior |
|---|---|
| `scheme = 0` | Default second-order finite differences |
| `scheme = 1` | Chebyshev differentiation matrices with the same physical equations |
| `jacobian = 2` | Default hybrid analytic Jacobian when bubble heat is enabled; otherwise sparsity-pattern fallback |
| `jacobian = 3` | Batched finite-difference Jacobian for all model configurations accepted by the checker |

The [output contract](../running/output.md) defines every returned array. [Examples](../getting-started/examples.md) show invocation forms. [Compatibility](../reference/compatibility.md) defines restrictions; in particular the nodal two-state stress-5 form cannot supply medium heating.

Implementation: `f_prepare_case` → `f_imr_system` → `f_odesolve` → output reconstruction. Nodal operators and closures live in `src/private/imr_nodal_system.m`. Read the [numerical workflow](../method/overview.md), [wall closure](../method/wall-boundaries.md), and [integration method](../method/time-integration.md). Historical comparisons to the old `f_imr_fd` and convergence measurements live in [validation](../validation/convergence.md).
