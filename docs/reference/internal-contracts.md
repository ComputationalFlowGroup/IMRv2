# Prepared cases and numerical systems

The common public entry point is [`imr.solve`](solve.md). `imr.solve_nodal` and `imr.solve_spectral` remain supported wrappers. Internally, each one prepares a named case, constructs a numerical system, integrates it, then reconstructs the original eight outputs and optional wall acceleration. These developer contracts replace positional option vectors inside the solvers. They are not additional simulation options.

```matlab
c = imr.configuration.prepare_case('nodal','radial',2,'Req',1,'Nt',12);
system = imr.solvers.build(c);
n = c.numerics;
[t,X] = imr.numerics.integrate(system.rhs,system.initial,n.method,n.divisions, ...
    c.time,system.jacobian,n.rtol,n.atol);
out = system.outputs(t,X);
R = out{2};
```

## Named case

| Field | Responsibility |
|---|---|
| `solver` | `nodal` or `spectral` |
| `model` | Radial, thermal, vapor, mass-transfer and stress selections |
| `numerics` | Integrator, spatial orders, tolerances and Jacobian selection |
| `initial` | Radius, velocity, total/gas/vapor pressure, equilibrium radius and stress |
| `acoustic` | Ambient sound speed and pressure-EoS coefficients |
| `forcing` | Analytic waveform settings or measured spline data |
| `stress`, `thermal`, `mass` | Named dimensionless coefficients and derived groups |
| `time`, `output` | Requested times and progress display |
| `inputs` | Public requested/default record, preserving the existing replay semantics |
| `layout` | Effective state blocks after compatibility normalization |

`imr.configuration.defaults` owns default values and derives physical groups through the SI converter. `imr.configuration.schema` owns public/internal names and unconditional choice/count rules. `imr.configuration.validate_case` owns physical ranges and coupled compatibility constraints. `imr.configuration.prepare_case` applies overrides, computes derived groups and initializes collapse when requested. The [parameter table](../model/parameters.md) remains the user-facing definition of each input.

`prepare_case('',...)` selects the solver from the routing option. It defaults to nodal. A fixed first argument rejects a conflicting public selector. The code parses the selector separately from the numerical/physical schema. This keeps `imr.default_inputs` solver-neutral. The normalized selection is recorded as `inputs.solver` for replay.

For example, `inputs.nv` records the requested count. `numerics.Nv` records the count selected for the active stress model. This distinction preserves the existing [input-record limits](../running/output.md).

## State and shape contract

`imr.solvers.state_layout` assigns disjoint row-vector index blocks, in this order: `R`, `Rdot`, `P`, `bubble`, `medium`, `vapor`, `stress1`, `stress2`, `polymer`. Empty blocks represent disabled fields. `count` is the total number of states. The code never uses a negative grid order as an inactivity marker.

| Quantity | Nodal | Spectral |
|---|---|---|
| State | Column, or one state per column in an RHS batch | Column |
| Bubble field | `Nt` Kirchhoff values | `Nt+1` even Chebyshev coefficients |
| Medium field | `Mt` temperature values | `Mt+1` Chebyshev coefficients |
| Vapor field | `Nt` mass fractions when active | Unsupported |
| Stored trajectory | One time per row | One time per row |
| Returned temperature | Time by node | Collocation point by time |

Stress blocks depend on the constitutive model; the [function call map](function-call-map.md#3-integration-one-checked-case-becomes-histories) describes them. Initialization, RHS, Jacobians and output reconstruction all use the same layout object.

## Numerical system

`imr.solvers.build` chooses `imr.solvers.build_nodal` or `imr.solvers.build_spectral`. Both implementations live in `src/+imr/+solvers/`. MATLAB resolves their fully qualified names through the parent `src` directory. There is no `private` directory. There is no recursive path setup. Package functions remain callable for developer tests. See the [namespace layout](namespaces.md).

| Member | Contract |
|---|---|
| `initial` | Column of length `layout.count` |
| `rhs(t,X)` | State derivative with matching shape; nodal supports batches |
| `output_rhs(t,X)` | `[dXdt,S,Sdot]`, derivatives plus the stress integral and its time derivative, with progress display suppressed; used for requested acceleration, `S` and `Sdot` |
| `reset_output_rhs(X)` | Optional nodal wall-cache reset at the first state of an output history |
| `jacobian` | Selected derivative handle, sparsity matrix or empty fallback |
| `outputs(t,X)` | Cell row containing the first eight public outputs |
| `operatorSizes` | Bubble and medium operator row counts; zero when disabled |
| `layout` | The prepared case's state indices |

Factories build only enabled operators. They retain local wall closures and scheme-specific discretization. `imr.physics.thermal_inverse` and `imr.physics.saturation_pressure` share pure expressions with the Jacobians and property evaluator. System construction does not integrate time. Output reconstruction recovers algebraic wall states, using the same closures as the RHS.

When you request output 9, 10 or 11, the shared runner does four things. It resets the nodal wall cache to the same starting guess used by field reconstruction. It evaluates `output_rhs` once per returned state. It extracts index `layout.Rdot` for `Rddot`. It reads `S`/`Sdot` straight from `output_rhs`'s other two outputs. One evaluation serves all three outputs, so requesting 10 or 11 costs the same as requesting 9 alone. The runner appends this history after the eight entries from `outputs`. Other calls perform no acceleration/stress postprocessing. Other calls preserve the numerical callbacks used for integration.

## Verification

`tests/suites/s_test_suite_structure.m` checks public-name/default consistency, struct versus pair overrides, legacy serialization, state activation and independent centered differences of the complete RHS against each available Jacobian. It covers nodal schemes, thermal/vapor combinations, constitutive laws and spectral stress fields. It also checks reconstructed wall temperature/saturation, the common polytropic ODE limit for all six radial models, and a small thermal nodal/spectral comparison.

The suite also rejects a zero-mode active spectral stress field before construction. If it did not, the collocation matrix would divide by zero. Disabled fields can still have zero grid counts.

These checks require MATLAB. This page does not claim that anyone ran a numerical pass in the cloud environment. The existing [test workflow](../verification/tests.md) includes this suite. The suite regenerates no trajectory reference. `toolchain/s_benchmark_structure.m` reports local construction, RHS and Jacobian timings for profiling. Compare runs on the same MATLAB release and the same machine before you draw performance conclusions.
