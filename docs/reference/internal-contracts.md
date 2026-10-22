# Prepared cases and numerical systems

The public entry points remain `f_imr_nodal` and `f_imr_spectral`. Internally they prepare a named case, construct a numerical system, integrate it, then reconstruct the eight public outputs. These developer contracts replace positional option vectors inside the solvers; they are not additional simulation options.

```matlab
c = f_prepare_case('nodal','radial',2,'Req',1,'Nt',12);
system = f_imr_system(c);
n = c.numerics;
[t,X] = f_odesolve(system.rhs,system.initial,n.method,n.divisions, ...
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

`default_case` owns default values and derives physical groups through the SI converter. `f_parameter_schema` owns public/internal names and unconditional choice/count rules. `f_checkers` owns physical ranges and coupled compatibility constraints. `f_prepare_case` applies overrides, computes derived groups and initializes collapse when requested. The [parameter table](../model/parameters.md) remains the user-facing definition of each input.

For example, `inputs.nv` records the requested count while `numerics.Nv` records the count selected for the active stress model. This distinction preserves the existing [input-record limits](../running/output.md). `f_call_params` remains a compatibility adapter that serializes the named case into its original twelve outputs; new solver code should use named fields.

## State and shape contract

`f_state_layout` assigns disjoint row-vector index blocks in the order `R`, `Rdot`, `P`, `bubble`, `medium`, `vapor`, `stress1`, `stress2`, `polymer`. Empty blocks represent disabled fields. `count` is the total number of states. No negative grid order is used as an inactivity marker.

| Quantity | Nodal | Spectral |
|---|---|---|
| State | Column, or one state per column in an RHS batch | Column |
| Bubble field | `Nt` Kirchhoff values | `Nt+1` even Chebyshev coefficients |
| Medium field | `Mt` temperature values | `Mt+1` Chebyshev coefficients |
| Vapor field | `Nt` mass fractions when active | Unsupported |
| Stored trajectory | One time per row | One time per row |
| Returned temperature | Time by node | Collocation point by time |

Stress blocks depend on the constitutive model; the [code map](code-map.md#3-integration-one-checked-case-becomes-histories) describes them. Initialization, RHS, Jacobians and output reconstruction all use the same layout object.

## Numerical system

`f_imr_system` chooses `src/private/imr_nodal_system.m` or `imr_spectral_system.m`. MATLAB resolves this `private` directory through the parent `src` function; users need only `imr_setup`, with no extra path arguments.

| Member | Contract |
|---|---|
| `initial` | Column of length `layout.count` |
| `rhs(t,X)` | State derivative with matching shape; nodal supports batches |
| `jacobian` | Selected derivative handle, sparsity matrix or empty fallback |
| `outputs(t,X)` | Cell row containing the eight unchanged public outputs |
| `operatorSizes` | Bubble and medium operator row counts; zero when disabled |
| `layout` | The prepared case's state indices |

Factories build only enabled operators. They retain local wall closures and scheme-specific discretization, while `f_thermal_inverse` and `f_saturation_pressure` share pure expressions with the Jacobians and property evaluator. System construction does not integrate time. Output reconstruction recovers algebraic wall states using the same closures as the RHS.

## Verification

`toolchain/s_test_suite_structure.m` checks public-name/default consistency, struct versus pair overrides, legacy serialization, state activation and independent centered differences of the complete RHS against each available Jacobian. It covers nodal schemes, thermal/vapor combinations, constitutive laws and spectral stress fields. It also checks reconstructed wall temperature/saturation, the common polytropic ODE limit for all six radial models, and a small thermal nodal/spectral comparison.

The suite also rejects a zero-mode active spectral stress field before construction; otherwise its collocation matrix would divide by zero. Disabled fields may still have zero grid counts.

These checks require MATLAB and are not a claim of an executed numerical pass in the cloud environment. The existing [test workflow](../validation/tests.md) includes this suite. No trajectory reference is regenerated. `toolchain/s_benchmark_structure.m` reports local construction, RHS and Jacobian timings for profiling; compare runs on the same MATLAB release and machine before drawing performance conclusions.
