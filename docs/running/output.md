# Outputs

Both solvers return dimensionless arrays. Both solvers print a model summary through `imr.solvers.describe_case`. Neither solver saves a simulation file automatically.

| Position | Name | Meaning | Restore SI units |
|---|---|---|---|
| 1 | `t` | Time vector | `t*ref.t0` |
| 2 | `R` | Bubble radius | `R*ref.R0` |
| 3 | `Rdot` | Wall velocity | `Rdot*ref.Uc` |
| 4 | `P` | Bubble pressure | `P*ref.P8` |
| 5 | `T` | Bubble temperature; empty without bubble heat transfer | `T*ref.T8` |
| 6 | `Tm` | Medium temperature; a column of ones without medium heat transfer | `Tm*ref.T8` |
| 7 | `kv` | Vapor mass fraction; empty without mass transfer and always empty for spectral solver | Already dimensionless |
| 8 | `inputs` | Requested/default inputs plus actual time vector and selected solver | A configuration record, not a physical output |
| 9 | `Rddot` | Optional bubble-wall acceleration, one column aligned with `t` | `Rddot*ref.Uc^2/ref.R0` |
| 10 | `S` | Optional stress integral, one column aligned with `t`; identically zero for `stress = 0` | `S*ref.P8` |
| 11 | `Sdot` | Optional time derivative of `S`, one column aligned with `t` | `Sdot*ref.P8/ref.t0` |

Radius, velocity and pressure histories use one row per returned time. Nodal thermal and concentration fields use one row per time, with columns for the spatial nodes. Spectral temperature fields use a different orientation: one row per collocation point, and one column per time. This difference is deliberate; the documentation keeps this existing output orientation. The nodal bubble grid runs center to wall. The nodal medium grid runs wall to infinity. `T(:,end)` and `Tm(:,1)` give the wall temperature when medium heat transfer is enabled. The wall closure recovers these values at each output time. The spectral solver reconstructs temperature values from its coefficients before it returns them.

`imr.configuration.prepare_case` records inputs before applying all internal checker adjustments to that record. Use [compatibility](../reference/compatibility.md) to understand automatic stress-field selection, mode counts and the `Ca = Inf` stress-9 rule. The [API](../reference/api.md) distinguishes input contracts from internal state.

## Wall acceleration, stress integral and its derivative

Request output 9, 10 or 11 from `imr.solve`, `imr.solve_nodal` or `imr.solve_spectral`:

```matlab
[t,R,Rdot,P,T,Tm,kv,inputs,Rddot,S,Sdot] = imr.solve('radial',2,'req',1);
```

The shared runner evaluates the coupled right-hand side once, at every returned time and state. This one evaluation extracts the derivative of wall velocity (`Rddot`). The right-hand side already computes the stress integral (`S`) and its time derivative (`Sdot`) to evaluate the radial equation. So this one evaluation extracts `S` and `Sdot` too, at no extra cost. This evaluation includes the selected radial equation, the forcing at the actual time, and the stress and thermal/vapor wall closures. The runner does not estimate any of these values by differencing a sampled history. The runner does not add an ODE state. The runner does not change the time integration. Postprocessing uses the same boundary reconstruction as the right-hand side. Postprocessing suppresses the right-hand side's progress lines.

A call that requests at most eight outputs performs no extra RHS evaluations. A call that requests output 9, 10 or 11 adds one RHS evaluation per returned time. This one extra evaluation is shared across all three outputs. So a call for output 10 or 11 costs the same as a call for output 9 alone. This extra cost can be significant for large field grids. For a two-time span, `Rddot`, `S` and `Sdot` follow the returned adaptive steps. For a longer span, they follow the requested output times. Discontinuous forcing uses the governing equations' pointwise convention at the discontinuity.

The dimensional acceleration is

$$
\ddot R_{\mathrm{SI}}=\frac{U_c^2}{R_0}\ddot R
=\frac{R_0}{t_0^2}\ddot R.
$$

It has units of meters per second squared. The first eight outputs retain their positions and field orientations.

## Printed diagnostics

`imr.solvers.describe_case` labels the solver, radial equation, constitutive law, viscosity, thermal/vapor/mass-transfer choices and stress representation. Its complete label rules live in the [display reference](../reference/functions/describe_case.md).

`progdisplay = 1` prints on every right-hand-side evaluation. This can slow runs. Collapse setup can print a recomputed equilibrium radius and a stress-history note. For stress 5–8 with medium heat transfer, the spectral solver prints the maximum difference between the local field's integral and `Z1/R^3`. This printed value is a consistency diagnostic. It is not an automatic pass/fail accuracy certificate.

## Save and replay a run

In MATLAB from the checkout root:

```matlab
% Save a dimensionless configuration and its radius history.
imr_setup;
p = imr.default_inputs();
p.radial = 2;
p.stress = 1;
[t, R, ~, ~, ~, ~, ~, used] = imr.solve(p,'solver','nodal');
save('run1.mat', 't', 'R', 'used');
% Replay the input record through the same parser.
[t2, R2] = imr.solve(used);
```

Nobody ran this interface example here. If your case starts with SI properties, save `ref` too. For reproducible comparisons, record the code revision and the MATLAB version. Keep generated files out of commits.
