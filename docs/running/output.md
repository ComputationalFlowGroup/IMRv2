# Outputs

Both solvers return dimensionless arrays and print a model summary through `imr.solvers.describe_case`. They do not automatically save a simulation file.

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

Radius, velocity and pressure histories use one row per returned time. Nodal thermal and concentration fields use one row per time and columns for spatial nodes. Spectral temperature fields instead use one row per collocation point and one column per time; this existing output orientation is preserved. The nodal bubble grid runs center to wall; the medium grid runs wall to infinity. `T(:,end)` and `Tm(:,1)` are the wall temperature when medium heat transfer is enabled; the wall closure recovers these values at each output time. The spectral solver reconstructs temperature values from its coefficients before returning them.

`imr.configuration.prepare_case` records inputs before applying all internal checker adjustments to that record. Use [compatibility](../reference/compatibility.md) to understand automatic stress-field selection, mode counts and the `Ca = Inf` stress-9 rule. The [API](../api/index.md) distinguishes input contracts from internal state.

## Wall acceleration, stress integral and its derivative

Request output 9, 10 or 11 from `imr.solve`, `imr.solve_nodal` or `imr.solve_spectral`:

```matlab
[t,R,Rdot,P,T,Tm,kv,inputs,Rddot,S,Sdot] = imr.solve('radial',2,'req',1);
```

The shared runner evaluates the coupled right-hand side once at every returned time and state and extracts the derivative of wall velocity (`Rddot`) and, since the right-hand side already computes them to evaluate the radial equation, the stress integral (`S`) and its time derivative (`Sdot`) at no extra cost beyond that one evaluation. This includes the selected radial equation, forcing at the actual time, stress and thermal/vapor wall closures. It does not estimate any of these by differencing a sampled history, add an ODE state or change the time integration. Postprocessing uses the same boundary reconstruction as the RHS and suppresses its progress lines.

Calls requesting at most eight outputs perform no extra RHS evaluations. Requesting output 9, 10 or 11 adds one RHS evaluation per returned time (shared across all three, so asking for 10 or 11 costs the same as asking for 9 alone), which can be significant for large field grids. `Rddot`, `S` and `Sdot` follow the returned adaptive steps for a two-time span, or the requested output times for a longer span. Discontinuous forcing uses the governing equations' pointwise convention at the discontinuity.

The dimensional acceleration is

$$
\ddot R_{\mathrm{SI}}=\frac{U_c^2}{R_0}\ddot R
=\frac{R_0}{t_0^2}\ddot R.
$$

It has units of meters per second squared. The first eight outputs retain their positions and field orientations.

## Printed diagnostics

`imr.solvers.describe_case` labels the solver, radial equation, constitutive law, viscosity, thermal/vapor/mass-transfer choices and stress representation. Its complete label rules live in the [display reference](../reference/functions/f_display.md).

`progdisplay = 1` prints on every right-hand-side evaluation and can slow runs. Collapse setup may print a recomputed equilibrium radius and a stress-history note. For stress 5–8 with medium heat transfer, the spectral solver prints the maximum difference between the local field's integral and `Z1/R^3`; this is a consistency diagnostic, not an automatic pass/fail accuracy certificate.

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

This interface example has not been run here. Save `ref` too when your case starts with SI properties, and record the code revision and MATLAB version for reproducible comparisons. Keep generated files out of commits.
