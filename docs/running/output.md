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
| 8 | `inputs` | Requested/default inputs plus actual time vector | A configuration record, not a physical output |

Radius, velocity and pressure histories use one row per returned time. Nodal thermal and concentration fields use one row per time and columns for spatial nodes. Spectral temperature fields instead use one row per collocation point and one column per time; this existing output orientation is preserved. The nodal bubble grid runs center to wall; the medium grid runs wall to infinity. `T(:,end)` and `Tm(:,1)` are the wall temperature when medium heat transfer is enabled; the wall closure recovers these values at each output time. The spectral solver reconstructs temperature values from its coefficients before returning them.

`imr.configuration.prepare_case` records inputs before applying all internal checker adjustments to that record. Use [compatibility](../reference/compatibility.md) to understand automatic stress-field selection, mode counts and the `Ca = Inf` stress-5 rule. The [API](../api/index.md) distinguishes input contracts from internal state.

## Printed diagnostics

`imr.solvers.describe_case` labels the solver, radial equation, constitutive law, viscosity, thermal/vapor/mass-transfer choices and stress representation. Its complete label rules live in the [display reference](../reference/functions/f_display.md).

`progdisplay = 1` prints on every right-hand-side evaluation and can slow runs. Collapse setup may print a recomputed equilibrium radius and a stress-history note. For stress 3–4 with medium heat transfer, the spectral solver prints the maximum difference between the local field's integral and `Z1/R^3`; this is a consistency diagnostic, not an automatic pass/fail accuracy certificate.

## Save and replay a run

In MATLAB from the checkout root:

```matlab
% Save a dimensionless configuration and its radius history.
imr_setup;
p = imr.default_inputs();
p.radial = 2;
p.stress = 1;
[t, R, ~, ~, ~, ~, ~, used] = imr.solve_nodal(p);
save('run1.mat', 't', 'R', 'used');
% Replay the input record through the same parser.
[t2, R2] = imr.solve_nodal(used);
```

This interface example has not been run here. Save `ref` too when your case starts with SI properties, and record the code revision and MATLAB version for reproducible comparisons. Keep generated files out of commits.
