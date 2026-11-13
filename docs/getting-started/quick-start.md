# Quick start

This example does three things. It converts physical inputs to IMR's dimensionless variables. It solves a Kelvin–Voigt bubble case. It plots the dimensional radius against time. Run this example in MATLAB from the checkout root.

```matlab
% Make the solver functions available.
imr_setup;
% Convert SI values and retain the reference scales.
[varin, ref] = imr.nondimensionalize('r0',50e-6,'req',5e-6,'mu',1e-2,'g',1e3, ...
    'tvector',linspace(0,12e-6,100));
% Use the pressure-form Keller–Miksis model and neo-Hookean Kelvin–Voigt stress.
[t, R] = imr.solve_nodal(varin{:}, 'radial', 2, 'stress', 1);
% Restore seconds and meters for the plot.
plot(t*ref.t0, R*ref.R0);
xlabel('Time (s)');
ylabel('Bubble radius (m)');
```

| Line or expression | Meaning |
|---|---|
| `imr_setup` | Loads this checkout's MATLAB functions |
| `r0`, `req` | Initial and equilibrium radii in meters |
| `mu`, `g` | Viscosity in Pa s and shear modulus in Pa |
| `tvector` | Requested times in seconds before conversion |
| `varin{:}` | Expands the converted name/value cell array |
| `radial = 2`, `stress = 1` | Selects the radial and stress equations; other options retain defaults |
| `ref.t0`, `ref.R0` | Restore dimensions to the solver outputs |

## Execution status

This example is the same example from the README. The authors checked it against `src/+imr/nondimensionalize.m` and the solver interface. Nobody ran MATLAB for this documentation change. This page does not show a captured transcript or a verified numerical output yet. You must run and check the example yourself. This page does not invent any sample numbers or plots.

After a successful run, check two things. Check that the returned histories cover the requested output times. Check that the histories contain only finite values. Use [outputs](../running/output.md) to interpret the arrays. Run the [test suites](../verification/tests.md) for stronger checks. See [examples](examples.md) for other workflows.
