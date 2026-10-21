# Quick start

This example converts physical inputs to IMR's dimensionless variables, solves a Kelvin–Voigt bubble case, and plots dimensional radius against time. Run it in MATLAB from the checkout root.

```matlab
% Make the solver functions available.
addpath('src');
% Convert SI values and retain the reference scales.
[varin, ref] = f_dim2nondim('r0',50e-6,'req',5e-6,'mu',1e-2,'g',1e3, ...
    'tvector',linspace(0,12e-6,100));
% Use the pressure-form Keller–Miksis model and neo-Hookean Kelvin–Voigt stress.
[t, R] = f_imr_nodal(varin{:}, 'radial', 2, 'stress', 1);
% Restore seconds and meters for the plot.
plot(t*ref.t0, R*ref.R0);
xlabel('Time (s)');
ylabel('Bubble radius (m)');
```

| Line or expression | Meaning |
|---|---|
| `addpath('src')` | Loads this checkout's MATLAB functions |
| `r0`, `req` | Initial and equilibrium radii in meters |
| `mu`, `g` | Viscosity in Pa s and shear modulus in Pa |
| `tvector` | Requested times in seconds before conversion |
| `varin{:}` | Expands the converted name/value cell array |
| `radial = 2`, `stress = 1` | Selects the radial and stress equations; other options retain defaults |
| `ref.t0`, `ref.R0` | Restore dimensions to the solver outputs |

## Execution status

This is the existing README example, checked against `src/f_dim2nondim.m` and the solver interface. MATLAB was not run for this documentation change. There is no captured transcript or verified numerical output to display yet; the user will test locally. No sample numbers or plots have been invented.

After a successful run, inspect that the returned histories cover the requested output times and contain finite values. Use [outputs](../running/output.md) to interpret the arrays, and run the [test suites](../validation/tests.md) for stronger checks. See [examples](examples.md) for other workflows.
