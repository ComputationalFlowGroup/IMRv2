# `f_imr_sweep.m` Documentation

## Overview
`f_imr_sweep` runs a solver over every combination of the given input values (a full grid) and collects the results. Each run starts from the defaults of [`default_case`](default_case.md), then the inputs shared by all runs (`base`), then the swept values.

```matlab
[res, grid] = f_imr_sweep(solver, base, name1, values1, name2, values2, ..., option, value, ...)
```

## Inputs

| Name | Description |
|------|-------------|
| `solver` | Function handle, `@f_imr_nodal` or `@f_imr_spectral` |
| `base` | Inputs shared by every run: a cell array of name/value pairs (`{}` for the defaults) or a struct from [`f_imr_inputs`](f_imr_inputs.md). A transform needs the cell array form. |
| `name_k`, `values_k` | A solver input and the values it takes, a numeric vector or a cell array. Every combination is run. |

Options, given after the swept inputs:

| Option | Default | Description |
|--------|---------|-------------|
| `'transform'` | none | Function handle applied to each run's name/value pairs before the solver, for example `@f_dim2nondim` to sweep SI values |
| `'outputs'` | 2 | Number of solver outputs to keep (`t`, `R`, `Rdot`, `P`, `T`, `Tm`, `kv`) |
| `'parallel'` | `false` | `true` runs the grid with `parfor` on a parallel pool (Parallel Computing Toolbox); `false` runs it in series in base MATLAB |
| `'quiet'` | `true` | Hides the solver's console output |

## Outputs

| Name | Description |
|------|-------------|
| `res` | Struct array of size `[numel(values1), numel(values2), ...]`; the first swept input varies fastest. Fields: `inputs` (the run's swept name/value pairs), `out` (cell of the kept solver outputs), `error` (the error message of a failed run, `''` otherwise) |
| `grid` | Struct with the swept `names` and their `values` |

A failed run records its error message and the sweep continues.

## Examples

```matlab
% Reynolds number and elasticity, dimensionless inputs
base = {'radial', 2, 'stress', 1, 'req', 0.2};
res = f_imr_sweep(@f_imr_nodal, base, 're', [100 300 1000], 'ca', [10 100]);
plot(res(2,1).out{1}, res(2,1).out{2});      % Re = 300, Ca = 10

% viscosity and shear modulus in SI units, on a parallel pool
base = {'r0', 50e-6, 'req', 10e-6, 'stress', 1, 'tvector', linspace(0, 20e-6, 200)};
res = f_imr_sweep(@f_imr_nodal, base, 'mu', logspace(-3, -1, 8), ...
    'g', logspace(2, 5, 8), 'transform', @f_dim2nondim, 'parallel', true);
failed = arrayfun(@(r) ~isempty(r.error), res);
```

`examples/f_multirun.m` runs a larger study with this function: Newtonian, Kelvin–Voigt, UCM and Oldroyd-B grids of 4096 runs each, from SI values.
