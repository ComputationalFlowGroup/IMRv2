# `f_imr_sweep`

`f_imr_sweep` runs a solver over every combination of the given input values (a full grid) and collects the results. Each run starts from the defaults of [`default_case`](../model/parameters.md), then the inputs shared by all runs (`base`), then the swept values.

```text
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
