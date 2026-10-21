# Parameter sweeps

`f_imr_sweep` runs a full grid of input combinations and records each result or error. Its [API page](../api/f_imr_sweep.md) defines options and array layout.

## Small dimensionless grid

In MATLAB from the checkout root:

```matlab
% Sweep Reynolds number and elasticity with shared dimensionless inputs.
imr_setup;
base = {'radial', 2, 'stress', 1, 'req', 0.2};
res = f_imr_sweep(@f_imr_nodal, base, 're', [100 300 1000], 'ca', [10 100]);
failed = arrayfun(@(r) ~isempty(r.error), res);
plot(res(2,1).out{1}, res(2,1).out{2});
```

Inspect `failed` and each `res(k).error` before using a result. A completed sweep can contain failed simulations.

## SI inputs and a parallel pool

```matlab
% Convert each physical case before running it on a parallel pool.
imr_setup;
base = {'r0', 50e-6, 'req', 10e-6, 'stress', 1, 'tvector', linspace(0, 20e-6, 200)};
res = f_imr_sweep(@f_imr_nodal, base, 'mu', logspace(-3, -1, 8), ...
    'g', logspace(2, 5, 8), 'transform', @f_dim2nondim, 'parallel', true);
failed = arrayfun(@(r) ~isempty(r.error), res);
```

A transform requires a cell array for `base`; a struct cannot be combined with a transform. Parallel execution requires Parallel Computing Toolbox. Omit `parallel` or use `false` for serial execution. These examples are source-checked and await local MATLAB execution.

## Large bundled studies

`examples/f_multirun.m` returns five studies: Newtonian, neo-Hookean, Kelvin–Voigt, UCM and Oldroyd-B. Most contain 4096 combinations; the UCM study has `ceil(4096^(1/3))^3 = 4913`. The relaxation times are built as factors of `mu/G`, then converted to dimensionless inputs. These large grids are not installation smoke tests.

`examples/matlab-submitjob.sh` is a site-specific Slurm template that calls `matlab-threaded` and saves `multirun.mat`. Adapt scheduler resources and MATLAB invocation to your cluster before using it. No cluster execution was validated here.
