# Parameter sweeps

`imr.sweep` runs a full grid of input combinations. It records each result or error. Its [API page](../reference/sweep.md) defines the options and the array layout.

## Small dimensionless grid

In MATLAB from the checkout root:

```matlab
% Sweep Reynolds number and elasticity with shared dimensionless inputs.
imr_setup;
base = {'radial', 2, 'stress', 1, 'req', 0.2};
res = imr.sweep(@imr.solve_nodal, base, 're', [100 300 1000], 'ca', [10 100]);
failed = arrayfun(@(r) ~isempty(r.error), res);
plot(res(2,1).out{1}, res(2,1).out{2});
```

Inspect `failed` and each `res(k).error` before using a result. A completed sweep can contain failed simulations.

## SI inputs and a parallel pool

```matlab
% Convert each physical case before running it on a parallel pool.
imr_setup;
base = {'r0', 50e-6, 'req', 10e-6, 'stress', 1, 'tvector', linspace(0, 20e-6, 200)};
res = imr.sweep(@imr.solve_nodal, base, 'mu', logspace(-3, -1, 8), ...
    'g', logspace(2, 5, 8), 'transform', @imr.nondimensionalize, 'parallel', true);
failed = arrayfun(@(r) ~isempty(r.error), res);
```

A transform requires a cell array for `base`. You cannot combine a struct with a transform. Parallel execution requires the Parallel Computing Toolbox. Omit `parallel`, or set it to `false`, for serial execution. The authors checked these examples against the source. Nobody has run these examples on local MATLAB yet.

## Large bundled studies

`examples/f_multirun.m` returns five studies: Newtonian, neo-Hookean, Kelvin–Voigt, UCM and Oldroyd-B. Most studies contain 4096 combinations. The UCM study has `ceil(4096^(1/3))^3 = 4913` combinations. The script builds the relaxation times as factors of `mu/G`. The script then converts them to dimensionless inputs. These large grids are not installation smoke tests.

`examples/matlab-submitjob.sh` is a site-specific Slurm template. It calls `matlab-threaded`. It saves `multirun.mat`. Before you use it, adapt the scheduler resources and the MATLAB invocation to your cluster. Nobody validated a cluster execution here.
