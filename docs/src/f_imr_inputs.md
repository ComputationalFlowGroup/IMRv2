# `f_imr_inputs.m` Documentation

## Overview
`f_imr_inputs` returns every solver input with its default value, as a struct whose field names are the input names used in name/value pairs: `p.re` is the input `'re'`, `p.nt` is `'nt'`. It is the struct form of [`default_case`](default_case.md), which keeps the internal names (`Re8`, `Nt`, ...).

```matlab
p = f_imr_inputs();
p.radial = 2;
p.stress = 1;
p.re = 300;
[t, R] = f_imr_nodal(p);
```

## Using a struct
- **Same meaning as name/value pairs.** Only the fields that differ from the defaults are passed on, so `f_imr_nodal(p)` runs exactly like the name/value call with those inputs.
- **Struct, then pairs.** Name/value pairs after the struct override its fields: `f_imr_nodal(p, 're', 1000)`.
- **Optional inputs.** `p0star` and `stress0` have no default; add them as fields when needed (`p.stress0 = ...`).
- **Time.** Set either `p.tvector` or `p.tfin`, as with name/value pairs; a `tvector` that ends at `tfin` is taken as the time vector.
- **Unknown fields** raise the same `INPUT ERROR` as an unknown name/value pair.

## The inputs of a run
Both solvers return the struct of the inputs they used as an optional 8th output: the defaults, the given changes, and the time vector of the run.

```matlab
[t, R, ~, ~, ~, ~, ~, used] = f_imr_nodal(p);
save('run1.mat', 't', 'R', 'used');     % a record of the run
[t2, R2] = f_imr_nodal(used);           % reproduces it
```

## Sweeps
[`f_imr_sweep`](f_imr_sweep.md) accepts the struct as its shared `base`: `f_imr_sweep(@f_imr_nodal, p, 'ca', [10 100])`.
