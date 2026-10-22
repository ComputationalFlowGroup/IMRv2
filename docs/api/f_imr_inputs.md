# `imr.default_inputs`

`imr.default_inputs` returns every solver input with its default value, as a struct whose field names are the input names used in name/value pairs: `p.re` is the input `'re'`, `p.nt` is `'nt'`. It is the struct form of [`imr.configuration.defaults`](../model/parameters.md), which keeps the internal names (`Re8`, `Nt`, ...).

```matlab
imr_setup;
p = imr.default_inputs();
p.radial = 2;
p.stress = 1;
p.re = 300;
[t, R] = imr.solve_nodal(p);
```

## Using a struct

- **Same meaning as name/value pairs.** Only the fields that differ from the defaults are passed on, so `imr.solve_nodal(p)` runs exactly like the name/value call with those inputs.
- **Struct, then pairs.** Name/value pairs after the struct override its fields: `imr.solve_nodal(p, 're', 1000)`.
- **Optional inputs.** `p0star` and `stress0` have no default; add them as fields when needed (`p.stress0 = ...`).
- **Time.** Set either `p.tvector` or `p.tfin`, as with name/value pairs; a `tvector` that ends at `tfin` is taken as the time vector.
- **Unknown fields** raise the same `INPUT ERROR` as an unknown name/value pair.

## The inputs of a run

The optional eighth solver output records requested/default inputs and the time vector. See [save and replay](../running/output.md#save-and-replay-a-run) for a complete example and the distinction between recorded inputs and internal checker adjustments.

## Sweeps

[`imr.sweep`](f_imr_sweep.md) accepts the struct as its shared `base`: `imr.sweep(@imr.solve_nodal, p, 'ca', [10 100])`.
