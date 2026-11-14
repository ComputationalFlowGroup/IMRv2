# `imr.default_inputs`

`imr.default_inputs` returns every solver input with its default value, as a struct. Each field name is the input name used in a name/value pair: `p.re` is the input `'re'`, and `p.nt` is `'nt'`. This struct is the struct form of [`imr.configuration.defaults`](../model/parameters.md). `imr.configuration.defaults` keeps the internal names (`Re8`, `Nt`, ...).

```matlab
imr_setup;
p = imr.default_inputs();
p.radial = 2;
p.stress = 1;
p.re = 300;
[t, R] = imr.solve_nodal(p);
```

## Using a struct

- **Same meaning as name/value pairs.** The function passes on only the fields that differ from the defaults. `imr.solve_nodal(p)` then runs exactly like the matching name/value call.
- **Struct, then pairs.** Name/value pairs after the struct override its fields: `imr.solve_nodal(p, 're', 1000)`.
- **Optional inputs.** `p0star` and `stress0` have no default. Add them as fields when you need them (`p.stress0 = ...`).
- **Time.** Set either `p.tvector` or `p.tfin`, the same as with name/value pairs. If `tvector` ends at `tfin`, the solver uses `tvector` as the time vector.
- **Unknown fields** raise the same `INPUT ERROR` as an unknown name/value pair.

## The inputs of a run

The optional eighth solver output records requested/default inputs and the time vector. See [save and replay](../running/output.md#save-and-replay-a-run) for a complete example. That page also explains the difference between recorded inputs and internal checker adjustments.

`imr.default_inputs` stays solver-neutral. Add `p.solver = 'spectral'` or pass `'solver','spectral'` to [`imr.solve`](solve.md) when needed. The returned record includes the selected solver. This lets `imr.solve(used)` replay either implementation.

## Sweeps

[`imr.sweep`](sweep.md) accepts the struct as its shared `base`: `imr.sweep(@imr.solve_nodal, p, 'ca', [10 100])`.
