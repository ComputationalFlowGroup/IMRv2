# `imr.solve`

`imr.solve` is the common entry point for nodal and spectral simulations. It accepts dimensionless name/value pairs or an input struct followed by overriding pairs.

```matlab
imr_setup;
p = imr.default_inputs();
p.radial = 2;
p.req = 1;
p.rdotzero = 0.01;
[t,R] = imr.solve(p);                     % nodal, the default
[ts,Rs] = imr.solve(p,'solver','spectral');
```

`solver` accepts `'nodal'` or `'spectral'`, as a character vector or string scalar, case-insensitively. Later pairs override earlier values and a struct's `solver` field. `imr.default_inputs` stays solver-neutral. You can pass it to either existing wrapper. The returned `inputs.solver` records the selected solver. This lets you replay through `imr.solve(inputs)`.

The code selects the solver once, before it constructs operators and integrates. Both implementations retain their own field representations and [compatibility restrictions](compatibility.md). Unsupported options raise input errors. Selecting a solver does not enable spectral mass transfer. Selecting a solver does not automatically change to another solver. The existing `spectral` option continues to select a stress representation. It does not choose the solver.

## Outputs

```matlab
[t,R,Rdot,P,T,Tm,kv,inputs,Rddot,S,Sdot] = imr.solve(p,'solver','nodal');
```

The code preserves the first eight output positions. Optional output 9, `Rddot`, is the dimensionless wall acceleration at each returned time, evaluated from the governing equations. Optional outputs 10 and 11, `S` and `Sdot`, are the dimensionless stress integral and its time derivative, evaluated from the same right-hand-side call as `Rddot` at no added cost. If you request only the first eight outputs, no extra postprocessing runs. See [output shapes, SI scaling and acceleration cost](../running/output.md).

`imr.solve_nodal` and `imr.solve_spectral` remain supported and accept the same optional outputs 9-11. They fix the solver choice and reject a conflicting `solver` option. All calls share preparation, integration and optional reconstruction of `Rddot`, `S` and `Sdot`, through `imr.solvers.run`. Each solver implements its numerical equations once.

See [configuration](../running/configuration.md), [parameters](../model/parameters.md) and [solver modes](../running/solvers.md). These examples require local MATLAB validation. This page claims no numerical execution.
