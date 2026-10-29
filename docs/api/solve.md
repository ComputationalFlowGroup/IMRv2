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

`solver` accepts `'nodal'` or `'spectral'`, as a character vector or string scalar, case-insensitively. Later pairs override earlier values and a struct's `solver` field. `imr.default_inputs` stays solver-neutral so it can be passed to either existing wrapper. The returned `inputs.solver` records the selected solver, allowing replay through `imr.solve(inputs)`.

The solver is selected once before constructing operators and integrating. Both implementations retain their own field representations and [compatibility restrictions](../reference/compatibility.md). Unsupported options raise input errors; selection does not enable spectral mass transfer or automatically change to another solver. The existing `spectral` option continues to select a stress representation and does not choose the solver.

## Outputs

```matlab
[t,R,Rdot,P,T,Tm,kv,inputs,Rddot,S,Sdot] = imr.solve(p,'solver','nodal');
```

The first eight output positions are preserved. Optional output 9, `Rddot`, is the dimensionless wall acceleration at each returned time, evaluated from the governing equations. Optional outputs 10 and 11, `S` and `Sdot`, are the dimensionless stress integral and its time derivative, evaluated from the same right-hand-side call as `Rddot` at no added cost. No extra postprocessing runs when only the first eight outputs are requested. See [output shapes, SI scaling and acceleration cost](../running/output.md).

`imr.solve_nodal` and `imr.solve_spectral` remain supported and accept the same optional outputs 9-11. They fix the solver choice and reject a conflicting `solver` option. All calls share preparation, integration and optional reconstruction of `Rddot`, `S` and `Sdot` through `imr.solvers.run`; the numerical equations are implemented once per solver.

See [configuration](../running/configuration.md), [parameters](../model/parameters.md) and [solver modes](../running/solvers.md). These examples require local MATLAB validation; no numerical execution is claimed here.
