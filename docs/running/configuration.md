# Configuration

`f_call_params` builds a run from `default_case`, optional input-struct changes and explicit name/value pairs. Both solvers accept dimensionless values only.

## Forms and precedence

| Form | Example | Behavior |
|---|---|---|
| Name/value pairs | `f_imr_nodal('radial',2,'stress',1)` | Overrides defaults in order; later repeated names win |
| Cell array | `args = {'radial',2}; f_imr_nodal(args{:})` | Expands to name/value pairs |
| Struct | `p = f_imr_inputs(); p.radial = 2; f_imr_nodal(p)` | Passes only fields differing from the defaults |
| Struct plus pairs | `f_imr_nodal(p,'radial',1)` | Explicit pairs override changed struct fields |

Names in the name/value parser are case-insensitive. Use the lowercase fields returned by [`f_imr_inputs`](../api/f_imr_inputs.md) for structs. Optional `p0star` and `stress0` fields have no default and can be added when needed.

The full accepted names, meanings and defaults live in [parameters](../model/parameters.md). Internal variables such as `Rzero`, `Pb_star`, `Szero`, `tspan`, `wave_poly` and `wave_dpoly` are not accepted public option names.

## Time inputs

- Give `tvector` or `tfin`, not both in explicit pairs.
- `tvector` must contain at least two finite, strictly increasing values and start at a nonnegative time.
- `tfin` alone constructs `[0, tfin]`, returning the integrator's own output steps.
- A struct whose changed `tvector` ends at its changed `tfin` passes only `tvector`; the parser removes that redundant `tfin`.

The optional eighth solver output records the requested/default inputs, with its final time vector. `f_checkers` can adjust internal settings such as `Nv`, `Ca` and `spectral`; the recorded struct is not a complete dump of those internal adjustments. See [outputs](output.md).

## Units and case files

Use [`f_dim2nondim`](../api/f_dim2nondim.md) for physical SI properties. Model flags and unrecognized conversion names pass through that helper to the solver parser.

There is no case-file loader or config-file format. The old `casefile` option raises `INPUT ERROR: casefile is no longer supported.` Keep cases as cell arrays or structs. [Troubleshooting](troubleshooting.md) explains input errors and warnings.
