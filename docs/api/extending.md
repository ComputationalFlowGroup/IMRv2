# Extending IMR

IMR is extended by editing its MATLAB functions directly. There is no plugin loader; each input must pass from defaults through parsing and validation to the equations that use it.

For a new fluid using an existing EoS family, start with a [MATLAB fluid profile](../model/fluids.md). Material coefficients and saturation data belong there; the shared evaluator and SI converter handle selection. The trace below applies when adding a numerical option or a mathematical model that is not already supported.

## Worked trace: the existing `scheme` input

| Step | Existing implementation | What an extension must preserve |
|---|---|---|
| Default | `src/default_case.m`: `scheme = 0`, then `d.scheme = scheme` | One authoritative default |
| Public struct | `src/f_imr_inputs.m` copies unchanged field names | Public field matches parser name, or belongs in the `renamed` table |
| Parser | `src/f_call_params.m` reads `'scheme'` | Explicit accepted name and last-pair precedence |
| Validation | `src/f_checkers.m` accepts 0/1 and rejects nonzero scheme for spectral solver | Reject incompatible configurations before integration |
| Transport | `scheme` is appended to `solve_opts` | Preserve every existing positional index |
| Use | `src/f_imr_nodal.m` passes it to `f_diff_mat` | The input actually changes the intended operator |
| Evidence | `toolchain/s_test_suite_diffmat.m` and `f_test_cases.m` exercise modes | Meaningful tests for the new behavior |
| Documentation | [Parameters](../model/parameters.md), [discretization](../method/discretization.md), [compatibility](../reference/compatibility.md) | One canonical definition, linked from API pages |

Follow that trace for a new input. Add entries to the ends of positional option vectors, never their middles. If the input affects a physical model, update its state layout, right-hand side, Jacobian and output reconstruction together; a parser change alone is insufficient.

A new stress or viscosity law also needs checker ranges, labels in `f_display`, compatibility decisions and tests of limiting behavior. Keep solver paths compatible with base MATLAB. Discuss intended numerical-result changes before updating stored references; see [contributing](../reference/contributing.md).

This worked example explains existing code. It does not add an input or modify solver behavior.
