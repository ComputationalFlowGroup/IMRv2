# Troubleshooting

The parser and `imr.configuration.validate_case` report invalid inputs before integration. Errors below are quoted from the implementation; placeholders such as `%d` depend on the case.

| Message or identifier | Meaning | Fix |
|---|---|---|
| `input error: unmatched inputs` | A name/value pair has no partner | Supply one value for every name |
| `INPUT ERROR: unrecognized input(s): ...` | Name not accepted by the solver | Use [accepted inputs](../model/parameters.md); convert SI inputs with `imr.nondimensionalize` |
| `INPUT ERROR: casefile is no longer supported.` | Removed case-file interface | Use a struct or cell array |
| `INPUT ERROR: Only tvector or tfin can be specified, not both` | Conflicting explicit time inputs | Supply one; see [configuration](configuration.md#time-inputs) |
| `INPUT ERROR: tvector must be finite and strictly increasing` | Invalid time vector | Use at least two finite ascending values |
| `INPUT ERROR: tvector must start at a nonnegative time` | Negative start time | Start at zero or later |
| `INPUT ERROR: medtherm = 1 requires bubtherm = 1` | Medium heat without bubble heat | Enable bubble heat or disable medium heat |
| `INPUT ERROR: masstrans = 1 requires vapor = 1` | Mass transfer without vapor | Enable vapor |
| `INPUT ERROR: masstrans = 1 requires bubtherm = 1` | Isothermal vapor diffusion is unavailable | Enable bubble heat |
| `INPUT ERROR: collapse must have full model with vapor` | Collapse prerequisites absent | Use the [collapse guide](collapse.md) |
| `INPUT ERROR: stress0 must have %d element(s) for stress = %d` | Wrong stress-history size | Use the size for the selected stress representation |
| `INPUT ERROR: De = %g is below 1e-4` | Relaxation model is in its small-De limit | Use the simpler stress model named by the error |
| `INPUT ERROR: nu_model = %d needs a viscous stress model` | Non-Newtonian choice with stress 0 | Use any nonzero stress model |
| `INPUT ERROR: stress = 6 (PTT) and stress = 7 (Giesekus) have no closed-form ...` | PTT or Giesekus requested in `imr.solve_nodal` | Use `imr.solve_spectral` |
| `INPUT ERROR: method = 3 (CVODES) needs the ode object` | CVODES unavailable | Use method 15 or 23, or MATLAB R2024a+ with CVODES |
| `INPUT ERROR: jacobian = 3 (batched finite differences) is available in imr.solve_nodal only` | Spectral solver cannot batch the Jacobian | Use Jacobian 2 there |
| `INPUT ERROR: masstrans = 1 is not available in the spectral solver` | Unsupported solver feature | Use nodal solver |
| `INPUT ERROR: scheme selects the differentiation matrices of imr.solve_nodal.` | Nodal option used for spectral solver | Omit scheme for spectral solver |
| `IMR:dimensional` | Values look like SI inputs | Verify scaling before ignoring this warning |
| `IMR:parameters` | Suspicious values or an automatic setting change | Read the message; check boiling, dashpot viscosity, inactive viscosity or waveform extrapolation |
| `no stored result, run s_generate_references` | Regression baseline missing | Obtain approved references; do not silently regenerate them |
| `output size differs from the stored result` | Changed regression output shape | Compare configuration and reference revision |

Model combinations can also fail the exact [compatibility rules](../reference/compatibility.md). The [checker reference](../reference/functions/f_checkers.md) lists allowed ranges and conditions, including model-specific viscosity exponents, grid sizes and positive transport constants.

## During integration

If the mass-transfer wall residual has no sign change, `imr.solve_nodal` stops with a message beginning `imr.solve_nodal: the mass-transfer wall balance does not change sign on kv_w`. Inspect the reported time, state and saturation range. A violent Rayleigh–Plesset collapse can leave the saturation model's admissible range; the compatibility page records this limitation. Do not bypass the bracket check to force a run.

If initial-stress shooting cannot reach maximum radius with a subsonic start, provide a physically justified `stress0` or change the case. See [initial stress](../method/initial-stress.md). Grid refinement and model validity both matter; a tighter time tolerance alone cannot repair an invalid state.

For an unlisted error, report the exact solver call, full message, MATLAB version and platform through the [contribution guide](../reference/contributing.md#reporting-bugs).
