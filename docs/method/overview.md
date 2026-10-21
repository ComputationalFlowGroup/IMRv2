# Numerical workflow

Both solvers turn a configured physical model into a coupled system of ordinary differential equations. They integrate that state and reconstruct the returned fields.

| Stage | Implementation | Result |
|---|---|---|
| Convert SI, if needed | `f_dim2nondim` | Dimensionless pairs and reference scales |
| Defaults and overrides | `default_case`, `f_imr_inputs`, `f_call_params` | Grouped option vectors and initial state |
| Validate | `f_checkers` | Compatible settings or an error |
| Spatial operators | `f_diff_mat` or spectral setup | Bubble/medium derivatives and boundary maps |
| Evaluate derivative | Nested solver right-hand side and physics helpers | Radial, thermal, vapor and stress derivatives |
| Integrate | `f_odesolve` and selected Jacobian | State history |
| Reconstruct | Output assembly and wall closure | Public arrays and input record |

## State layout

The nodal state follows `[R; Rdot; P; theta; Tm; kv; Z]`, with inactive fields omitted. Stress 3–4 with medium heat add the local polymer stress after `Z1`. The spectral solver stores expansion coefficients for its thermal and local stress fields. `theta` is the [Kirchhoff variable](time-integration.md#thermal-variables), not the returned temperature.

The local polymer field supplies heating, while its stress integral drives radial motion. The [stress-heating theory](../theory/stress-heating.md) states that coupling and the collapse initialization approximation. The [output guide](../running/output.md) describes field reconstruction and wall values.

## Implementation considerations

Constant medium and thermal-map factors are precomputed; constant forcing skips `f_pinfinity`. Nodal nested functions share names with their parent workspace, so local loop names must avoid collisions; the output loop uses `it_out`.

The nodal solver was formerly `f_imr_fd`. Historical reproduction and closure changes are recorded in [reported comparisons](../validation/convergence.md), rather than asserted as current numerical validation. Grid resolution needs attention for large relaxation time scales, sharp waveforms and collapse.

See [discretization](discretization.md), [wall boundaries](wall-boundaries.md), [integration](time-integration.md), and [initial stress](initial-stress.md). The [code map](../reference/code-map.md) shows function calls.
