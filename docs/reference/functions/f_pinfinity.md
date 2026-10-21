# `f_pinfinity`

`f_pinfinity` calculates the time-dependent far-field pressure `p8` and its time derivative `p8dot` that drive the bubble. All quantities are dimensionless.

## Inputs

| Name        | Description                                     | Units   |
|-------------|-------------------------------------------------|---------|
| `t`         | Time, scalar or array                          | –       |
| `vararg`    | Cell `{[om, ee, tw, dt, mn, wave_type], wave_poly, wave_dpoly}` | –       |

## Outputs

| Name     | Description                                 | Units   |
|----------|---------------------------------------------|---------|
| `p8`     | Far-field pressure at time `t`, `p/P8`      | –       |
| `p8dot`  | Time derivative of `p8`                     | –       |

For the algorithm and assumptions, see the [topic guide](../../running/waveforms.md). Numerical measurements are retained as [historical reports](../../validation/convergence.md).
