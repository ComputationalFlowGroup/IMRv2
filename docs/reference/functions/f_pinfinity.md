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
| `p8`     | Far-field pressure perturbation at time `t`; absolute dimensionless pressure is `1 + p8` | – |
| `p8dot`  | Time derivative of `p8`                     | –       |

The [waveform guide](../../running/waveforms.md) defines every analytic and measured choice, signs, support, derivatives and endpoint conventions. The function evaluates each waveform elementwise without changing the time-array shape. Numerical measurements are retained as [historical reports](../../validation/convergence.md).
