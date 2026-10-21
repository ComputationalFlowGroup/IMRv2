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

For histotripsy (`wave_type = 2`), inside $|t-\mathtt{dt}|<\pi/\mathtt{om}$ the pressure is $\mathtt{ee}\cos^{2\mathtt{mn}}(\mathtt{om}(t-\mathtt{dt})/2)$. It vanishes outside that interval. A zero exponent gives a rectangular pulse, including its endpoints. The returned derivative is zero at the endpoints by a pointwise convention, as for the Heaviside waveform; this does not represent the distributional derivative of a step. For $0<\mathtt{mn}<1/2$, the interior derivative has an unbounded one-sided limit, and at $\mathtt{mn}=1/2$ it jumps. Use $\mathtt{mn}>1/2$ when a continuously differentiable pulse is required. The half-angle evaluation avoids cancellation near the edges.

For measured waveforms, see the [topic guide](../../running/waveforms.md). Numerical measurements are retained as [historical reports](../../validation/convergence.md).
