# `f_pinfinity.m` Documentation

## Overview
`f_pinfinity` calculates the time-dependent far-field pressure `p8` and its time derivative `p8dot` that drive the bubble. All quantities are dimensionless.

## Major Sections
1. **Input Unpacking**: Extracts waveform parameters from `vararg`.
2. **Measured Waveforms**: If `wave_type < 0`, evaluates the piecewise polynomials `wave_poly` and `wave_dpoly` built from `d_hn`, `d_mn` or `d_ml`. These are the measured waveforms of [Hersey, Rodriguez and Johnsen (2023)](https://doi.org/10.1121/10.0017342); see [`docs/waveform`](../waveform/README.md).
3. **Predefined Waveforms**:
    - **Impulse** (`wave_type = 0`)
    - **Gaussian** (`wave_type = 1`)
    - **Histotripsy** (`wave_type = 2`)
    - **Heaviside Impulse** (`wave_type = 3`)
Each waveform is evaluated inline, without nested function calls, since the solvers call `f_pinfinity` on every right-hand-side evaluation; the Gaussian evaluates its exponential once for `p8` and `p8dot`. For the impulse (`wave_type = 0`) the solvers skip the call, because the pressure is the constant `ee`.

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

## Notes
- The solvers pass `wave_poly` and `wave_dpoly` in `vararg`; `f_call_params` builds them for `wave_type < 0`.
- Gaussian and histotripsy waveforms depend on `ee` (amplitude `pA/P8`), `om` (frequency `omega*t0`), `tw` (width `TW/t0`), `dt` (delay `DT/t0`), and `mn` (histotripsy power).
- The default `ee = 0` means no driving pressure.
