# `f_pvsat.m` Documentation

## Overview
`[Pv, c1, c2] = f_pvsat(T)` returns the saturated water vapor pressure as a function of temperature using an empirical exponential fit, `Pv = c1*exp(c2/T)`. The constants are returned as optional outputs so that this file is their single source. `f_dim2nondim` uses them for the dimensionless vapor pressure `pv` and its exponent `thetav = -c2/T8`.

## Major Sections
1. **Empirical Formula**: Applies `Pv = 1.17e11 * exp(-5200/T)`.

## Inputs

| Name | Description        | Units |
|------|--------------------|-------|
| `T`  | Temperature       | K     |

## Outputs

| Name | Description                   | Units |
|------|-------------------------------|-------|
| `Pv` | Saturated vapor pressure      | Pa    |
| `c1` | Prefactor, `1.17e11`          | Pa    |
| `c2` | Exponent constant, `-5200`    | K     |

## Notes
- Source: A. Preston, 2004 (Thesis).
- Valid for typical biological temperature ranges (290–310 K).
