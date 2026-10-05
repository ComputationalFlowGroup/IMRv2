# `s_test_suite_equilibrium.m`

## Overview

`s_test_suite_equilibrium.m` checks that a bubble started at equilibrium stays there: `R0 = Req`, no driving, `t` from 0 to 10 µs. It needs no stored results, and it always runs in CI.

## Cases

From [`f_test_cases('equilibrium')`](f_test_cases.md), in the three solver modes. The cases are:

- Radial models 1–6.
- The polytropic bubble, bubble heat transfer, and bubble and medium heat transfer.
- With and without vapor.
- Stress models 0–5.

Only the combinations the compatibility matrix allows are kept: 144 cases per mode.

The physical inputs are `Req = R0 = 100 µm` and `T8 = 298.15 K`.

## Check

A case passes when `max |R - 1| <= 1e-13` over the 100 output times, with `R` scaled by `R0`. Only rounding should move the bubble.
