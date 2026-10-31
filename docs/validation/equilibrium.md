# Equilibrium

`s_test_suite_equilibrium.m` checks that a bubble started at equilibrium stays there: `R0 = Req`, no driving, `t` from 0 to 10 µs. It needs no stored results, and it always runs in CI.

MATLAB checks were not run for this documentation change. The descriptions below are checked against the suite implementation; numerical outcomes remain for local validation.


## Cases

From [`f_test_cases('equilibrium')`](tests.md), in the three solver modes. The cases are:

- Radial models 1–6.
- The polytropic bubble, bubble heat transfer, and bubble and medium heat transfer.
- With and without vapor.
- Stress models 0–11 (6 and 7, PTT and Giesekus, with a nonzero `eps3` for 7, run only in `imr.solve_spectral`; see [stress and viscosity](../model/stress-viscosity.md#ptt-and-giesekus-stress-6-7); 8, 9 are 2nd-order Fung, 10, 11 are Mooney-Rivlin).

Only the combinations the compatibility matrix allows are kept. This was 144 cases per mode before stress 6 and 7 were added, and has not been rerun since 8–11 were added (see the disclaimer above).

The physical inputs are `Req = R0 = 100 µm` and `T8 = 298.15 K`.

## Check

A case passes when `max |R - 1| <= 1e-13` over the 100 output times, with `R` scaled by `R0`. Only rounding should move the bubble.
