# Equilibrium

`s_test_suite_equilibrium.m` checks that a bubble started at equilibrium stays there: `R0 = Req`, no driving, `t` from 0 to 10 µs. This suite needs no stored results. This suite always runs in CI.

Nobody ran MATLAB for this documentation change. The author checked the descriptions below against the suite implementation. You must verify the numerical outcomes on your own MATLAB installation.


## Cases

From [`f_test_cases('equilibrium')`](tests.md), in the three solver modes. The cases are:

- Radial models 1–6.
- The polytropic bubble, bubble heat transfer, and bubble and medium heat transfer.
- With and without vapor.
- Stress models 0–11 (10 and 11, PTT and Giesekus, with a nonzero `eps3` for 11, run only in `imr.solve_spectral`; see [stress and viscosity](../model/stress-viscosity.md#ptt-and-giesekus-stress-10-11); 3, 7 are 2nd-order Fung, 4, 8 are Mooney-Rivlin).

The generator keeps only the combinations the compatibility matrix allows. Before the generator added PTT/Giesekus and the 2nd-order Fung/Mooney-Rivlin springs, this suite had 144 cases per mode (see the disclaimer above).

The physical inputs are `Req = R0 = 100 µm` and `T8 = 298.15 K`.

## Check

A case passes when `max |R - 1| <= 1e-13` over the 100 output times, with `R` scaled by `R0`. Only rounding should move the bubble.
