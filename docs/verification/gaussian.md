# Gaussian pulse

`s_test_suite_gaussian.m` drives a bubble at equilibrium with a Gaussian pressure pulse. The suite checks that the bubble returns to its equilibrium radius. This suite needs no stored results.

Nobody ran MATLAB for this documentation change. The author checked the descriptions below against the suite implementation. You must verify the numerical outcomes on your own MATLAB installation.


## Cases

From [`f_test_cases('gaussian')`](tests.md), in the three solver modes. The cases are:

- Radial models 1–6.
- The polytropic bubble, bubble heat transfer, and bubble and medium heat transfer.
- With and without vapor.
- The viscous stress models 1–11 (10 and 11, PTT and Giesekus, with a nonzero `eps3` for 11, run only in `imr.solve_spectral`; see [stress and viscosity](../model/stress-viscosity.md#ptt-and-giesekus-stress-10-11); 3, 7 are 2nd-order Fung, 4, 8 are Mooney-Rivlin).

The generator keeps only the combinations the compatibility matrix allows. Before the generator added PTT/Giesekus and the 2nd-order Fung/Mooney-Rivlin springs, this suite had 120 cases per mode (see the disclaimer above).

The physical inputs are:

- `Req = R0 = 10 µm`, with `mu = 0.5 Pa s`.
- A Gaussian pulse (`wave_type = 1`) of amplitude `pA = 1e5 Pa`, width `tw = 0.1 ms`, centered at `dt = 0.3 ms`.
- `t` from 0 to 2 ms.

## Check

A case passes when `|R(end) - R(1)| <= 1e-2`.
