# Gaussian pulse

`s_test_suite_gaussian.m` drives a bubble at equilibrium with a Gaussian pressure pulse and checks that it returns to its equilibrium radius. It needs no stored results.

MATLAB checks were not run for this documentation change. The descriptions below are checked against the suite implementation; numerical outcomes remain for local validation.


## Cases

From [`f_test_cases('gaussian')`](tests.md), in the three solver modes. The cases are:

- Radial models 1–6.
- The polytropic bubble, bubble heat transfer, and bubble and medium heat transfer.
- With and without vapor.
- The viscous stress models 1–11 (6 and 7, PTT and Giesekus, with a nonzero `eps3` for 7, run only in `imr.solve_spectral`; see [stress and viscosity](../model/stress-viscosity.md#ptt-and-giesekus-stress-6-7); 8, 9 are 2nd-order Fung, 10, 11 are Mooney-Rivlin).

Only the combinations the compatibility matrix allows are kept. This was 120 cases per mode before stress 6 and 7 were added, and has not been rerun since 8–11 were added (see the disclaimer above).

The physical inputs are:

- `Req = R0 = 10 µm`, with `mu = 0.5 Pa s`.
- A Gaussian pulse (`wave_type = 1`) of amplitude `pA = 1e5 Pa`, width `tw = 0.1 ms`, centered at `dt = 0.3 ms`.
- `t` from 0 to 2 ms.

## Check

A case passes when `|R(end) - R(1)| <= 1e-2`.
