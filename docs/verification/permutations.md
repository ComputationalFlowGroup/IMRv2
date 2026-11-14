# Option permutations

`s_test_suite_permutations.m` checks that its enumerated option combinations run, when `imr.configuration.validate_case` accepts them. Consider an input that passes the checks, but then fails inside a solver: this failure is a bug. Either the solver must handle this input, or `imr.configuration.validate_case` must reject this input with an `INPUT ERROR`. This suite needs no stored results. This suite always runs in CI.

Nobody ran MATLAB for this documentation change. The author checked the descriptions below against the suite implementation. You must verify the numerical outcomes on your own MATLAB installation.


## Cases

From [`f_test_cases('permutations')`](tests.md), in the three solver modes. The equilibrium and Gaussian suites already cross the radial, thermal, vapor and stress models. This suite crosses the remaining options:

- **Non-Newtonian viscosity:** every model (`nu_model` 1–7) with every viscous stress model (1–11, Giesekus with a nonzero `eps3`), with and without medium heat transfer.
- **Driving waveforms:** every `wave_type` (−3 to 3, the measured waveforms included) with every radial model.
- **Integrators and Jacobians:** `method` 3 (CVODES), 15 and 23 with `jacobian` 0–3, every thermal configuration (mass transfer included) and stress 0, 1, 3, 4, 5, 7, 8, 9, 10 and 11.
- **Limits:** no viscosity, surface tension or elasticity (`Re`, `We` or `Ca` = `Inf`), with every stress model. This suite skips one combination: PTT and Giesekus (10, 11) with no viscosity (`Re = Inf`). The suite skips this combination because their nonlinear terms multiply `Re8` directly, which gives `0*Inf` from rest. The 2nd-order Fung and Mooney-Rivlin springs (3, 4, 7, 8) only divide by `Re8`, so these springs stay in the sweep.
- **Collapse:** radial 1 and 2, stress 0–6 and 9 (stress 9 with `stress0`), alone and with each limit above. This suite leaves out stress 7, 8, 10 and 11, because every mode rejects their collapse. Widening this loop to include them would add no cases.
- **Other inputs:** the stress 9 field (`spectral = 1`), `eps3`, a moving start, a given initial pressure and stress, the end time as `tfin`, and the integrator given by name.

The suite automatically leaves out combinations that `imr.configuration.validate_case` rejects. This includes CVODES before MATLAB R2024a, which lacks the `ode` object's `cvodesstiff` solver. This also includes `jacobian = 3` in `imr.solve_spectral`.

Most runs use 0 to 2 µs, for a 50 µm bubble, with five output times and coarse grids. Collapse cases override this geometry and time span. Earlier notes compared this suite's cost to one other suite. Nobody reran that timing here.

## Check

The suite requires at least two finite radius values. The suite requires matching finite times. The suite requires completion at the requested final time, for both `tvector` and `tfin`. With `tvector`, the suite also checks the requested output count. So an integrator warning with a truncated history cannot silently pass as a completed `tfin` run. Empty suites fail; the suite does not report zero passed cases as success.
