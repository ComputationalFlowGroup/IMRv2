# Option permutations

`s_test_suite_permutations.m` checks that its enumerated option combinations run when `f_checkers` accepts them. An input that passes the checks and then fails inside a solver is a bug: either the solver must handle it, or `f_checkers` must reject it with an `INPUT ERROR`. The suite needs no stored results, and it always runs in CI.

MATLAB checks were not run for this documentation change. The descriptions below are checked against the suite implementation; numerical outcomes remain for local validation.


## Cases

From [`f_test_cases('permutations')`](tests.md), in the three solver modes. The equilibrium and Gaussian suites already cross the radial, thermal, vapor and stress models; this suite crosses the remaining options:

- **Non-Newtonian viscosity:** every model (`nu_model` 1–7) with every viscous stress model (1–5), with and without medium heat transfer.
- **Driving waveforms:** every `wave_type` (−3 to 3, the measured waveforms included) with every radial model.
- **Integrators and Jacobians:** `method` 3 (CVODES), 15 and 23 with `jacobian` 0–3, every thermal configuration (mass transfer included) and stress 0, 1, 3 and 5.
- **Limits:** no viscosity, surface tension or elasticity (`Re`, `We` or `Ca` = `Inf`) with every stress model.
- **Collapse:** radial 1 and 2, every stress model (stress 5 with `stress0`), alone and with each limit above.
- **Other inputs:** the stress 5 field (`spectral = 1`), `eps3`, a moving start, a given initial pressure and stress, the end time as `tfin`, and the integrator given by name.

Combinations that `f_checkers` rejects are left out automatically. This includes CVODES before MATLAB R2024a, which lacks the `ode` object's `cvodesstiff` solver, and `jacobian = 3` in `f_imr_spectral`.

Most runs use 0 to 2 µs for a 50 µm bubble with five output times and coarse grids. Collapse cases override that geometry and time span. Earlier notes compared its cost to one other suite; that timing was not rerun here.

## Check

The suite requires at least two finite radius values, matching finite times, and completion at the requested final time for both `tvector` and `tfin`. With `tvector`, it also checks the requested output count. An integrator warning with a truncated history therefore cannot silently pass as a completed `tfin` run. Empty suites fail rather than reporting zero passed cases as success.
