# Numerical bug regressions

`tests/suites/s_test_suite_numerics.m` covers the corrections below using base MATLAB, with no stored baseline files. Run from the checkout root:

```bash
cd tests
matlab -batch "imr_test_setup; s_test_suite_numerics"
```

## Driven enthalpy equations

The previous radial 3–6 branches computed wall enthalpy but subtracted the external pressure perturbation directly. This mixed enthalpy and pressure when the external pressure moved away from ambient. The corrected [equations](../theory/radial-equations.md#enthalpy-forms) evaluate the wall and far-field terms with the same EoS, including the far-field density factor in the time derivative.

The suite checks zero acceleration when wall and far-field pressures and their derivatives match. It also compares moving, driven states against independent adaptive quadrature of $1/\rho(p)$ for both Tait and Hugoniot-based models. These tests distinguish the correction from the old approximation. Driven enthalpy results can change; no stored reference was regenerated.

## Maximum time step

`divisions` now limits a step to the **elapsed interval** divided by that count, rather than final absolute time divided by the count. A run from 1000 to 1001 with 40 divisions therefore has `MaxStep = 0.025`, instead of `25.025`. The suite integrates a constant solution so accuracy control does not hide an oversized maximum step. It checks `ode15s`, `ode23tb` and, when available, CVODES on R2024a or later.

## Thermal input limits

The implemented temperature-transform inverse divides by the gas conductivity slope. Unsupported `alpha_g = 0` thermal inputs are rejected before initialization or wall reconstruction. The checks preserve existing `kappa = 1` limits and allow zero vapor conductivity when that component is unused. Profile preparation retains its stricter SI transport requirements.

## Test completion

`f_test_run` now returns the actual time vector as an optional third output. The permutation suite checks the final time for `tfin` and `tvector`, catching an early integrator stop even when a finite radius history is returned. Empty generated suites raise an error. The focused suite checks that an early-return fixture exposes its actual times and that an empty report fails.

## Execution status

The current review adds these regressions to the existing suites:

- **Tait limit:** radial 3 and 4 at `nstate = 1` and nearby exponents, checked against independent quadrature near ambient pressure and exact pressure balance.
- **Histotripsy edges:** zero, fractional and integer exponents; endpoint conventions, zero amplitude, scalar/array shape and derivatives checked against centered differences inside the pulse.
- **Inactive grids:** both solvers, both nodal schemes, and inactive grid counts 0–2. Previously, initialization could build an invalid stencil or divide by a zero collocation order even though the field was disabled. A thermal bubble with no medium grid is also checked.
- **Input limits:** reject nonfinite, complex, fractional or nonscalar grid/division counts and malformed time vectors before allocation or integration. Increasing real row and column time vectors remain accepted.
- **Table orientation:** `s_test_suite_fluids` checks identical saturation results from row and column tables and compares derivatives with the analytic inverse-temperature interpolation, including interior-knot conventions.

MATLAB and Octave were not installed or run in this cloud task. The numerical assertions are prepared for local MATLAB execution and the existing CI workflow; static syntax checks and independent reference calculations do not establish that these MATLAB suites pass. The regression suite still requires its separately approved baseline files.
