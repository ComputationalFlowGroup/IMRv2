# Numerical bug regressions

`tests/suites/s_test_suite_numerics.m` covers the corrections below using base MATLAB, with no stored baseline files. Run from the checkout root:

```bash
cd tests
matlab -batch "imr_test_setup; s_test_suite_numerics"
```

## Driven enthalpy equations

The previous radial 3–6 branches computed wall enthalpy correctly, but subtracted the external pressure perturbation directly. This direct subtraction mixed enthalpy and pressure, when the external pressure moved away from ambient. The corrected [equations](../model/radial-dynamics.md#enthalpy-forms) evaluate the wall and far-field terms with the same EoS. The corrected equations include the far-field density factor in the time derivative.

The suite checks for zero acceleration, when the wall and far-field pressures and their derivatives match. The suite also compares moving, driven states against independent adaptive quadrature of $1/\rho(p)$, for both Tait and Hugoniot-based models. These tests distinguish the correction from the old approximation. Driven enthalpy results can change. Nobody regenerated a stored reference for this correction.

## Maximum time step

`divisions` now limits a step to the **elapsed interval** divided by that count. The previous version divided the final absolute time by that count instead. For example, a run from 1000 to 1001 with 40 divisions now has `MaxStep = 0.025`. The previous version gave `MaxStep = 25.025` for the same run. The suite integrates a constant solution, so accuracy control cannot hide an oversized maximum step. The suite checks `ode15s`, `ode23tb` and, when available, CVODES on R2024a or later.

## Thermal input limits

The implemented temperature-transform inverse divides by the gas conductivity slope. The suite rejects unsupported `alpha_g = 0` thermal inputs before initialization or wall reconstruction. The checks preserve the existing `kappa = 1` limits. The checks allow zero vapor conductivity when the vapor component is unused. Profile preparation keeps its stricter SI transport requirements.

## Test completion

`f_test_run` now returns the actual time vector as an optional third output. The permutation suite checks the final time for `tfin` and `tvector`. This check catches an early integrator stop, even when the solver returns a finite radius history. Empty generated suites raise an error. The focused suite checks two things: that an early-return fixture exposes its actual times, and that an empty report fails.

## Execution status

The current review adds these regressions to the existing suites:

- **Tait limit:** radial 3 and 4 at `nstate = 1` and nearby exponents, checked against independent quadrature near ambient pressure and exact pressure balance.
- **Histotripsy edges:** zero, fractional and integer exponents. Endpoint conventions, zero amplitude, scalar/array shape and derivatives checked against centered differences inside the pulse.
- **Inactive grids:** both solvers, both nodal schemes, and inactive grid counts 0–2. Previously, initialization could build an invalid stencil, or divide by a zero collocation order, even though the field was disabled. The suite also checks a thermal bubble with no medium grid.
- **Input limits:** the suite rejects nonfinite, complex, fractional or nonscalar grid/division counts, and malformed time vectors, before allocation or integration. The suite still accepts increasing real row and column time vectors.
- **Table orientation:** `s_test_suite_fluids` checks identical saturation results from row and column tables and compares derivatives with the analytic inverse-temperature interpolation, including interior-knot conventions.

Nobody installed or ran MATLAB or Octave in this cloud task. The numerical assertions are ready for local MATLAB execution and for the existing CI workflow. Static syntax checks and independent reference calculations do not prove that these MATLAB suites pass. The regression suite still needs its separately approved baseline files.
