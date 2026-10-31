# Regression


MATLAB checks were not run for this documentation change. The descriptions below are checked against the suite implementation; numerical outcomes remain for local validation.

`s_test_suite_regression.m` runs the regression cases in the three solver modes and compares `R(t)` with the stored results in `tests/reference/<mode>.mat`. It catches unintended changes of the results.

## Cases

From [`f_test_cases('regression')`](tests.md). Only the combinations the compatibility matrix allows are kept: this was 230 cases each for `nodal_fd` and `nodal_cheb`, and 155 for `spectral`, before the 2nd-order Fung and Mooney-Rivlin springs (stress 3, 4, 7, 8) were added to this generator; the stored references in `tests/reference/` and the new count have not been regenerated here (see the disclaimer above) and need a local MATLAB run to refresh.

- **Radial, thermal and stress models, for two materials.** The cases are:

    - Radial models 1–6.
    - Five thermal configurations: polytropic, bubble, bubble and medium, bubble with mass transfer, and bubble and medium with mass transfer.
    - Stress models 0–9 (3, 7 are 2nd-order Fung, 4, 8 are Mooney-Rivlin, see [stress and viscosity](../model/stress-viscosity.md#2nd-order-fung-stress-3-7)).

    The two materials are a soft, weakly viscous one (`mu = 1e-4 Pa s`, `G = 100 Pa`, `alphax = 1e-3`, `C2 = 0.1`, `lambda1 = 1e-7 s`) and a stiff, viscous one (`mu = 0.1 Pa s`, `G = 2.5 kPa`, `alphax = 1`, `C2 = 0.5`, `lambda1 = 1 ms`).

    The runs are:

    - Without mass transfer: a 50 µm bubble released from four times its equilibrium radius and followed for 12 µs.
    - The polytropic bubble: released from twice its equilibrium radius. It loses no heat, so its collapse is far sharper.
    - With mass transfer: a 200 µm bubble released from four times its equilibrium radius and followed for 30 µs.

    These runs exercise every term without being so violent that the result mostly tests the integrator. From ten times the equilibrium radius, as the earlier suites used, about half of the cases stop at the first collapse.

- **Other features**, on Keller–Miksis with bubble and medium heat transfer:

    - The seven non-Newtonian viscosity models.
    - The Gaussian, histotripsy and Heaviside waveforms.
    - The UCM / Oldroyd-B stress field (`spectral = 1`, spectral solver only).
    - Collapse with stress 0, 1, 3, 4, 5 and 6, with the soft material, up to the first collapse and rebound (20 µs). This is nodal solver only, since collapse needs mass transfer.

## Check

A case fails if:

- its run fails;
- it has no stored result;
- its output size differs from the stored one; or
- `norm(R./R_ref - 1, 2) > 1e-4` over its 100 output times.

The suite also notes stored results that no case uses any more.

## Generating the stored results

The shared `tests/helpers/generate_references.m` runs every regression case and writes `tests/reference/nodal_fd.mat`, `nodal_cheb.mat` and `spectral.mat`. The controller and the existing `toolchain/s_generate_references.m` script call this same implementation. Each file holds:

- `keys`: the case keys.
- `R`: the results.
- `info`: the date and program version.

Regenerate them only when a change of results is intended and agreed (see Contributing).

The maintainers will recreate the missing baselines locally in MATLAB. From the checkout root, run independent prechecks, generate the baselines only if they pass, and then run regression:

```bash
matlab -batch "cd tests; run_validation('generate-and-test')"
```

For generation only, use `run_validation('generate')` from `tests/`. To replace existing files deliberately, pass `true` as the third argument, for example `run_validation('generate-and-test','all',true)`. See [all runner modes and suite selection](../../tests/README.md#one-command-for-all-tests). Ordinary test mode never invokes the generator.

The existing script remains available from the checkout root:

```bash
cd toolchain
matlab -batch "s_generate_references"
```

The default `overwrite = false` allows creating absent files. Set `overwrite = true` only when intentionally replacing existing baselines. Review the cross-mode report and record the code revision alongside the generated files; `info` records the generation time and MATLAB version. Then, from the checkout root, run the comparison:

```bash
cd tests
matlab -batch "imr_test_setup; s_test_suite_regression"
```

- **Overwrite guard.** It refuses to overwrite existing files unless `overwrite = true` is set in the script.
- **Parallel runs.** It runs the cases with `parfor`, which uses the default parallel pool when the Parallel Computing Toolbox is available and runs in series otherwise.
- **Completed simulations.** If any case fails, stops early or returns mismatched output times, it lists the failures and publishes nothing. Otherwise it stages all files in a temporary folder and then moves them into place. Temporary staging files are cleaned up on failure, and overwrite protection is checked again before publishing.
- **Cross-mode report.** For the cases the modes share, it reports the largest relative difference between them, and lists any case above `1e-2`. This is a sanity check on the new results, not a test: the grids differ, and violent runs may differ more.

A new case added to `f_test_cases` fails the suite until its result is stored.

## Baseline availability

The current checkout does not contain `tests/reference/nodal_fd.mat`, `nodal_cheb.mat` or `spectral.mat`. The regression suite therefore reports missing stored results; it has not passed here. The unused older `.mat` histories have been removed, and no legacy archive is retained.

Generating new results from the current code creates a new reference; it cannot establish agreement with an independently approved prior version. Review the generated results using the analytic suites and cross-mode comparisons before accepting them as baselines. A regression run against freshly generated files checks storage and repeatability, not independent numerical accuracy.
