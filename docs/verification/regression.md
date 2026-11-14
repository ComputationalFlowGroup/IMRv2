# Regression


Nobody ran MATLAB for this documentation change. The author checked the descriptions below against the suite implementation. You must verify the numerical outcomes locally.

`s_test_suite_regression.m` runs the regression cases in the three solver modes. It compares `R(t)` with the stored results in `tests/reference/<mode>.mat`. It catches unintended changes in the results.

## Cases

From [`f_test_cases('regression')`](tests.md). The generator keeps only the combinations the compatibility matrix allows. Before the generator added the 2nd-order Fung and Mooney-Rivlin springs (stress 3, 4, 7, 8), this count was 230 cases each for `nodal_fd` and `nodal_cheb`, and 155 for `spectral`. The stored references in `tests/reference/` do not reflect the new count yet. Nobody has regenerated them (see the disclaimer above). A local MATLAB run must refresh them.

- **Radial, thermal and stress models, for two materials.** The cases are:

    - Radial models 1–6.
    - Five thermal configurations: polytropic, bubble, bubble and medium, bubble with mass transfer, and bubble and medium with mass transfer.
    - Stress models 0–9 (3, 7 are 2nd-order Fung, 4, 8 are Mooney-Rivlin, see [stress and viscosity](../model/stress-viscosity.md#2nd-order-fung-stress-3-7)).

    The two materials are a soft, weakly viscous one (`mu = 1e-4 Pa s`, `G = 100 Pa`, `alphax = 1e-3`, `C2 = 0.1`, `lambda1 = 1e-7 s`) and a stiff, viscous one (`mu = 0.1 Pa s`, `G = 2.5 kPa`, `alphax = 1`, `C2 = 0.5`, `lambda1 = 1 ms`).

    The runs are:

    - Without mass transfer: a 50 µm bubble released from four times its equilibrium radius and followed for 12 µs.
    - The polytropic bubble: released from twice its equilibrium radius. It loses no heat, so its collapse is far sharper.
    - With mass transfer: a 200 µm bubble released from four times its equilibrium radius and followed for 30 µs.

    These runs exercise every term. They are not so violent that the result mostly tests the integrator instead of the physics. The earlier suites released bubbles from ten times the equilibrium radius. At that radius, about half of the cases stop at the first collapse.

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

Regenerate them only when the team intends a change to the results and agrees to it (see Contributing).

The maintainers will recreate the missing baselines locally in MATLAB. From the checkout root, run the independent prechecks first. Generate the baselines only if the prechecks pass. Then run the regression suite:

```bash
matlab -batch "cd tests; run_verification('generate-and-test')"
```

For generation only, use `run_verification('generate')` from `tests/`. To replace existing files deliberately, pass `true` as the third argument, for example `run_verification('generate-and-test','all',true)`. See [all runner modes and suite selection](../../tests/README.md#one-command-for-all-tests). Ordinary test mode never invokes the generator.

The existing script remains available from the checkout root:

```bash
cd toolchain
matlab -batch "s_generate_references"
```

The default `overwrite = false` allows creating absent files. Set `overwrite = true` only when you intend to replace existing baselines. Review the cross-mode report. Record the code revision alongside the generated files. The `info` field records the generation time and the MATLAB version. Then, from the checkout root, run the comparison:

```bash
cd tests
matlab -batch "imr_test_setup; s_test_suite_regression"
```

- **Overwrite guard.** It refuses to overwrite existing files unless `overwrite = true` is set in the script.
- **Parallel runs.** It runs the cases with `parfor`, which uses the default parallel pool when the Parallel Computing Toolbox is available and runs in series otherwise.
- **Completed simulations.** If any case fails, stops early, or returns mismatched output times, the script lists the failures and publishes nothing. Otherwise, it stages all files in a temporary folder and then moves them into place. If it fails, it cleans up the temporary staging files. It checks overwrite protection again before it publishes the files.
- **Cross-mode report.** For the cases the modes share, it reports the largest relative difference between them. It lists any case above `1e-2`. This is a sanity check on the new results, not a test. The grids differ between modes. Violent runs can differ more.

A new case added to `f_test_cases` fails the suite until its result is stored.

## Baseline availability

The current checkout does not contain `tests/reference/nodal_fd.mat`, `nodal_cheb.mat`, or `spectral.mat`. The regression suite therefore reports missing stored results. It has not passed here. The project has removed the unused older `.mat` histories. No legacy archive remains.

Generating new results from the current code creates a new reference. This process cannot establish agreement with an independently approved prior version. Review the generated results with the analytic suites and the cross-mode comparisons before you accept them as baselines. A regression run against freshly generated files checks storage and repeatability. It does not check independent numerical accuracy.
