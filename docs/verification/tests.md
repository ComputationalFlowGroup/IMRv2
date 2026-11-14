# Test suites


Nobody ran MATLAB in this cloud task. The author checked the descriptions below against the suite implementation. You must verify the numerical outcomes on your own MATLAB installation.

The test suites live in `tests/suites/`. From `tests/`, call `imr_test_setup` once. This call exposes the suites, helpers, reference functions and source. The [test directory guide](../../tests/README.md) explains the layout. It also explains the difference between frozen references and golden trajectories. The `tests.yml` workflow runs each suite in its own MATLAB process on every pull request. The workflow fails if any suite fails.

| Suite | Checks | Stored results |
|---|---|---|
| [`s_test_suite_equilibrium`](equilibrium.md) | A bubble at equilibrium stays there | None |
| [`s_test_suite_gaussian`](gaussian.md) | A bubble returns to equilibrium after a Gaussian pressure pulse | None |
| [`s_test_suite_permutations`](permutations.md) | Every accepted combination of the solver options runs | None |
| [`s_test_suite_regression`](regression.md) | `R(t)` matches the stored results | `tests/reference/` |
| [`s_test_suite_diffmat`](operators.md) | The differentiation operators of `imr.numerics.differentiation` | None |
| [`s_test_suite_polytropic`](polytropic.md) | Constant-vapor pressure derivative and histories | None |
| [`s_test_suite_numerics`](numerics.md) | Driven enthalpy, time-step limits and input guards | None |
| [`s_test_suite_settings`](settings.md) | Explicit acceptance/rejection, phase limits and initial-state guards | None |
| [`s_test_suite_fluids`](../reference/defaults.md) | Water properties, defaults, overrides and custom profiles | Frozen old default values |
| [`s_test_suite_namespace`](../reference/namespaces.md) | Namespace resolution, legacy calls and measured data | None |
| [`s_test_suite_solver`](../reference/solve.md) | Unified selection, replay and optional acceleration | None |
| [`s_test_suite_verification`](../../tests/README.md#one-command-for-all-tests) | Runner modes, workspace isolation and golden-data safeguards | Temporary synthetic fixtures only |
| [`s_test_suite_structure`](../reference/internal-contracts.md#verification) | Prepared cases, states, independent Jacobians and output boundaries | None |

`s_test_suite_solver` checks struct/pair selection precedence, conflicting wrapper selections, existing compatibility errors, wrapper equivalence and replay. The suite compares acceleration with independent Rayleigh–Plesset and pressure-form Keller–Miksis equations, under Gaussian forcing at nonzero absolute times. The suite also compares acceleration with differences of sampled velocity. The suite also covers the available time integrators, adaptive output times, equilibrium for all six radial models, coupled thermal/vapor/stress fields and quiet RHS evaluation. The suite does not generate or need golden data.

The generated equilibrium, Gaussian, permutation and regression suites use three solver modes. The separate polytropic and fluid suites also exercise these modes. The operator suite checks matrices directly:

| Mode | Solver |
|---|---|
| `nodal_fd` | `imr.solve_nodal`, finite differences (`scheme = 0`) |
| `nodal_cheb` | `imr.solve_nodal`, Chebyshev collocation (`scheme = 1`) |
| `spectral` | `imr.solve_spectral` |

## `f_test_cases(suite)`

This function returns the cases of `'equilibrium'`, `'gaussian'`, `'regression'` or `'permutations'`, for all three modes. The function builds each case from the model features. The function keeps a case only if the solver's own input checks accept it (`imr.configuration.prepare_case` and [`imr.configuration.validate_case`](../reference/functions/validate_case.md)). The generator enumerates selected combinations. The generator filters these combinations through the rules behind the [compatibility matrix](../reference/compatibility.md). The generator does not test every possible parameter value. The generator does not prove that every pairwise-compatible set of features is valid.

Each case has these fields:

- `key`: a name for the case within its mode, such as `r2_bubble+medium_s3_mat1`. Stored results are looked up by it.
- `mode`: the solver mode.
- `solver`: the solver function handle.
- `args`: the dimensionless name/value pairs, grid included.

You give the inputs in SI units. `f_test_cases` converts them with [`imr.nondimensionalize`](../reference/nondimensionalize.md). The grids are:

| Suite | `nodal_fd` | `nodal_cheb` | `spectral` |
|---|---|---|---|
| Equilibrium | `Nt = Mt = 150` | 20 | 12 |
| Gaussian | 50 | 20 | 12 |
| Regression | 70 | 20 | 12 |
| Permutations | 20 | 12 | 8 (`Nv = 8`) |

Finite differences need finer grids, because finite differences are only second-order accurate. The Chebyshev and spectral expansions converge spectrally, so they need coarser grids. Except in permutations, the spectral mode sets `Nv = 32` stress modes. This mode count applies to the stress 9–11 field and to the Zener field of stress 5–8 with medium heat transfer. The default `Nv` is 150.

**Historically reported case counts** (`nodal_fd`, `nodal_cheb`, `spectral`):

| Suite | Cases |
|---|---|
| Equilibrium | 144, 144, 144 (before PTT/Giesekus and the 2nd-order Fung/Mooney-Rivlin springs were added; not rerun here) |
| Gaussian | 120, 120, 120 (before PTT/Giesekus and the 2nd-order Fung/Mooney-Rivlin springs were added; not rerun here) |
| Regression | 230, 230, 155 (before the 2nd-order Fung/Mooney-Rivlin springs were added; stored references need a local MATLAB regeneration) |

## Helpers

- `f_test_run(c)` runs one case. It hides the solver's console output. It checks all seven physical outputs for real, finite values. It returns `[R, msg, t]`. On success, `msg` is empty, and the optional third output `t` holds the actual time history. If the run fails, `msg` holds the error message instead, so one failing case does not stop a suite.
- `f_test_report(name, cases, failed, msg)` prints the passes for each mode. It prints every failure with its mode and key. It then raises an error if any case failed, or if no cases were generated. In a batch run (`matlab -batch`, as in CI), this error exits MATLAB with a nonzero code. In an interactive session, this error only stops the script.

## Run the suites

From the checkout root, use a shell with licensed MATLAB on its path. The combined runner discovers all suites. The runner isolates each suite's function workspace within one MATLAB process. The runner continues after a failure. The runner returns a nonzero batch status if any stage fails:

```bash
matlab -batch "run('tests/s_run_verification.m')"
```

From `tests/`, use `run_verification('tests',{'diffmat','solver'})` to run only a subset. Use `run_verification('generate-and-test')` for independent prechecks, followed by explicit generation and regression. See [runner modes, selection and overwrite protection](../../tests/README.md#one-command-for-all-tests).

For separate MATLAB processes, use the following shell loop instead. This loop also discovers the available suites. This loop also preserves failure status:

```bash
cd tests
status=0
for file in suites/s_test_suite_*.m; do
    suite=${file##*/}
    suite=${suite%.m}
    matlab -batch "imr_test_setup; $suite" || status=1
done
exit "$status"
```

Read the per-mode counts and the final report. When a case generator rejects an incompatible combination, that rejection is an exclusion. An exclusion is not a passing simulation. Trust the counts printed by the current run. The static counts kept on these pages are only historical reports. These static counts are not verified counts for every MATLAB release.

The regression suite needs [approved baseline files](regression.md#baseline-availability). This checkout does not contain these files. Other suites do not depend on these files. These tests check numerical consistency. No bundled external-paper comparison proves scientific accuracy to a stated number of digits.
