# Test suites


MATLAB checks were not run for this documentation change. The descriptions below are checked against the suite implementation; numerical outcomes remain for local validation.

The test suites live in `toolchain/` and run from that folder. The tests workflow `tests.yml` runs each one in its own MATLAB process on every pull request, and fails if any suite fails.

| Suite | Checks | Stored results |
|---|---|---|
| [`s_test_suite_equilibrium`](equilibrium.md) | A bubble at equilibrium stays there | None |
| [`s_test_suite_gaussian`](gaussian.md) | A bubble returns to equilibrium after a Gaussian pressure pulse | None |
| [`s_test_suite_permutations`](permutations.md) | Every accepted combination of the solver options runs | None |
| [`s_test_suite_regression`](regression.md) | `R(t)` matches the stored results | `tests/reference/` |
| [`s_test_suite_diffmat`](operators.md) | The differentiation operators of `f_diff_mat` | None |

All but the operator suite run every case in three solver modes:

| Mode | Solver |
|---|---|
| `nodal_fd` | `f_imr_nodal`, finite differences (`scheme = 0`) |
| `nodal_cheb` | `f_imr_nodal`, Chebyshev collocation (`scheme = 1`) |
| `spectral` | `f_imr_spectral` |

## `f_test_cases(suite)`

Returns the cases of `'equilibrium'`, `'gaussian'`, `'regression'` or `'permutations'` for all three modes. It builds each case from the model features, then keeps it only if the solver's own input checks (`f_call_params` and [`f_checkers`](../reference/functions/f_checkers.md)) accept it. The generator enumerates selected combinations and filters them through the rules behind the [compatibility matrix](../reference/compatibility.md). It does not exhaust every possible parameter value or prove that every pairwise-compatible collection of features is valid.

Each case has these fields:

- `key`: a name for the case within its mode, such as `r2_bubble+medium_s3_mat1`. Stored results are looked up by it.
- `mode`: the solver mode.
- `solver`: the solver function handle.
- `args`: the dimensionless name/value pairs, grid included.

The inputs are given in SI units and converted with [`f_dim2nondim`](../api/f_dim2nondim.md). The grids are:

| Suite | `nodal_fd` | `nodal_cheb` | `spectral` |
|---|---|---|---|
| Equilibrium | `Nt = Mt = 150` | 20 | 12 |
| Gaussian | 50 | 20 | 12 |
| Regression | 70 | 20 | 12 |
| Permutations | 20 | 12 | 8 (`Nv = 8`) |

Finite differences need finer grids because they are second-order accurate. The Chebyshev and spectral expansions converge spectrally. Except in permutations, the spectral mode sets `Nv = 32` stress modes, for the stress 5 field and the Zener field of stress 3–4 with medium heat transfer. The default is 150.

**Historically reported case counts** (`nodal_fd`, `nodal_cheb`, `spectral`):

| Suite | Cases |
|---|---|
| Equilibrium | 144, 144, 144 |
| Gaussian | 120, 120, 120 |
| Regression | 230, 230, 155 |

## Helpers

- `f_test_run(c)` runs one case with the solver's console output hidden. It returns `R`, or the error message if the run fails, so one failing case does not stop a suite.
- `f_test_report(name, cases, failed, msg)` prints the passes per mode and every failure with its mode and key. It then raises an error if any case failed. In a batch run (`matlab -batch`, as in CI) the error exits MATLAB with a nonzero code; interactively it only stops the script.

## Run the suites

From the checkout root, use a shell with licensed MATLAB on its path. Each suite starts a fresh MATLAB process. The loop preserves failure status while allowing all suites to finish.

```bash
cd toolchain
status=0
for suite in s_test_suite_equilibrium s_test_suite_gaussian \
    s_test_suite_permutations s_test_suite_diffmat s_test_suite_regression; do
    matlab -batch "$suite" || status=1
done
exit "$status"
```

Read the per-mode counts and final report. A case generator rejecting an incompatible combination is an exclusion, not a passing simulation. Counts printed by the current run are authoritative; static counts retained on these pages are historical reports, not verified counts for every MATLAB release.

The regression suite needs [approved baseline files](regression.md#baseline-availability), which are absent from this checkout. Other suites do not depend on those files. These tests exercise numerical consistency; no bundled external-paper comparison establishes scientific accuracy to a stated number of digits.
