# Test suites and `f_test_cases.m`

## Overview

The test suites live in `toolchain/` and run from that folder. The tests workflow `tests.yml` runs each one in its own MATLAB process on every pull request, and fails if any suite fails.

| Suite | Checks | Stored results |
|---|---|---|
| [`s_test_suite_equilibrium`](s_test_suite_equilibrium.md) | A bubble at equilibrium stays there | None |
| [`s_test_suite_gaussian`](s_test_suite_gaussian.md) | A bubble returns to equilibrium after a Gaussian pressure pulse | None |
| [`s_test_suite_regression`](s_test_suite_regression.md) | `R(t)` matches the stored results | `tests/reference/` |
| [`s_test_suite_diffmat`](s_test_suite_diffmat.md) | The differentiation operators of `f_diff_mat` | None |

The first three run every case in three solver modes:

| Mode | Solver |
|---|---|
| `nodal_fd` | `f_imr_nodal`, finite differences (`scheme = 0`) |
| `nodal_cheb` | `f_imr_nodal`, Chebyshev collocation (`scheme = 1`) |
| `spectral` | `f_imr_spectral` |

## `f_test_cases(suite)`

Returns the cases of `'equilibrium'`, `'gaussian'` or `'regression'` for all three modes. It builds each case from the model features, then keeps it only if the solver's own input checks (`f_call_params` and [`f_checkers`](../src/f_checkers.md)) accept it. So the suites cover exactly the combinations of the [compatibility matrix](../compatibility.md), and they follow it when the checks change.

Each case has these fields:

- `key`: a name for the case within its mode, such as `r2_bubble+medium_s3_mat1`. Stored results are looked up by it.
- `mode`: the solver mode.
- `solver`: the solver function handle.
- `args`: the dimensionless name/value pairs, grid included.

The inputs are given in SI units and converted with [`f_dim2nondim`](../src/f_dim2nondim.md). The grids are:

| Suite | `nodal_fd` | `nodal_cheb` | `spectral` |
|---|---|---|---|
| Equilibrium | `Nt = Mt = 150` | 20 | 12 |
| Gaussian | 50 | 20 | 12 |
| Regression | 70 | 20 | 12 |

Finite differences need finer grids because they are second-order accurate. The Chebyshev and spectral expansions converge spectrally.

**Case counts** (`nodal_fd`, `nodal_cheb`, `spectral`):

| Suite | Cases |
|---|---|
| Equilibrium | 144, 144, 144 |
| Gaussian | 120, 120, 120 |
| Regression | 230, 230, 155 |

## Helpers

- `f_test_run(c)` runs one case with the solver's console output hidden. It returns `R`, or the error message if the run fails, so one failing case does not stop a suite.
- `f_test_report(name, cases, failed, msg)` prints the passes per mode and every failure with its mode and key. It then raises an error if any case failed. In a batch run (`matlab -batch`, as in CI) the error exits MATLAB with a nonzero code; interactively it only stops the script.
