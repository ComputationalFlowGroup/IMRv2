# MATLAB and batch commands

IMR exposes MATLAB functions, not an application command-line interface. Use MATLAB interactively or its `-batch` option from a shell.

All commands below start from the checkout root unless a directory is stated. MATLAB commands have not been executed for this documentation change.

| Task | Command | Working directory |
|---|---|---|
| Minimal solver call | `matlab -batch "imr_setup; [t,R]=imr.solve_nodal('radial',2,'stress',1);"` | Root |
| Function help | `matlab -batch "imr_setup; help imr.solve_nodal"` | Root |
| Bundled Gaussian case | `matlab -batch "run('examples/s_run_full_model.m')"` | Root |
| Operator tests | `matlab -batch "imr_test_setup; s_test_suite_diffmat"` | `tests/` |
| All numerical suites | `matlab -batch "run('tests/s_run_validation.m')"` | Root |
| Selected suites | `matlab -batch "run_validation('tests',{'diffmat','solver'})"` | `tests/` |
| Generate golden data | `matlab -batch "run_validation('generate')"` | `tests/` |
| Precheck, generate and compare | `matlab -batch "run_validation('generate-and-test')"` | `tests/` |
| Numerical test commands | See [tests](../validation/tests.md#run-the-suites) | `tests/` |
| Source documentation checks | `bash toolchain/docs/check.sh` | Root |
| Documentation checker failure tests | `bash toolchain/docs/test.sh` | Root |
| Spelling, lint, formatting | See [CI](../reference/ci.md) | As listed there |

| Option | Consumer | Meaning |
|---|---|---|
| `-batch` | MATLAB executable | Evaluates MATLAB code noninteractively and returns a failure status for an uncaught error |

Model choices such as `radial` and `scheme` are MATLAB name/value inputs, not shell flags. See [configuration](configuration.md) and [parameters](../model/parameters.md).

The website is dormant. Future website commands are described under [reactivating MkDocs](../reference/ci.md#reactivating-mkdocs).
