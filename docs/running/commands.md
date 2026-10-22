# MATLAB and batch commands

IMR exposes MATLAB functions, not an application command-line interface. Use MATLAB interactively or its `-batch` option from a shell.

All commands below start from the checkout root unless a directory is stated. MATLAB commands have not been executed for this documentation change.

| Task | Command | Working directory |
|---|---|---|
| Minimal solver call | `matlab -batch "imr_setup; [t,R]=imr.solve_nodal('radial',2,'stress',1);"` | Root |
| Function help | `matlab -batch "imr_setup; help imr.solve_nodal"` | Root |
| Bundled Gaussian case | `matlab -batch "run('examples/s_run_full_model.m')"` | Root |
| Operator tests | `matlab -batch "s_test_suite_diffmat"` | `toolchain/` |
| Numerical test commands | See [tests](../validation/tests.md#run-the-suites) | `toolchain/` |
| Strict documentation build | `mkdocs build --strict` | Root |
| Local documentation server | `mkdocs serve` | Root |
| Spelling, lint, formatting | See [CI](../reference/ci.md) | As listed there |

| Option | Consumer | Meaning |
|---|---|---|
| `-batch` | MATLAB executable | Evaluates MATLAB code noninteractively and returns a failure status for an uncaught error |
| `--strict` | MkDocs | Makes documentation warnings fail the build |
| `--site-dir` | `mkdocs build` | Overrides the generated documentation output directory |

Model choices such as `radial` and `scheme` are MATLAB name/value inputs, not shell flags. See [configuration](configuration.md) and [parameters](../model/parameters.md).
