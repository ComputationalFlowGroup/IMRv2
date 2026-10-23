# Repository map

The checkout separates solver functions, examples, developer tooling and documentation. Generated run outputs are not source files and should not be committed.

| Path | Purpose |
|---|---|
| `README.md` | Short front door with installation, examples and documentation links |
| `CONTRIBUTING.md` | Pointer to the contributor guide |
| `LICENSE` | GNU GPL version 3 or later |
| `imr_setup.m` | One-session path setup for source and examples |
| `tests/` | [Test setup and organization](../../tests/README.md) |
| `tests/suites/` | Executable MATLAB test suites |
| `tests/helpers/` | Shared case generation, output checking and reports |
| `tests/+imrtest/` | Independent frozen data functions; not part of the solver runtime |
| `src/` | MATLAB solver, parser, conversion and physics modules; [code map](code-map.md) |
| `src/+imr/` | Public namespace and implementation packages; [organization](namespaces.md) |
| `src/+imr/+solvers/` | Shared runner and scheme-specific builders |
| `src/data/d_hn.mat`, `d_mn.mat`, `d_ml.mat` | Runtime measured-waveform spline data |
| `examples/` | Simulation scripts, large sweeps and legacy postprocessing; [inventory](../getting-started/examples.md#bundled-script-inventory) |
| `toolchain/` | Formatter, linter, profiling aids and guarded reference generator |
| `toolchain/docs/` | Bash checks for Markdown destinations and LaTeX equations; no MATLAB dependency |
| `archive/legacy-tests/` | [2,328 preserved older histories](../../archive/legacy-tests/README.md), unused by current suites |
| `tests/reference/` | Expected regression baselines; absent from the audited checkout |
| `benchmarks/integrators/` | Reproducible integrator comparison driver; [method](../method/integrator-benchmarks.md) |
| `benchmarks/scaling/` | Legacy resolution scripts, cluster templates and stored results; [provenance](../validation/convergence.md#legacy-scaling-data) |
| `docs/` | Topic-based Markdown documentation |
| `docs/waveform/` | Raw measured data and MATLAB spline rebuild script |
| `docs/assets/logos/` | Adopted TikZ logo source and compiled SVG/PDF exports |
| `docs/imr.jfif` | Retained original branding image; rebuild source unavailable |
| `docs/index.md` | GitHub-readable documentation index and future website home |
| `docs/requirements.txt` | Dormant website dependencies, independent of MATLAB |
| `mkdocs.yml` | Dormant website theme, navigation, math markup and strict validation settings |
| `.github/workflows/` | Tests, lint, format, spelling and source-documentation checks |
| `.gitignore` | Generated/local-file patterns |
| `site/` | Ignored local documentation build output |

For local development commands, see [CI](ci.md). Some historical benchmark scripts write into their own directory; preserve supplied data and use a separate output location when rerunning them.
