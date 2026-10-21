# Repository map

The checkout separates solver functions, examples, developer tooling and documentation. Generated run outputs are not source files and should not be committed.

| Path | Purpose |
|---|---|
| `README.md` | Short front door and source of the local documentation home |
| `CONTRIBUTING.md` | Pointer to the contributor guide |
| `LICENSE` | GNU GPL version 3 or later |
| `src/` | MATLAB solver, parser, conversion and physics modules; [code map](code-map.md) |
| `src/d_hn.mat`, `d_mn.mat`, `d_ml.mat` | Runtime measured-waveform spline data |
| `examples/` | Simulation scripts, large sweeps and legacy postprocessing; [inventory](../getting-started/examples.md#bundled-script-inventory) |
| `toolchain/` | Seven numerical suites, helpers, formatter, linter and reference generator |
| `tests/*.mat` | Older stored data, unused by the current suites |
| `tests/reference/` | Expected regression baselines; absent from the audited checkout |
| `benchmarks/integrators/` | Reproducible integrator comparison driver; [method](../method/integrator-benchmarks.md) |
| `benchmarks/scaling/` | Legacy resolution scripts, cluster templates and stored results; [provenance](../validation/convergence.md#legacy-scaling-data) |
| `docs/` | Topic-based Markdown documentation |
| `docs/waveform/` | Raw measured data and MATLAB spline rebuild script |
| `docs/assets/workflow.svg` | Editable source of the README overview diagram; no simulation data |
| `docs/imr.jfif` | Retained original branding image; rebuild source unavailable |
| `docs/index.md` | Includes README through standard MkDocs plugin configuration, with no custom Python source |
| `docs/requirements.txt` | Documentation dependencies, independent of MATLAB |
| `docs/javascripts/mathjax.js` | MathJax configuration for equations |
| `mkdocs.yml` | Local build theme, navigation, math and strict validation |
| `.github/workflows/` | Tests, lint, format, spelling and build-only docs workflows |
| `.gitignore` | Generated/local-file patterns |
| `site/` | Ignored local documentation build output |

For local development commands, see [CI](ci.md). Some historical benchmark scripts write into their own directory; preserve supplied data and use a separate output location when rerunning them.
