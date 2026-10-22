# Installation

IMRv2 runs in MATLAB R2021a or newer. The forward solvers use base MATLAB. Read documentation directly on GitHub; no website build is needed.

| Component | Requirement | Scope |
|---|---|---|
| MATLAB | R2021a+ | Forward solvers and numerical tests |
| MATLAB `ode` with `cvodesstiff` | R2024a+ | Optional CVODES integration |
| Parallel Computing Toolbox | Optional | Parallel sweeps and parallel reference generation |
| Optimization Toolbox | Optional | Parameter fitting outside the forward solvers |
| Image Processing Toolbox | Optional | Video analysis outside the forward solvers |
| Curve Fitting Toolbox | Optional | Rebuilding measured waveform derivatives with `fnder` |
| Symbolic Math Toolbox | Check `which hypergeom` | Standalone approximate initial-stress helper, not the solver path |
| Pandoc, jq and pdfLaTeX | Optional developer tools | Markdown destinations and equation syntax checks |
| Python | Optional developer tool | Codespell; future MkDocs activation |

## Get the code

From a shell:

```bash
git clone https://github.com/mrodrig6/imrv2.git
cd imrv2
```

In MATLAB, from the checkout root:

```matlab
imr_setup;
which imr.solve_nodal
which imr.solve_spectral
```

`which` must resolve both functions inside this checkout's `src/`. `imr_setup` also exposes the example scripts and the root setup function, so these remain callable after changing folders. It does not use `genpath`, add test fixtures or save paths permanently. Run it once per session, or persist the resulting paths with MATLAB's `savepath` if desired. The [quick start](quick-start.md) exercises the solver; resolving a function alone does not validate a simulation.

## Documentation tools

Reading the [documentation index](../index.md) requires no installation. For source-documentation checks, install the tools listed in [CI](../reference/ci.md#source-documentation-checks), then run from the checkout root:

```bash
bash toolchain/docs/check.sh
```

Website configuration is retained but inactive. See [reactivating MkDocs](../reference/ci.md#reactivating-mkdocs) when a website is wanted, and [troubleshooting](../running/troubleshooting.md) for solver setup errors.
