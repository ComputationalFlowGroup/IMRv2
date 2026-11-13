# Installation

IMRv2 runs in MATLAB R2021a or a later release. The forward solvers need only base MATLAB. You can read this documentation directly on GitHub. You do not need to build a website.

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

Each `which` command must resolve to a path inside this checkout's `src/` folder. `imr_setup` also adds the example scripts and the root setup function to the path. These stay callable after you change folders. `imr_setup` does not use `genpath`. It does not add test fixtures. It does not save the path permanently. Run `imr_setup` once in each MATLAB session. You can also save the path with MATLAB's `savepath` command.

Resolving a function name is not a full check. Run the [quick start](quick-start.md) to check that the solver runs.

## Documentation tools

You do not need to install anything to read the [documentation index](../index.md). To check the source documentation, first install the tools listed in [CI](../reference/ci.md#source-documentation-checks). Then run this command from the checkout root:

```bash
bash toolchain/docs/check.sh
```

The website configuration stays in the repository, but it is not active. See [reactivating MkDocs](../reference/ci.md#reactivating-mkdocs) if you want to build the website. See [troubleshooting](../running/troubleshooting.md) for solver setup errors.
