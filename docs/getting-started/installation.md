# Installation

IMRv2 runs in MATLAB R2021a or newer. The forward solvers use base MATLAB; Python is needed only for documentation and spelling tools.

| Component | Requirement | Scope |
|---|---|---|
| MATLAB | R2021a+ | Forward solvers and numerical tests |
| MATLAB `ode` with `cvodesstiff` | R2024a+ | Optional CVODES integration |
| Parallel Computing Toolbox | Optional | Parallel sweeps and parallel reference generation |
| Optimization Toolbox | Optional | Parameter fitting outside the forward solvers |
| Image Processing Toolbox | Optional | Video analysis outside the forward solvers |
| Curve Fitting Toolbox | Optional | Rebuilding measured waveform derivatives with `fnder` |
| Symbolic Math Toolbox | Check `which hypergeom` | Standalone approximate initial-stress helper, not the solver path |
| Python | 3.12 used by documentation CI | MkDocs and codespell only |

## Get the code

From a shell:

```bash
git clone https://github.com/mrodrig6/imrv2.git
cd imrv2
```

In MATLAB, from the checkout root:

```matlab
imr_setup;
which f_imr_nodal
which f_imr_spectral
```

`which` must resolve both functions inside this checkout's `src/`. `imr_setup` also exposes the example scripts and the root setup function, so these remain callable after changing folders. It does not use `genpath`, add test fixtures or save paths permanently. Run it once per session, or persist the resulting paths with MATLAB's `savepath` if desired. The [quick start](quick-start.md) exercises the solver; resolving a function alone does not validate a simulation.

## Documentation tools

Use a virtual environment from the checkout root:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r docs/requirements.txt codespell
mkdocs build --strict
```

The generated `site/` directory is local output. CI builds documentation without deploying it. See [CI](../reference/ci.md) for all check commands and [troubleshooting](../running/troubleshooting.md) for common errors.
