<p align="center">
  <a href="https://github.com/InertialMicrocavitationRheometry/IMRv2">
    <img src="docs/imr.jfif" alt="IMR Banner" width="250"/>
  </a>
</p>

<p align="center">
  <a href="https://doi.org/10.1016/j.jmps.2017.12.006"><img src="https://zenodo.org/badge/doi/10.1016/j.jmps.2017.12.006.svg" alt="DOI"/></a>
  <a href="https://www.gnu.org/licenses/gpl-3.0.html"><img src="https://img.shields.io/badge/License-GPLv3-blue.svg" alt="License: GPLv3"/></a>
  <a href="docs/compatibility.md"><img src="https://img.shields.io/badge/Model%20compatibility-grid-2ea44f" alt="Model compatibility grid"/></a>
  <a href="CONTRIBUTING.md"><img src="https://img.shields.io/badge/Contributions-welcome-brightgreen" alt="Contributions welcome"/></a>
</p>

<p align="center">
  <a href="https://github.com/InertialMicrocavitationRheometry/IMRv2/actions/workflows/test_equilibrium.yml"><img src="https://github.com/InertialMicrocavitationRheometry/IMRv2/actions/workflows/test_equilibrium.yml/badge.svg" alt="Base tests"/></a>
  <a href="https://github.com/InertialMicrocavitationRheometry/IMRv2/actions/workflows/test_comparison.yml"><img src="https://github.com/InertialMicrocavitationRheometry/IMRv2/actions/workflows/test_comparison.yml/badge.svg" alt="Comparison tests"/></a>
  <a href="https://github.com/InertialMicrocavitationRheometry/IMRv2/actions/workflows/test_masstrans.yml"><img src="https://github.com/InertialMicrocavitationRheometry/IMRv2/actions/workflows/test_masstrans.yml/badge.svg" alt="Mass transfer tests"/></a>
  <a href="https://github.com/InertialMicrocavitationRheometry/IMRv2/actions/workflows/format.yml"><img src="https://github.com/InertialMicrocavitationRheometry/IMRv2/actions/workflows/format.yml/badge.svg" alt="Format check"/></a>
  <a href="https://github.com/InertialMicrocavitationRheometry/IMRv2/actions/workflows/spellcheck.yml"><img src="https://github.com/InertialMicrocavitationRheometry/IMRv2/actions/workflows/spellcheck.yml/badge.svg" alt="Spellcheck"/></a>
</p>

# Inertial Microcavitation Rheometry (IMR)

IMR characterizes soft materials at high strain rates. It simulates a cavitation bubble in a viscoelastic material and matches the simulated bubble radius to high-speed video of the bubble, which gives the bubble pressure and the stress field in the material.

## Contents

- [Features](#features)
- [Getting started](#getting-started)
- [Documentation](#documentation)
- [Testing](#testing)
- [Contributing](#contributing)
- [Citation](#citation)
- [Developers and contact](#developers-and-contact)
- [License and acknowledgments](#license-and-acknowledgments)

## Features

### Solvers

| Solver | Discretization | Notes |
|---|---|---|
| `f_imr_nodal`, `scheme = 0` (default) | Finite differences on node values | All models, including vapor mass transfer and collapse |
| `f_imr_nodal`, `scheme = 1` | Chebyshev collocation on node values | All models, including vapor mass transfer and collapse |
| `f_imr_spectral` | Chebyshev spectral collocation | No mass transfer or collapse; the only solver with the spectral UCM / Oldroyd-B stress field |

### Models

- **Radial dynamics:** Rayleigh–Plesset (incompressible), Keller–Miksis in pressure form (constant sound speed), Keller–Miksis in enthalpy form and Gilmore with the Tait equation of state, and the same two with the Mie–Grüneisen equation of state.
- **Heat transfer:** polytropic bubble, heat transfer in the bubble, and heat transfer in the bubble and the surrounding medium. Medium heat transfer needs radial 1 or 2, since the enthalpy forms (radial 3–6) take the liquid at the far-field temperature.
- **Stress:** neo-Hookean and quadratic neo-Hookean Kelvin–Voigt; linear Maxwell, Jeffreys and Zener; UCM and Oldroyd-B.
- **Viscosity:** Newtonian and seven non-Newtonian models (Carreau, Carreau–Yasuda, Powell–Eyring, modified Powell–Eyring, Cross, simplified Cross, modified Cross).
- **Driving pressure:** impulse, Gaussian, histotripsy, Heaviside, and three measured waveforms.
- **Time integration:** the stiff integrators `ode15s` and `ode23tb`, with an analytic Jacobian where one is implemented.

Not every combination is available. The [model compatibility grid](docs/compatibility.md) shows which features run together in each solver.

### Methods

The numerical methods are described in:
- Estrada et al., "High Strain-rate Soft Material Characterization via Inertial Cavitation," [Journal of the Mechanics and Physics of Solids, 2017](https://doi.org/10.1016/j.jmps.2017.12.006).
- Warnez and Johnsen, "Numerical modeling of bubble dynamics in viscoelastic media with relaxation," [Physics of Fluids, 2015](https://doi.org/10.1063/1.4928860).

## Getting started

### Requirements

- **MATLAB** R2021a or newer. The solvers use base MATLAB only.
- **Optional toolboxes** for work outside the solvers: the Optimization Toolbox for parameter fitting, and the Image Processing Toolbox for analyzing high-speed video.
- **GNU Octave** runs most of the solvers. It lacks `ode23tb`, the default integrator, and differs from MATLAB in a few places; see the notes at the end of [`docs/compatibility.md`](docs/compatibility.md).

### Installation

Clone the repository and add it to your MATLAB path:

```bash
git clone https://github.com/InertialMicrocavitationRheometry/IMRv2.git
cd IMRv2
```

```matlab
addpath(genpath('IMRv2'))
savepath
```

### Dimensionless inputs

`f_imr_nodal` and `f_imr_spectral` take and return **dimensionless** quantities only:

| Quantity | Scale |
|---|---|
| Length | Initial radius `R0` |
| Pressure | Far-field pressure `P8` |
| Velocity | `sqrt(P8/rho8)` |
| Time | `R0/sqrt(P8/rho8)` |
| Temperature | Far-field temperature `T8` |

The inputs and their defaults are listed in [`docs/src/default_case.md`](docs/src/default_case.md), and `f_checkers` validates every run.

### First run

Convert physical (SI) values with `toolchain/f_dim2nondim.m`, run the solver, and scale the result back:

```matlab
addpath('src'); addpath('toolchain');
[varin, ref] = f_dim2nondim('r0',50e-6,'req',5e-6,'mu',1e-2,'g',1e3, ...
    'tvector',linspace(0,12e-6,100));
[t, R] = f_imr_nodal(varin{:}, 'radial', 2, 'stress', 1);
plot(t*ref.t0, R*ref.R0);
```

`f_imr_nodal` uses finite differences by default. Add `'scheme', 1` to use Chebyshev differentiation matrices instead; the other inputs stay the same.

### Examples

- `s_run_fdmodel.m`, `s_run_full_model.m` and `s_run_mercury.m` in the repository root run example simulations.
- The `examples/` folder has more.

## Documentation

| Page | Contents |
|---|---|
| [`docs/compatibility.md`](docs/compatibility.md) | Which features run together in each solver |
| [`docs/src/default_case.md`](docs/src/default_case.md) | Every input, its dimensionless definition and its default |
| [`docs/src/f_imr_nodal.md`](docs/src/f_imr_nodal.md), [`docs/src/f_imr_spectral.md`](docs/src/f_imr_spectral.md) | The two solvers; `docs/src/` has one page per physics function |
| [`docs/toolchain/`](docs/toolchain) | The unit conversion `f_dim2nondim` and the test suites |
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | How to contribute: fork, change, check, pull request |

## Testing

The test suites live in `toolchain/` and are run from that folder:

| Suite | Checks |
|---|---|
| `s_test_suite_compatibility` | Both solvers against stored results, over radial, thermal and stress models |
| `s_test_suite_equilibrium` | A bubble at equilibrium stays there, in both solvers |
| `s_test_suite_gaussian` | A bubble returns to equilibrium after a Gaussian pressure pulse |
| `s_test_suite_masstrans` | The nodal solver with mass transfer against stored results |
| `s_test_suite_diffmat` | The differentiation operators of `f_diff_mat` |

Each suite prints whether it passed. All but `s_test_suite_masstrans` also exit with code 0 on success.

`toolchain/s_format.m` formats the code and `toolchain/s_lint.m` runs the MATLAB linter. The CI workflows in `.github/workflows` run these, the spelling check and the suites on every pull request. [`CONTRIBUTING.md`](CONTRIBUTING.md#3-run-the-checks) shows how to run them locally.

## Contributing

Contributions are welcome, from bug reports to new models. [`CONTRIBUTING.md`](CONTRIBUTING.md) covers:
- forking the repository,
- the code conventions,
- running the spelling, format, lint and test checks,
- opening a pull request.

For bugs and feature requests, open an issue on GitHub.

## Citation

If you use IMR in your research, please cite:

```bibtex
@article{Estrada_2017,
  title   = {High Strain-rate Soft Material Characterization via Inertial Cavitation},
  author  = {J. B. Estrada and C. Barajas and D. L. Henann and E. Johnsen and C. Franck},
  journal = {Journal of the Mechanics and Physics of Solids},
  year    = {2017},
  volume  = {112},
  pages   = {291–317},
  doi     = {10.1016/j.jmps.2017.12.006}
}
```

For the numerical methods, also see the two papers under [Methods](#methods).

## Developers and contact

IMR is developed by researchers at several institutions (alphabetical):
- [Spencer Bryngelson](https://comp-physics.group/) (Georgia Tech)
- [Jon Estrada](https://me.engin.umich.edu/people/faculty/jon-estrada/) (University of Michigan)
- [Christian Franck](https://directory.engr.wisc.edu/me/Faculty/Franck_Christian/) (University of Wisconsin-Madison)
- [David Henann](https://vivo.brown.edu/display/dhenann) (Brown University)
- [Eric Johnsen](https://me.engin.umich.edu/people/faculty/eric-johnsen/) (University of Michigan)
- [Mauro Rodriguez](https://vivo.brown.edu/display/mrodri97) (Brown University), **lead developer**
- [Jin Yang](https://sites.utexas.edu/yang) (University of Texas at Austin)

For questions, contact [Mauro Rodriguez](mailto:mrodri97@brown.edu) or ask to join the IMR Slack workspace.

## License and acknowledgments

IMR is released under the [GNU General Public License v3.0](LICENSE).

IMR development has been supported by the National Science Foundation (NSF), the Department of Defense (DOD), and other research institutions.
