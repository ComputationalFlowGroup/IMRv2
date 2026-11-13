# Inertial Microcavitation Rheometry (IMR)

[![Tests](https://github.com/ComputationalFlowGroup/IMRv2/actions/workflows/tests.yml/badge.svg)](https://github.com/ComputationalFlowGroup/IMRv2/actions/workflows/tests.yml)
[![Docs checks](https://github.com/ComputationalFlowGroup/IMRv2/actions/workflows/docs.yml/badge.svg)](https://github.com/ComputationalFlowGroup/IMRv2/actions/workflows/docs.yml)
[![Lint](https://github.com/ComputationalFlowGroup/IMRv2/actions/workflows/lint.yml/badge.svg)](https://github.com/ComputationalFlowGroup/IMRv2/actions/workflows/lint.yml)
[![Formatting](https://github.com/ComputationalFlowGroup/IMRv2/actions/workflows/format.yml/badge.svg)](https://github.com/ComputationalFlowGroup/IMRv2/actions/workflows/format.yml)
[![Spelling](https://github.com/ComputationalFlowGroup/IMRv2/actions/workflows/spellcheck.yml/badge.svg)](https://github.com/ComputationalFlowGroup/IMRv2/actions/workflows/spellcheck.yml)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/License-GPL--3.0--or--later-blue.svg)](https://github.com/ComputationalFlowGroup/IMRv2/blob/main/LICENSE)

![IMR monogram with an open bubble ring. Inertial Microcavitation Rheometry.](docs/assets/logos/imr-monogram.svg)

The Inertial Microcavitation Technique (IMR) forward solver simulates a spherical cavitation bubble in a viscoelastic material. Researchers use the bubble-radius history to study soft materials at high strain rates. Researchers can compare the simulation results with high-speed measurements.

IMR gives you:
- Bubble radius, wall velocity, optional wall acceleration and pressure, with Rayleigh–Plesset, Keller–Miksis or Gilmore dynamics.
- Bubble and medium temperatures, with optional vapor mass transfer in the nodal solver.
- Viscoelastic stress models and Newtonian or shear-dependent viscosity.
- A water ambient fluid, with an extensible profile structure. You can add a custom fluid with the data schema.
- Parameter sweeps, measured driving waveforms and numerical test suites.

You can extend IMR. Add inputs or constitutive equations with the [documented MATLAB extension path](docs/reference/extending.md).

## Install

Use MATLAB R2021a or a later release. The forward solvers need only base MATLAB. Clone this repository. Run `imr_setup` once in each MATLAB session. See [installation](docs/getting-started/installation.md) for the optional toolboxes and documentation tools.

## Quick start

```matlab
imr_setup;
p = imr.default_inputs();
p.radial = 2;
p.stress = 1;
[t, R] = imr.solve(p,'solver','nodal');
```

See the [complete quick start](docs/getting-started/quick-start.md) for SI conversion and plotting. Nobody has run these examples on local MATLAB yet, so this page does not claim any numerical output. Set `'solver','spectral'` to select the spectral solver. See the [unified solver API](docs/reference/solve.md).

## Verification

```matlab
run('tests/s_run_verification.m');
```

This command runs all available MATLAB test suites from the checkout root. See the [verification runner](tests/README.md#one-command-for-all-tests) to select suites, to generate golden data, or to generate golden data and then run the regression suite. The default test mode never generates or overwrites a baseline.

## Documentation

Read the [documentation index](docs/index.md) directly on GitHub for the complete page list, organized by topic. The MkDocs website is deferred. [Reactivation instructions](docs/reference/ci.md#reactivating-mkdocs) are available for later use.

IMR builds on the inertial-cavitation methods of Estrada et al. (2017) and Warnez and Johnsen (2015). See [references and citation](docs/reference/references.md) for the project BibTeX entry and the papers to cite for measured waveforms.

## License

[GNU GPL version 3 or later](https://github.com/mrodrig6/imrv2/blob/main/LICENSE), copyright the IMR developers.
