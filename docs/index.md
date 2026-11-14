# IMR documentation

Read these Markdown pages directly on GitHub. GitHub's own renderer shows the equations and Mermaid diagrams. You do not need any local documentation tool. For installation and a first run, start with the [repository README](../README.md).

**New here?** Read these pages in order:
1. [Quick start](getting-started/quick-start.md): a complete runnable example.
2. [Assumptions](model/overview.md): what IMR models, and what it does not model.
3. [Running a simulation](running/configuration.md): how to run your own case.

Come back to the table below for all other topics. The table does not follow reading order.

| Topic | Start here |
|---|---|
| Getting started | [Installation](getting-started/installation.md), [quick start](getting-started/quick-start.md), [examples](getting-started/examples.md) |
| Running simulations | [Configuration](running/configuration.md), [solvers](running/solvers.md), [outputs](running/output.md), [sweeps](running/sweeps.md) |
| Physical models and equations | [Assumptions](model/overview.md), [radial dynamics](model/radial-dynamics.md), [stress and viscosity](model/stress-viscosity.md), [heat and mass transfer](model/heat-mass-transfer.md) |
| Numerical methods | [Workflow](method/overview.md), [discretization](method/discretization.md), [integration and Jacobians](method/time-integration.md), [wall boundaries](method/wall-boundaries.md) |
| Verification | [Test suites](verification/tests.md), [settings audit](verification/settings.md), [reported comparisons](verification/convergence.md) |
| Reference | [Public API](reference/api.md), [function call map](reference/function-call-map.md), [compatibility](reference/compatibility.md), [contributing](reference/contributing.md), [CI](reference/ci.md) |

## Finding an equation by model choice

| Choosing `radial = ` or `stress = ` | Its equation |
|---|---|
| `radial` (1–6): Rayleigh–Plesset, Keller–Miksis, Gilmore | [Radial dynamics](model/radial-dynamics.md) |
| `stress` 1–8: Kelvin-Voigt and Zener springs (neo-Hookean, quadratic, 2nd-order Fung, Mooney-Rivlin) | [Medium heating](model/stress-viscosity.md#medium-heating), with each spring's own section above it |
| `stress` 9–11: UCM/Oldroyd-B, PTT, Giesekus | [Spectral stress field](model/stress-viscosity.md#spectral-stress-field) |
| `nu_model` 1–7: non-Newtonian viscosity integral | [Non-Newtonian viscosity](model/stress-viscosity.md#non-newtonian-viscosity) |

The [stress and viscosity](model/stress-viscosity.md) page gives every model's equation by its `stress` number. The [radial dynamics](model/radial-dynamics.md) page gives every model's equation by its `radial` number. The [Method](method/overview.md) pages keep implementation details, such as equation-of-state evaluation and spectral/nodal storage. Each relevant section links to its Method page.

The MkDocs website is deferred. The [reactivation instructions](reference/ci.md#reactivating-mkdocs) keep the MkDocs configuration ready for future use.
