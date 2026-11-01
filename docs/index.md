# IMR documentation

Read these Markdown pages directly on GitHub. Equations and Mermaid diagrams use GitHub's renderer; no local documentation tools are needed. For installation and a first run, start with the [repository README](../README.md).

**New here?** Read, in order: the [quick start](getting-started/quick-start.md) for a complete runnable example, [assumptions](model/overview.md) for what IMR models and does not, then [running a simulation](running/configuration.md) for your own case. Come back to the table below for everything else, by topic rather than by reading order.

| Topic | Start here |
|---|---|
| Getting started | [Installation](getting-started/installation.md), [quick start](getting-started/quick-start.md), [examples](getting-started/examples.md) |
| Running simulations | [Configuration](running/configuration.md), [solvers](running/solvers.md), [outputs](running/output.md), [sweeps](running/sweeps.md) |
| MATLAB API | [Public functions](api/index.md), [unit conversion](api/nondimensionalize.md), [extension guide](api/extending.md) |
| Physical models | [Assumptions](model/overview.md), [fluid profiles](model/fluids.md), [stress and viscosity](model/stress-viscosity.md), [heat and mass transfer](model/heat-mass-transfer.md) |
| Numerical methods | [Workflow](method/overview.md), [discretization](method/discretization.md), [integration and Jacobians](method/time-integration.md), [wall boundaries](method/wall-boundaries.md) |
| Equations | [Radial dynamics](theory/radial-equations.md), [stress heating](theory/stress-heating.md), [viscosity integrals](theory/viscosity-integrals.md) |
| Validation | [Test suites](validation/tests.md), [settings audit](validation/settings.md), [reported comparisons](validation/convergence.md) |
| Development | [Function call map](reference/function-call-map.md), [compatibility](reference/compatibility.md), [contributing](reference/contributing.md), [CI](reference/ci.md) |

## Finding an equation by model choice

| Choosing `radial = ` or `stress = ` | Its equation |
|---|---|
| `radial` (1–6): Rayleigh–Plesset, Keller–Miksis, Gilmore | [Radial equations and equation of state](theory/radial-equations.md) |
| `stress` 1–8: Kelvin-Voigt and Zener springs (neo-Hookean, quadratic, 2nd-order Fung, Mooney-Rivlin) | [Heating per model](theory/stress-heating.md#heating-per-model), with the elastic stress difference just below it |
| `stress` 9–11: UCM/Oldroyd-B, PTT, Giesekus | [Spectral stress field](theory/stress-heating.md#spectral-stress-field) |
| `nu_model` 1–7: non-Newtonian viscosity integral | [Viscosity integrals](theory/viscosity-integrals.md) |

The [stress and viscosity](model/stress-viscosity.md) page lists every model by its `stress` number with the same structure; the [radial dynamics](model/radial-dynamics.md) page does the same for `radial`. Each links back to the theory page above for the derivation.

The MkDocs website is deferred. Its configuration is retained for [future activation](reference/ci.md#reactivating-mkdocs).
