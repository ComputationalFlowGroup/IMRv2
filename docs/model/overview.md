# Problem and assumptions

IMR models a spherical cavitation bubble in a surrounding material. Researchers use its radius history to study material response at high strain rates and compare simulations with high-speed bubble measurements.

The forward solvers compute radial motion coupled to selected heat, vapor and stress models. Parameter fitting and video analysis are separate workflows; this repository's documented entry points are forward solvers and parameter sweeps.

| Part | Choices | Canonical detail |
|---|---|---|
| Radial dynamics | Rayleigh–Plesset, Keller–Miksis and Gilmore variants | [Radial models](radial-dynamics.md) |
| Constitutive response | Kelvin–Voigt, Maxwell/Jeffreys/Zener, UCM/Oldroyd-B | [Stress and viscosity](stress-viscosity.md) |
| Transport | Polytropic bubble, heat and vapor mass transfer | [Heat and mass transfer](heat-mass-transfer.md) |
| Forcing | Constant, Gaussian, histotripsy, Heaviside and measured waveforms | [Waveforms](../running/waveforms.md) |
| Discretization | Nodal finite differences, nodal Chebyshev or spectral expansions | [Methods](../method/discretization.md) |

The spherical formulation does not resolve nonspherical shapes or jets. Model selection, admissible pressure/temperature ranges, spatial resolution and time tolerances limit interpretation of results. A successful integration is not evidence that the chosen physics is appropriate. Use [compatibility](../reference/compatibility.md) before combining features and [validation](../validation/tests.md) to distinguish executable checks from reported measurements.

For mercury-specific ambient properties and material limits, see the [mercury model](mercury.md).

The original methods are described by Estrada et al. (2017) and Warnez and Johnsen (2015); see [references](../reference/references.md). Original project attribution and contacts remain in [about](../reference/about.md).
