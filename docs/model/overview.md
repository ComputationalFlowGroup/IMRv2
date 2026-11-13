# Problem and assumptions

IMR models a spherical cavitation bubble in a surrounding material. Researchers use the bubble radius history to study material response at high strain rates. Researchers can compare the simulation results with high-speed bubble measurements.

The forward solvers compute radial motion. The radial motion is coupled to selected heat, vapor and stress models. Parameter fitting and video analysis are separate workflows. This repository's documented entry points are forward solvers and parameter sweeps.

| Part | Choices | Canonical detail |
|---|---|---|
| Radial dynamics | Rayleigh–Plesset, Keller–Miksis and Gilmore variants | [Radial models](radial-dynamics.md) |
| Constitutive response | Kelvin–Voigt, Maxwell/Jeffreys/Zener, UCM/Oldroyd-B, PTT, Giesekus | [Stress and viscosity](stress-viscosity.md) |
| Transport | Polytropic bubble, heat and vapor mass transfer | [Heat and mass transfer](heat-mass-transfer.md) |
| Forcing | Constant, Gaussian, histotripsy, Heaviside and measured waveforms | [Waveforms](../running/waveforms.md) |
| Discretization | Nodal finite differences, nodal Chebyshev or spectral expansions | [Methods](../method/discretization.md) |

The spherical formulation does not resolve nonspherical shapes or jets. Model selection, pressure/temperature ranges, spatial resolution and time tolerances all limit how you can interpret the results. A successful integration run is not proof that the chosen physics model is correct for the case. Check [compatibility](../reference/compatibility.md) before you combine features. Check [verification](../verification/tests.md) to see which checks the code runs and which results are only reported measurements.

Estrada et al. (2017) and Warnez and Johnsen (2015) describe the original methods; see [references](../reference/references.md). The [about](../reference/about.md) page keeps the original project attribution and contacts.
