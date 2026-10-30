# References and citation

Cite the source papers for the models and measured waveforms you use. The project citation is given below.

## Project citation

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

## Numerical methods

- Estrada et al., “High Strain-rate Soft Material Characterization via Inertial Cavitation,” *Journal of the Mechanics and Physics of Solids* (2017). <https://doi.org/10.1016/j.jmps.2017.12.006>
- Warnez and Johnsen, “Numerical modeling of bubble dynamics in viscoelastic media with relaxation,” *Physics of Fluids* 27, 063103 (2015). <https://doi.org/10.1063/1.4922598>. Their Eq. 13 is the PTT and Giesekus closure of the spectral stress field, `stress = 6` and `7` (see [stress and viscosity](../model/stress-viscosity.md#ptt-and-giesekus-stress-6-7)).
- Yang and Church, “A model for the dynamics of gas bubbles in soft tissue,” *Journal of the Acoustical Society of America* 118, 3595–3606 (2005). <https://doi.org/10.1121/1.2118307>. Their linear Kelvin-Voigt spring is `stress = 8` (see [stress and viscosity](../model/stress-viscosity.md#linear-kelvin-voigt-stress-8)).

The original project is [InertialMicrocavitationRheometry/IMRv2](https://github.com/InertialMicrocavitationRheometry/IMRv2). This checkout and its contribution instructions refer to the `mrodrig6/imrv2` fork.

## Measured waveforms

If you use the measured pressure waveforms (`wave_type = -1, -2, -3`), please also cite the paper they come from:

```bibtex
@article{Hersey_2023,
  title   = {Dynamics of an oscillating microbubble in a blood-like Carreau fluid},
  author  = {E. Hersey and M. Rodriguez and E. Johnsen},
  journal = {Journal of the Acoustical Society of America},
  year    = {2023},
  volume  = {153},
  number  = {3},
  pages   = {1836–1845},
  doi     = {10.1121/10.0017342}
}
```
