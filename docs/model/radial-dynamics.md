# Radial dynamics

`src/f_radial_eq.m` evaluates the wall acceleration for six choices of `radial`. The selector changes the liquid compressibility model and how pressure enters the radial equation.

| `radial` | Equation | Liquid description |
|---|---|---|
| 1 | Rayleigh–Plesset | Incompressible |
| 2 | Keller–Miksis, pressure form | Constant far-field sound speed |
| 3 | Keller–Miksis, enthalpy form | Tait EoS |
| 4 | Gilmore | Tait EoS, local sound speed |
| 5 | Keller–Miksis, enthalpy form | Mie–Grüneisen EoS |
| 6 | Gilmore | Mie–Grüneisen EoS, local sound speed |

An equation of state (EoS) relates pressure to the liquid density and acoustic quantities. The enthalpy forms use a pressure-only relation at the far-field temperature; they cannot be combined with medium heat transfer in this implementation. Radial 5–6 are further limited to the polytropic bubble. See [compatibility](../reference/compatibility.md) for complete requirements.

`cstar`, `gama`, `nstate` and `hugoniot_s` are defined in [parameters](parameters.md#acoustic-options). Choosing radial 5 or 6 does not supply a validated material parameter set. In particular, the bundled “mercury” example uses helper defaults rather than a verified mercury calibration.

The [radial theory](../theory/radial-equations.md) states implemented equations and the enthalpy construction. [Reported comparisons](../validation/convergence.md) preserve earlier numerical measurements without claiming they were rerun.
