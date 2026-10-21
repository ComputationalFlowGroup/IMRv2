# Radial equations and equation of state

`src/f_radial_eq.m` computes `Rddot` using the selected radial model. These are dimensionless equations; the [parameter table](../model/parameters.md) defines their scales and names.

## Rayleigh–Plesset and pressure-form Keller–Miksis

Write `Q = P - 1 - Pf8 - iWe/R + S`. The Rayleigh–Plesset branch is

$$\ddot R = \frac{Q - 3\dot R^2/2}{R}.$$

The pressure-form Keller–Miksis branch adds sound-speed terms to the numerator and denominator:

$$\ddot R = \frac{(1+\dot R/C)Q + (R/C)(\dot P+iWe\dot R/R^2+\dot S-\dot P_{f8})
-\tfrac32(1-\dot R/(3C))\dot R^2}
{(1-\dot R/C)R+JdotA/C+6R\,ddintfnu\,iDRe/C}.$$

`Pf8` is the far-field pressure above ambient. The `iDRe` definition, including the solvent share for stress 3–5, lives in [stress and viscosity](../model/stress-viscosity.md). The stress derivative is split so that its acceleration coefficient moves to the denominator.

## Enthalpy forms

Radial 3–6 replace the pressure difference with liquid enthalpy and use `hH = 1/rho` for its pressure derivative. Tait branches use the pressure shift `GAMa`, with local sound speed for Gilmore and far-field sound speed for Keller–Miksis. Mie–Grüneisen branches compute density strain from a stable quadratic root.

For `s = hugoniot_s`, the implemented Mie–Grüneisen relation is

$$A = p/Cstar^2 = g(\mu)=\frac{\mu(1+nog\,\mu)}{(1-s\mu)^2},\qquad \rho=1+\mu.$$

The enthalpy difference is formed by integration by parts:

$$h(P)=\int_1^P\frac{dp}{1+\mu}
=\left[\frac{p}{1+\mu}\right]_1^P+Cstar^2\int_{\mu(1)}^{\mu(P)}\frac{g(\mu)}{(1+\mu)^2}\,d\mu.$$

Partial fractions at `mu = -1` and `mu = 1/s` give the closed form in the nested `f_mie_gruneisen_eos_scalar` and `f_mie_mu`. The differences use stable logarithms near the reference state. This is the implemented pressure-based liquid model, not a full arbitrary-temperature thermodynamic EoS.

See [compatibility](../reference/compatibility.md) for thermal restrictions and [references](../reference/references.md) for the published physical background. Historical numerical agreement of this enthalpy formula is preserved in [reported comparisons](../validation/convergence.md#reported-radial-enthalpy-comparison).
