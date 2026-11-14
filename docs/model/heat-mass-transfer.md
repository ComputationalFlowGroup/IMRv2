# Heat and mass transfer

`bubtherm`, `medtherm`, `vapor` and `masstrans` choose the thermal and vapor models. The solvers always work with the dimensionless scales in [parameters](parameters.md).

| Input | Role | Implementation |
|---|---|---|
| `bubtherm = 0` | Polytropic pressure law | Solver pressure equation |
| `bubtherm = 1` | Bubble thermal PDE | `imr.solve_nodal` or `imr.solve_spectral` |
| `medtherm = 1` | Surrounding-medium thermal PDE | Solver diffusion/advection plus `imr.physics.stress_power` |
| `vapor = 1` | Saturated far-field vapor pressure | `imr.configuration.prepare_case` |
| `masstrans = 1` | Vapor diffusion and wall coupling | Nodal solver and wall closure |

The [compatibility matrix](../reference/compatibility.md) owns the prerequisites and exclusions. Without medium heat transfer, the wall stays at the far-field temperature. This treats the liquid as an infinite heat reservoir. This also neglects wall cooling and heating from phase change.

`imr.nondimensionalize` obtains saturation from the selected [fluid profile](fluids.md) through [`imr.fluid_properties`](../reference/fluid_properties.md#saturation). Without a profile, it uses the legacy water exponential fit from Preston's 2004 thesis. The historically quoted 290–310 K interval is not a validated accuracy guarantee. If you override `pv` alone, you do not replace the temperature-law exponent.

The wall law uses dimensionless temperature $\tau=T/T_\infty$ and pressure:

$$
\frac{P_{\mathrm{sat}}(T)}{P_\infty}
=p_v\exp\left[\Theta_v\left(1-\frac{1}{\tau}\right)\right],\qquad
p_v=\frac{P_{\mathrm{sat}}(T_\infty)}{P_\infty},\quad
\Theta_v=-\frac{c_2}{T_\infty}.
$$

Thus $\tau=1$ recovers the ambient vapor pressure. A negative dimensional $c_2$ gives a positive `thetav`. The mass-transfer wall closure analytically inverts this law. A saturation table can supply an ambient value without supplying that inverse law. A profile like this must disable dynamic mass transfer and must disable collapse. The solvers do not support either feature without the inverse law.

| Symbol | Meaning | Scale / domain |
|---|---|---|
| *T*, *T∞* | Physical and ambient absolute temperatures | K, both positive |
| *τ* | Temperature passed to the dynamic wall law | *T*/*T∞* > 0 |
| *p_v* | Enabled ambient saturation pressure | `Pv_star`, pressure divided by *P∞* |
| *Θ_v* | Inverse-temperature exponent | `thetav`, positive for active mass transfer |
| *k_v* | Vapor mass fraction | Between zero and one |
| *r* | Vapor to noncondensable gas-constant ratio | `rv/ra`, dimensionless and positive |

For total bubble pressure $P>0$, the ideal-mixture wall fraction is

$$
k_{v,w}=\left[1+r\left(\frac{P}{p_{\mathrm{sat}}(\tau_w)}-1\right)\right]^{-1},
\qquad
\frac{dp_{\mathrm{sat}}}{d\tau}
=\frac{\Theta_v}{\tau^2}p_{\mathrm{sat}}(\tau).
$$

This result follows from the partial vapor pressure of the ideal gas mixture. It assumes the local vapor is saturated. Admissible mixtures require $0<p_{\mathrm{sat}}(\tau_w)<P$. This gives $0<k_{v,w}<1$. As the saturation pressure tends to zero at fixed total pressure, the wall vapor fraction also tends to zero. The shared `imr.physics.saturation_pressure` function evaluates the pressure and its derivative. The nodal factory imposes the wall fraction. The [fluid suite](../verification/tests.md) checks the saturation derivatives. [System checks](../reference/internal-contracts.md#verification) check the reconstructed wall saturation.

Bubble conductivity is linear in temperature. The [wall method](../method/wall-boundaries.md) explains the Kirchhoff transform, saturated vapor closure and the known stored-mixture-coefficient limitation. The [medium heating](stress-viscosity.md#medium-heating) section states the full-stress-power convention, including elastic heating and cooling.
