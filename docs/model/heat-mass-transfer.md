# Heat and mass transfer

`bubtherm`, `medtherm`, `vapor` and `masstrans` choose the thermal and vapor models. The solvers always work with the dimensionless scales in [parameters](parameters.md).

| Input | Role | Implementation |
|---|---|---|
| `bubtherm = 0` | Polytropic pressure law | Solver pressure equation |
| `bubtherm = 1` | Bubble thermal PDE | `imr.solve_nodal` or `imr.solve_spectral` |
| `medtherm = 1` | Surrounding-medium thermal PDE | Solver diffusion/advection plus `imr.physics.stress_power` |
| `vapor = 1` | Saturated far-field vapor pressure | `imr.configuration.prepare_case` |
| `masstrans = 1` | Vapor diffusion and wall coupling | Nodal solver and wall closure |

The [compatibility matrix](../reference/compatibility.md) owns prerequisites and exclusions. Without medium heat transfer the wall remains at the far-field temperature; this treats the liquid as an infinite heat reservoir and neglects wall cooling/heating from phase change.

`imr.nondimensionalize` obtains saturation from the selected [fluid profile](fluids.md) through [`imr.fluid_properties`](../api/f_fluid_properties.md#saturation). Without a profile, it retains the legacy water exponential fit attributed to Preston's 2004 thesis. Its historically quoted 290–310 K interval is not a validated accuracy guarantee. Overriding `pv` alone does not replace the temperature-law exponent.

The wall law uses dimensionless temperature $\tau=T/T_\infty$ and pressure:

$$
\frac{P_{\mathrm{sat}}(T)}{P_\infty}
=p_v\exp\left[\Theta_v\left(1-\frac{1}{\tau}\right)\right],\qquad
p_v=\frac{P_{\mathrm{sat}}(T_\infty)}{P_\infty},\quad
\Theta_v=-\frac{c_2}{T_\infty}.
$$

Thus $\tau=1$ recovers the ambient vapor pressure, and negative dimensional $c_2$ gives positive `thetav`. This law is analytically inverted by the mass-transfer wall closure. A saturation table can supply an ambient value without supplying that inverse law. Such profiles must disable unsupported dynamic mass transfer and collapse.

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

This follows by writing partial vapor pressure from the ideal gas mixture, with the local vapor assumed saturated. Admissible mixtures require $0<p_{\mathrm{sat}}(\tau_w)<P$, hence $0<k_{v,w}<1$. As saturation pressure tends to zero at fixed total pressure, the wall vapor fraction tends to zero. The shared `imr.physics.saturation_pressure` evaluates pressure and derivative; the nodal factory imposes the wall fraction. The [fluid suite](../validation/tests.md) checks saturation derivatives and [system checks](../reference/internal-contracts.md#verification) check reconstructed wall saturation.

Bubble conductivity is linear in temperature. The [wall method](../method/wall-boundaries.md) explains the Kirchhoff transform, saturated vapor closure and the known stored-mixture-coefficient limitation. The [stress-heating theory](../theory/stress-heating.md) states the full-stress-power convention, including elastic heating and cooling.
