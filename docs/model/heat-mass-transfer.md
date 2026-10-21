# Heat and mass transfer

`bubtherm`, `medtherm`, `vapor` and `masstrans` choose the thermal and vapor models. The solvers always work with the dimensionless scales in [parameters](parameters.md).

| Input | Role | Implementation |
|---|---|---|
| `bubtherm = 0` | Polytropic pressure law | Solver pressure equation |
| `bubtherm = 1` | Bubble thermal PDE | `f_imr_nodal` or `f_imr_spectral` |
| `medtherm = 1` | Surrounding-medium thermal PDE | Solver diffusion/advection plus `f_stress_dissipation` |
| `vapor = 1` | Saturated far-field vapor pressure | `f_call_params` |
| `masstrans = 1` | Vapor diffusion and wall coupling | Nodal solver and wall closure |

The [compatibility matrix](../reference/compatibility.md) owns prerequisites and exclusions. Without medium heat transfer the wall remains at the far-field temperature; this treats the liquid as an infinite heat reservoir and neglects wall cooling/heating from phase change.

`f_dim2nondim` obtains saturation from the selected [fluid profile](fluids.md) through [`f_fluid_properties`](../api/f_fluid_properties.md#saturation). Without a profile, it retains the legacy water exponential fit attributed to Preston's 2004 thesis. Its historically quoted 290–310 K interval is not a validated accuracy guarantee. Overriding `pv` alone does not replace the temperature-law exponent.

The wall law uses dimensionless temperature $\tau=T/T_\infty$ and pressure:

$$
\frac{P_{\mathrm{sat}}(T)}{P_\infty}
=p_v\exp\left[\Theta_v\left(1-\frac{1}{\tau}\right)\right],\qquad
p_v=\frac{P_{\mathrm{sat}}(T_\infty)}{P_\infty},\quad
\Theta_v=-\frac{c_2}{T_\infty}.
$$

Thus $\tau=1$ recovers the ambient vapor pressure, and negative dimensional $c_2$ gives positive `thetav`. This law is analytically inverted by the mass-transfer wall closure. A saturation table can supply an ambient value without supplying that inverse law. Such profiles must disable unsupported dynamic mass transfer and collapse.

Bubble conductivity is linear in temperature. The [wall method](../method/wall-boundaries.md) explains the Kirchhoff transform, saturated vapor closure and the known stored-mixture-coefficient limitation. The [stress-heating theory](../theory/stress-heating.md) states the full-stress-power convention, including elastic heating and cooling.
