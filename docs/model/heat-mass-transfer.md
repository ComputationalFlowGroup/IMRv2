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

`f_dim2nondim` uses the water-vapor fit in `f_pvsat`, attributed to Preston's 2004 thesis. The former function page stated a typical biological-temperature interval of 290–310 K, but the code neither checks it nor supplies a validation dataset. That interval remains an unverified historical statement, not an established accuracy guarantee. Overriding `pv` does not replace the temperature-law exponent; supply consistent parameters for another material.

Bubble conductivity is linear in temperature. The [wall method](../method/wall-boundaries.md) explains the Kirchhoff transform, saturated vapor closure and the known stored-mixture-coefficient limitation. The [stress-heating theory](../theory/stress-heating.md) states the full-stress-power convention, including elastic heating and cooling.
