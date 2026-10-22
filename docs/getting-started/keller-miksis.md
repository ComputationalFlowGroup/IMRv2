# Worked Keller–Miksis pressure case

This example follows water properties from SI input to the pressure-form radial equation. It uses a constant ambient density and sound speed, a polytropic gas and a Newtonian liquid. It does not evaluate a changing wall-temperature EoS.

```matlab
imr_setup
[args,ref] = imr.nondimensionalize('fluid',@imr.fluids.water, ...
    'r0',50e-6,'req',50e-6,'p8',101325,'t8',298.15, ...
    'rho8',998,'c8',1480,'mu',1e-3,'surft',0.072,'g',0, ...
    'radial',2,'stress',1,'bubtherm',0,'medtherm',0, ...
    'vapor',0,'masstrans',0,'u0',0.1,'tfin',20e-6);
[t,R,U,P] = imr.solve_nodal(args{:});
plot(t*ref.t0*1e6,R*ref.R0*1e6)
xlabel('Time (microseconds)')
ylabel('Radius (micrometers)')
```

The explicit `rho8` and `c8` override the legacy water reference values for this example; they are not a new water fit. The zero shear modulus gives `ca = Inf`. `stress = 1` then supplies Newtonian viscosity without elasticity. `Req = R0` sets the initial mechanical equilibrium pressure; the nonzero initial velocity perturbs that equilibrium and starts the oscillation.

## Conversion and interpretation

For ambient pressure $P_\infty$, density $\rho_\infty$ and initial radius $R_0$, define

$$
U_c=\sqrt{P_\infty/\rho_\infty},\qquad
t_0=R_0/U_c,\qquad
C=C_\infty/U_c,\qquad
\mathrm{Re}=\frac{P_\infty R_0}{\mu U_c},\qquad
\mathrm{We}=\frac{P_\infty R_0}{2\sigma}.
$$

Here $\mu$ is dynamic viscosity in Pa s and $\sigma$ is surface tension in N/m. All solver arguments and returned quantities are dimensionless; physical time, radius, velocity and pressure are `t*ref.t0`, `R*ref.R0`, `U*ref.Uc` and `P*ref.P8`.

| Converted quantity | Approximate value | Where it enters |
|---|---|---|
| `ref.Uc` | 10.0761 m/s | Velocity and acoustic scaling |
| `ref.t0` | 4.96223 microseconds | Time scaling |
| `cstar` | 146.882 | Finite sound-speed terms in `imr.physics.radial_acceleration`, radial 2 |
| `re` | 502.798 | Stress $S=-4\dot R/(\mathrm{Re}R)$ and `JdotA = 4/re` |
| `we` | 35.1823 | Laplace pressure $1/(\mathrm{We}R)$ |
| `req`, `ca` | 1, `Inf` | Equilibrium radius and zero elasticity |
| `tfin` | 4.03045 | End of the dimensionless time interval |

`imr.nondimensionalize` computes these groups once. `imr.configuration.prepare_case` derives initial pressure and acceleration-dependent stress coefficients. The nodal factory passes them by name into the RHS, which calls `imr.physics.material_stress` and `imr.physics.radial_acceleration`. The [pressure-form equation](../theory/radial-equations.md#pressure-forms) defines every term and its sign.

The equation requires positive radius and a nonsingular acceleration denominator. The Keller–Miksis approximation describes finite sound-speed effects through first-order wall-Mach terms; it is not a model for arbitrary supersonic wall motion. Increasing `cstar` should approach the [Rayleigh–Plesset limit](../theory/radial-equations.md#pressure-forms). See [numerical checks](../validation/numerics.md) and [system checks](../reference/internal-contracts.md#verification) for executable consistency tests rather than a claimed experimental accuracy.
