# Stress power and local fields

`imr.physics.stress_power` uses the full stress power, including elastic work, for medium heating.
The following derivation states that convention and the local Zener field chosen to reproduce the stress integral.

`imr.physics.stress_power` computes the medium heating by the stress, `taugradu = Br*Phi`, when `medtherm = 1`. The medium temperature equation adds it: `Tmdot = advection + diffusion + taugradu`. Both solvers call it.

The heating is the full stress power, `Phi = tau:D`. This is the entropy-elastic (rubber-like) convention: elastic work is released as heat under loading and absorbed under unloading.

With the local strain rate `a = R^2*Rdot/r^3 = (Rdot/R)*yT^-3` and `s = tau_rr - tau_thetatheta`, the stress power is `tau:D = -2*a*s`. Only this normal-stress difference enters; the trace of the stress does not.

## Heating per model

| `stress` | `taugradu / Br` |
|---|---|
| 0 | `0` |
| 1 | `12*mu_KV*a^2 - 2*a*se`, neo-Hookean `se` |
| 2 | the same with the quadratic `se` |
| 3, 4 | `12*mu_s*a^2 - 2*a*sp`, with the Zener polymer field `sp` of [`imr.physics.zener_field`](../reference/functions/f_zener_field.md) |
| 5, 6, 7 | `12*mu_s*a^2 - 2*a*(tau_rr - tau_tt)`, from the spectral stress field (UCM/Oldroyd-B, PTT or Giesekus) |
| 8 | `12*mu_KV*a^2 - 2*a*se_wall`, with `se_wall` a wall-level (not locally resolved) approximation; see [below](#kelvin-voigt-yang-church-stress-8) |

**Viscosities**, with `fnu` the non-Newtonian factor at the medium points:

- `mu_KV = 1/Re8 + fnu*iDRe`: stress 1, 2 and 8, where the whole viscosity acts instantly.
- `mu_s = LAM/Re8 + fnu*iDRe`: stress 3–7, where only the solvent part acts instantly.

`iDRe` is `1/DRe` for stress 1, 2 and 8 and `LAM/DRe` for stress 3–7, the same factor as in the radial equation (set in `imr.configuration.prepare_case`).

**Elastic stress difference.** With the hoop stretch `lam = r/r0`, where `r0^3 = r^3 - R^3 + Req^3`:

- neo-Hookean: `se = (lam^-4 - lam^2)/Ca`;
- quadratic: the same times `1 + alphax*(lam^-4 + 2*lam^2 - 3)`.

In the code, `lam^2 = yT2./x2` and `lam^-4 = x2.^2.*iyT4`, with `x2 = (yT3 - 1 + (Req/R)^3).^(2/3)`.

## Total heating

Integrating over the medium gives `4*pi*R^3*int_1^Inf y^2*Phi dy = -4*pi*R^2*Rdot*S`, with `S` the stress integral of the radial equation. This holds for every model. For stress 3 and 4 it equals `16*pi*LAM/Re8*R*Rdot^2 + 24*pi*R^2*Rdot*intfnu*iDRe - 4*pi*Rdot*Z1/R`.

## Equation

On the grid `y = r/R`, the field obeys

    d(sp)/dt|_y = -ydot*d(sp)/dy - (sp - sKV)/De - 3*(Rdot/R)*(1 - y^-3)*sp
    ydot = (Rdot/R)*(y^-2 - y)
    sKV  = se - 6*etap*a,   a = (Rdot/R)*y^-3,   etap = (1 - LAM)/Re8

- `se` is the neo-Hookean elastic stress difference for stress 3, or the quadratic one for stress 4 (see [`imr.physics.stress_power`](../reference/functions/f_stress_dissipation.md)).
- The caller computes the advection `-ydot*d(sp)/dy` on its own grid.
- **Boundary conditions:** none at the wall, where `ydot = 0`. `sp` is 0 at infinity.
- **Initial condition:** `sp = 0` from rest.

## Why this law

It is the local law whose integral reproduces the solver's stress state exactly: `2*int_1^Inf sp/y dy = Z1/R^3`, where `Z1` obeys `De*Z1dot + Z1 = Ze - 4*etap*R^2*Rdot` with `Ze = R^3*Se`.

**Proof.** Written for the weighted field `w = y^-3*sp`, the law reads `Dw/Dt = -(w - y^-3*sKV)/De` along a material point. From rest its solution is

    y^-3*sp(t) = int_0^t K(t-t')*y'^-3*sKV(t') dt',   K(t) = exp(-t/De)/De

where primes denote the same material point at time `t'`.

Incompressibility gives `r^3 - r'^3 = R^3 - R'^3`, so `dr/r = r'^2*dr'/r^3`. Hence

    2*int_R^Inf sp/r dr = (1/R^3)*int_0^t K*R'^3*[2*int_R'^Inf sKV/r' dr'] dt'
                        = (1/R^3)*int_0^t K*R'^3*S_KV(t') dt'

with `S_KV = Se - 4*etap*Rdot/R`, the stress integral of `sKV`. This is `Z1/R^3`, the solution of the `Z1` equation from rest.

The `Z1` equation fixes only this weighted integral. The law above is the simplest local field consistent with it, not the only one. It sets the radial distribution of the heating, but neither `R(t)` nor the total heating.

So:

- `Z1` still drives `R(t)`, unchanged.
- The field only distributes the heating in `r`.
- The total heating is fixed by `Z1`.

## In the solvers

- **`imr.solve_spectral`.** `sp` is expanded in the stress 5 basis `T_n(z) - 1` with `Nv` modes, stored after `Z1`, and evaluated at the medium points with `ZZT`. After a run, the solver prints the largest mismatch of `2*cdd'*e` against `Z1/R^3` as a check.
- **`imr.solve_nodal`.** `sp` is stored at the medium nodes `yT` after `Z1`: the first node is the wall, and the last, at infinity, is held at 0. The advection reuses the `Rdot` part of the medium temperature advection, `xi1_sq./(Lt*R).*(Rdot./yT2.*om_yT3/2).*(dmatrix_Tm*sp)`, which equals `-ydot*d(sp)/dy`.
- **Collapse (nodal only).** The history of the field before `t = 0` is not computed. As a stopgap, the initial field takes the shape of the elastic stress difference of the initial radius, scaled so that the integral identity holds for the initial `Z1` from `imr.physics.initial_stress`. Without elasticity (`Ca = Inf`), it takes the viscous shape `y^-3` instead.

## Spectral stress field

For stress 5, 6 and 7, `imr.solve_spectral` expands `tau_rr` and `tau_tt` in `T_n(z) - 1`, `n = 1..Nv`, with
`z = 1 - 2/(1 + (r/R - 1)/Lv)` and collocation points `z_j = cos(pi*j/Nv)`.
If `c` and `d` are their coefficients, the stress integral is
`S = 2*cdd'*(c - d) - 4*LAM/Re8*Rdot/R - 6*intfnu*iDRe`, the same for all three models: only the
local evolution of `c` and `d` differs (see [stress and viscosity](../model/stress-viscosity.md#ptt-and-giesekus-stress-6-7) for the PTT and Giesekus terms added to the UCM/Oldroyd-B closure, from [Warnez and Johnsen (2015)](../reference/references.md)).

The weights are
`cdd(n) = 2*Lv*int_{-1}^{1} (T_n(x) - 1)/((Lv*(1 + x) + 1 - x)*(1 - x)) dx`.
One Gauss–Legendre rule with `Nv + 64` points in `theta = acos(x)` evaluates them, where the integrands are smooth. The local field supplies the polymer part of medium stress power. See [reported comparisons](../validation/convergence.md#reported-spectral-stress-comparisons) for historical quadrature, timing and resolution measurements.

Stress 3–4 keep their integral `Z1` and use the same basis for the local field when medium heat is on. For stress 5 without that field, only `Z1` and `Z2` are available; the [compatibility rules](../reference/compatibility.md) explain why these cannot supply local heating. Stress 6 and 7 always run with the field (`spectral = 1` is forced), so this restriction does not apply to them.

## Kelvin-Voigt (Yang-Church, stress 8)

Stress 1 and 2's elastic heating term above is the fully local one: `se` evaluated at the stretch of every medium point `y`, the rigorous choice for a hyperelastic material. Stress 8's linear spring does not (yet) have that treatment.

A first attempt derived a candidate local field by matching stress 8's own radial-equation integral, `S = 2*int_1^Inf se(lam(y))/y dy`, the same relation that reproduces stress 1's `se = (lam^-4-lam^2)/Ca` exactly. Linearizing the incompressible kinematics for small strain gives `se = 6/Ca*(1-lam)`, but a numerical check of the integral against stress 8's own closed form (`S = -4/(3*Ca)*(1-Rst^3)`) shows the two disagree away from `lam = 1` — the leading-order linearization is not the exact local field, and no exact one has been found.

Instead, `imr.physics.stress_power`'s stress 8 branch uses the wall-level (not locally resolved beyond the shared `1/y^6` geometric decay every model's viscous term already has) approximation of `funIMRsolver.m` (`InertialMicrocavitationRheometry/IMR_v1`) and an older spectral reference in this project's own history, both of which compute this model's medium heating the same way without a local field:

```text
taugradu = 12*Br*mu_KV*(Rdot/R)^2/y^6 + 4*Br*(Rdot/R)*(1 - Rst^3)/(Ca*y^6)
```

The viscous term's `Br`/`Re8` normalization matches `imr.physics.stress_power`'s own exactly (confirmed by expanding both forms), which is why this elastic term is trusted enough to port despite not being derived from first principles here. Upgrading this to a fully local field, if one exists for this spring, is future work.
