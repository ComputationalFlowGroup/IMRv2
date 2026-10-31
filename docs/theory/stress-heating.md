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
| 8 | the same with the 2nd-order Fung `se` |
| 10 | the same with the Mooney-Rivlin `se` |
| 3, 4, 9, 11 | `12*mu_s*a^2 - 2*a*sp`, with the Zener polymer field `sp` of [`imr.physics.zener_field`](../reference/functions/f_zener_field.md) |
| 5, 6, 7 | `12*mu_s*a^2 - 2*a*(tau_rr - tau_tt)`, from the spectral stress field (UCM/Oldroyd-B, PTT or Giesekus) |

**Viscosities**, with `fnu` the non-Newtonian factor at the medium points:

- `mu_KV = 1/Re8 + fnu*iDRe`: stress 1, 2, 8 and 10, where the whole viscosity acts instantly.
- `mu_s = LAM/Re8 + fnu*iDRe`: stress 3–7, 9 and 11, where only the solvent part acts instantly.

`iDRe` is `1/DRe` for stress 1, 2, 8, 10 and `LAM/DRe` for stress 3–7, 9, 11, the same factor as in the radial equation (set in `imr.configuration.prepare_case`).

**Elastic stress difference.** With the hoop stretch `lam = r/r0`, where `r0^3 = r^3 - R^3 + Req^3`, and `x = lam^-4 + 2*lam^2 - 3` (the invariant term `I1 - 3`):

- neo-Hookean (1, 3): `se = (lam^-4 - lam^2)/Ca`;
- quadratic (2, 4): the same times `1 + alphax*x`, the first-order term of the series `exp(alphax*x)`;
- 2nd-order Fung (8, 9): the same times `1 + alphax*x + (alphax^2/2)*x^2`, the next term of the same series -- unlike stress 8's own wall (radial-equation) formula, which was cross-checked against an independent codebase, this local field is derived here (the natural next term of the series) and verified only by its own consistency: it reproduces the independently-sourced wall formula exactly under `S = 2*int_1^Inf se/y dy`, and its `Rst`-derivative matches the wall formula's own `Sdot` (see `imr.physics.material_stress`);
- Mooney-Rivlin (10, 11): `se = (1/Ca)*[(lam^-4 - lam^2) - C2*(lam^4 - lam^-2)]`, the neo-Hookean term plus one new term scaled by the second-modulus ratio `C2`; this field and its wall integral are both derived here from the two-invariant Mooney-Rivlin strain energy, with no external source, and cross-checked the same way (integral and derivative consistency).

In the code, `lam^2 = yT2./x2` and `lam^-4 = x2.^2.*iyT4`, with `x2 = (yT3 - 1 + (Req/R)^3).^(2/3)`; `lam^4 = (yT2.*ix2).^2` and `lam^-2 = x2./yT2` for the Mooney-Rivlin term.

## Total heating

Integrating over the medium gives `4*pi*R^3*int_1^Inf y^2*Phi dy = -4*pi*R^2*Rdot*S`, with `S` the stress integral of the radial equation. This holds for every model. For stress 3, 4, 9 and 11 it equals `16*pi*LAM/Re8*R*Rdot^2 + 24*pi*R^2*Rdot*intfnu*iDRe - 4*pi*Rdot*Z1/R`.

## Equation

On the grid `y = r/R`, the field obeys

    d(sp)/dt|_y = -ydot*d(sp)/dy - (sp - sKV)/De - 3*(Rdot/R)*(1 - y^-3)*sp
    ydot = (Rdot/R)*(y^-2 - y)
    sKV  = se - 6*etap*a,   a = (Rdot/R)*y^-3,   etap = (1 - LAM)/Re8

- `se` is the neo-Hookean elastic stress difference for stress 3, the quadratic one for stress 4, the 2nd-order Fung one for stress 9, or the Mooney-Rivlin one for stress 11 (see [`imr.physics.stress_power`](../reference/functions/f_stress_dissipation.md)).
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

Stress 3, 4, 9 and 11 keep their integral `Z1` and use the same basis for the local field when medium heat is on. For stress 5 without that field, only `Z1` and `Z2` are available; the [compatibility rules](../reference/compatibility.md) explain why these cannot supply local heating. Stress 6 and 7 always run with the field (`spectral = 1` is forced), so this restriction does not apply to them.
