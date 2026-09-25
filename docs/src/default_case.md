# `default_case.m` Documentation

## Overview
`default_case` is a function that returns the default value of every solver input as a struct. `f_call_params` starts from these values, and any of them can be changed with a name/value pair, for example `f_imr_nodal('radial', 2, 'stress', 1)`. **Every quantity is dimensionless.** The solvers take no dimensional inputs and return no dimensional outputs; convert physical values beforehand, for example with [`toolchain/f_dim2nondim.m`](../toolchain/f_dim2nondim.md).

### Scalings
| Quantity | Scale |
|----------|-------|
| length | initial radius `R0` |
| pressure, stress | far-field pressure `P8` |
| velocity | `Uc = sqrt(P8/rho8)` |
| time | `t0 = R0/Uc` |
| temperature | far-field temperature `T8` |

The default values correspond to a 50 µm bubble in water at 1 atm and 298.15 K.

---

## Parameters

The *input* column is the name used in the name/value call (case-insensitive).

### Equation options
| Input | Variable | Description | Default |
|-------|----------|-------------|---------|
| `radial` | `radial` | 1 Rayleigh–Plesset, 2 Keller–Miksis (pressure), 3 Keller–Miksis (enthalpy, Tait), 4 Gilmore (Tait), 5 Keller–Miksis (Mie–Grüneisen), 6 Gilmore (Mie–Grüneisen). Radial 3–6 do not allow medium heat transfer, and radial 5 and 6 support only the polytropic bubble. | 1 |
| `bubtherm` | `bubtherm` | Thermal PDE in the bubble (0/1) | 0 |
| `medtherm` | `medtherm` | Thermal PDE in the medium (0/1, requires `bubtherm = 1` and `radial` 1 or 2) | 0 |
| `stress` | `stress` | 0 none, 1 NHKV, 2 quadratic NHKV, 3 Maxwell/Jeffreys/Zener (NH), 4 same (quadratic NH), 5 UCM/Oldroyd-B | 0 |
| `eps3` | `eps3` | Second-order strain term, in [0, 0.5]; currently unused by the stress models | 0 |
| `vapor` | `vapor` | Vapor pressure (0/1) | 0 |
| `masstrans` | `masstrans` | Vapor mass transfer (0/1, requires `bubtherm = vapor = 1`) | 0 |

### Solver options
| Input | Variable | Description | Default |
|-------|----------|-------------|---------|
| `tfin` | `TFin` | Final time, `TFin/t0` | 3.9034 |
| `tvector` | `TVector` | Output times, `t/t0` (give `tvector` or `tfin`, not both) | `linspace(0,TFin,128)` |
| `method` | `method` | Stiff time integrator: 15 `ode15s`, 23 `ode23tb` | 23 |
| `spectral` | `spectral` | Spectral stress representation for `stress = 5`, `f_imr_spectral` only (0/1) | 0 |
| `divisions` | `divisions` | Minimum number of time steps (0 automatic) | 0 |
| `nt`, `mt`, `nv` | `Nt`, `Mt`, `Nv` | Grid points: bubble, medium, stress. `Nt, Mt >= 3` when the thermal models are on; `Nv` is used only by the spectral stress representation | 25, 25, 150 |
| `lv`, `lt` | `Lv`, `Lt` | Grid stretching lengths | 3, 2 |
| `reltol`, `abstol` | `rtol`, `atol` | Integrator tolerances | 1e-8, 1e-8 |
| `jacobian` | `jacobian` | Jacobian for the stiff integrators, both solvers: 0 dense numerical, 1 sparsity pattern (nodal solver; numerical in the spectral solver), 2 analytic | 2 |
| `scheme` | `scheme` | Differentiation matrices of `f_imr_nodal`: 0 finite difference, 1 Chebyshev. `f_imr_spectral` takes no scheme | 0 |
| `progdisplay` | `progdisplay` | Print progress (every right-hand-side call; slows runs) | 0 |

### Initial conditions
| Input | Variable | Definition | Default |
|-------|----------|------------|---------|
| `collapse` | `collapse` | Out-of-equilibrium collapse initial condition (0/1). Requires `bubtherm = medtherm = vapor = masstrans = 1`, so it is available only in `f_imr_nodal`. | 0 |
| `req` | `Req` | `Req/R0` | 0.2 |
| `rdotzero` | `Rdotzero` | `Rdot(0)/Uc` | 0 |
| `p0star` | – | `P0/P8`; computed from equilibrium at `Req` if not given | computed |
| `stress0` | `Szero` | Initial stress state: empty for stress 0–2, one value for stress 3–4, `2*Nv` values for stress 5 with `spectral = 1` and 2 values with `spectral = 0`; computed if not given (from the collapse for stress 3–4 with `collapse = 1`, zero otherwise) | computed |

### Acoustic options
| Input | Variable | Definition | Default |
|-------|----------|------------|---------|
| `cstar` | `Cstar` | `C8/Uc` | 152.071 |
| `gama` | `GAMa` | Tait parameter `GAM/P8` | 3009.26 |
| `nstate` | `nstate` | Tait exponent | 7.15 |
| `hugoniot_s` | `hugoniot_s` | Hugoniot slope (Mie–Grüneisen) | 1.65 |

### Pressure waveform
| Input | Variable | Definition | Default |
|-------|----------|------------|---------|
| `wave_type` | `wave_type` | -1/-2/-3 measured waveforms (`d_hn`, `d_mn`, `d_ml`, already dimensionless, from [Hersey et al. 2023](https://doi.org/10.1121/10.0017342)), 0 impulse, 1 Gaussian, 2 histotripsy, 3 Heaviside | 0 |
| `ee` | `ee` | Amplitude `pA/P8` | 0 |
| `om` | `om` | Frequency `omega*t0` | 0 |
| `tw` | `tw` | Width `TW/t0` | 0 |
| `dt` | `dt` | Delay `DT/t0` | 0 |
| `mn` | `mn` | Histotripsy power | 0 |

### Stress and viscosity
| Input | Variable | Definition | Default |
|-------|----------|------------|---------|
| `we` | `We` | Weber `P8*R0/(2*S)` (Inf: no surface tension) | 36.1875 |
| `re` | `Re8` | Reynolds `P8*R0/(mu8*Uc)` (Inf: inviscid) | 623.366 |
| `ca` | `Ca` | Cauchy `P8/G` (Inf: no elasticity) | 101.325 |
| `de` | `De` | Deborah `lambda1*Uc/R0` | 0.019517 |
| `lam` | `LAM` | `lambda2/lambda1` | 0.1 |
| `alphax` | `alphax` | Strain-stiffening parameter | 0.25 |
| `dre` | `DRe` | `P8*R0/((mu0-mu8)*Uc)`, 0 for no non-Newtonian contribution | 0 |
| `v_a`, `v_nc` | `v_a`, `v_nc` | Exponents of the non-Newtonian models (see `f_viscosity`) | 0, 0 |
| `v_lambda` | `v_lambda_star` | Viscosity time scale `v_lambda/t0` | 0 |
| `nu_model` | `nu_model` | 0 Newtonian, 1 Carreau, 2 Carreau–Yasuda, 3 Powell–Eyring, 4 modified Powell–Eyring, 5 Cross, 6 simplified Cross, 7 modified Cross (see `f_viscosity`) | 0 |

### Thermal transport
| Input | Variable | Definition | Default |
|-------|----------|------------|---------|
| `kappa` | `kappa` | Ratio of specific heats | 1.47 |
| `foh` | `Foh` | Thermal Fourier `Dm/(Uc*R0)` | 2.5339e-4 |
| `br` | `Br` | Brinkman `Uc^2/(Cp*T8)` | 7.6394e-5 |
| `alpha_g`, `beta_g` | | Gas conductivity `K = K8*(alpha_g*T + beta_g)` | 0.57610, 0.42634 |
| `alpha_v`, `beta_v` | | Vapor conductivity coefficients | 0.36006, 0.63750 |
| `chi` | `chi` | `T8*K8/(P8*R0*Uc)` | 0.16479 |
| `kmstar` | `Km_star` | Medium to gas conductivity ratio `Km/K8` | 20.1276 |

### Mass transfer and vapor pressure
| Input | Variable | Definition | Default |
|-------|----------|------------|---------|
| `fom` | `Fom` | Mass Fourier `D0/(Uc*R0)` | 0.049597 |
| `lheat` | `L_heat_star` | Latent heat `L_heat/Uc^2` | 23781.9 |
| `rv`, `ra` | `Rv_star`, `Ra_star` | Gas constants `R/(P8/(rho8*T8))` | 1444.95, 898.681 |
| `pv` | `Pv_sat` | Saturated vapor pressure at `T8`, `Pv/P8` | 0.030760 |
| `thetav` | `thetav` | Vapor pressure exponent: `Pv(T)/P8 = pv*exp(thetav*(1 - 1/T))` | 17.4409 |

---

## Notes
- All inputs are validated by [`f_checkers`](f_checkers.md); inputs that look dimensional raise the `IMR:dimensional` warning.
- There are no case files: a case is the list of inputs that differ from these defaults, kept for example as a cell array `varin` and passed as `f_imr_nodal(varin{:})`. The former `casefile` input now raises an error that says so. To run many cases, see [`f_imr_sweep`](f_imr_sweep.md).
- See [`compatibility`](../compatibility.md) for which features run together in each solver.
