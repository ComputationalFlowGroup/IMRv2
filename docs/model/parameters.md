# Parameters and units

`default_case` is a function that returns the default value of every solver input as a struct. `f_prepare_case` starts from these values, and any of them can be changed with a name/value pair, for example `f_imr_nodal('radial', 2, 'stress', 1)`. **Every quantity is dimensionless.** The solvers take no dimensional inputs and return no dimensional outputs; convert physical values beforehand, for example with [`f_dim2nondim.m`](../api/f_dim2nondim.md).

### Scalings

| Quantity | Scale |
|----------|-------|
| length | initial radius `R0` |
| pressure, stress | far-field pressure `P8` |
| velocity | `Uc = sqrt(P8/rho8)` |
| time | `t0 = R0/Uc` |
| temperature | far-field temperature `T8` |

The default values use the legacy water reference set for a 50 µm bubble at 1 atm and 298.15 K.


## Parameters

The *input* column is the accepted name in a name/value call (case-insensitive). All public solver inputs are dimensionless, including time and temperature. Values in the tables are rounded for readability; Model and numerical defaults come from `src/default_case.m`; physical groups are derived from the shared SI converter and water profile. See [defaults and units](../reference/defaults.md). SI names, defaults and units are listed in the [conversion API](../api/f_dim2nondim.md#inputs-si).

### Equation options

| Input | Variable | Description | Default | Units |
|-------|----------|-------------|------------|---|
| `radial` | `radial` | 1 Rayleigh–Plesset, 2 Keller–Miksis (pressure), 3 Keller–Miksis (enthalpy, Tait), 4 Gilmore (Tait), 5 Keller–Miksis (Mie–Grüneisen), 6 Gilmore (Mie–Grüneisen). Radial 3–6 do not allow medium heat transfer, and radial 5 and 6 support only the polytropic bubble. | 1 | 1 (dimensionless) |
| `bubtherm` | `bubtherm` | Thermal PDE in the bubble (0/1) | 0 | 1 (dimensionless) |
| `medtherm` | `medtherm` | Thermal PDE in the medium (0/1, requires `bubtherm = 1` and `radial` 1 or 2) | 0 | 1 (dimensionless) |
| `stress` | `stress` | 0 none, 1 NHKV, 2 quadratic NHKV, 3 Maxwell/Jeffreys/Zener (NH), 4 same (quadratic NH), 5 UCM/Oldroyd-B | 0 | 1 (dimensionless) |
| `eps3` | `eps3` | Second-order strain term, in [0, 0.5]; currently unused by the stress models | 0 | 1 (dimensionless) |
| `vapor` | `vapor` | Vapor pressure (0/1) | 0 | 1 (dimensionless) |
| `masstrans` | `masstrans` | Vapor mass transfer (0/1, requires `bubtherm = vapor = 1`) | 0 | 1 (dimensionless) |

### Solver options

| Input | Variable | Description | Default | Units |
|-------|----------|-------------|------------|---|
| `tfin` | `TFin` | Final time, `TFin/t0` | 3.9034 | 1 (dimensionless) |
| `tvector` | `TVector` | Output times, `t/t0` (give `tvector` or `tfin`, not both) | `linspace(0,TFin,128)` | 1 (dimensionless) |
| `method` | `method` | Stiff time integrator: 23 `ode23tb`, 15 `ode15s`, 3 CVODES (MATLAB R2024a or later). The names `'ode23tb'`, `'ode15s'` and `'cvodes'` are also accepted | 23 | 1 (dimensionless) |
| `spectral` | `spectral` | Spectral stress representation for `stress = 5`, `f_imr_spectral` only (0/1) | 0 | 1 (dimensionless) |
| `divisions` | `divisions` | Minimum number of time steps (0 automatic) | 0 | 1 (dimensionless) |
| `nt`, `mt`, `nv` | `Nt`, `Mt`, `Nv` | Grid points: bubble, medium, stress. `Nt, Mt >= 3` when the thermal models are on; `Nv` is the number of stress modes of `f_imr_spectral`, for the stress 5 field (`spectral = 1`) and for the stress 3–4 polymer stress field with medium heat transfer | 25, 25, 150 | 1 (dimensionless) |
| `lv`, `lt` | `Lv`, `Lt` | Grid stretching lengths | 3, 2 | 1 (dimensionless) |
| `reltol`, `abstol` | `rtol`, `atol` | Integrator tolerances | 1e-8, 1e-8 | 1 (dimensionless) |
| `jacobian` | `jacobian` | Jacobian for the stiff integrators, both solvers: 0 dense numerical, 1 sparsity pattern (nodal solver; numerical in the spectral solver), 2 analytic, 3 batched finite differences (nodal solver only) | 2 | 1 (dimensionless) |
| `scheme` | `scheme` | Differentiation matrices of `f_imr_nodal`: 0 finite difference, 1 Chebyshev. `f_imr_spectral` takes no scheme | 0 | 1 (dimensionless) |
| `progdisplay` | `progdisplay` | Print progress (every right-hand-side call; slows runs) | 0 | 1 (dimensionless) |

### Initial conditions

| Input | Variable | Definition | Default | Units |
|-------|----------|------------|------------|---|
| `collapse` | `collapse` | Out-of-equilibrium collapse initial condition (0/1). Requires `bubtherm = medtherm = vapor = masstrans = 1`, so it is available only in `f_imr_nodal`. With stress 5 it also requires `stress0`. | 0 | 1 (dimensionless) |
| `req` | `Req` | `Req/R0` | 0.2 | 1 (dimensionless) |
| `rdotzero` | `Rdotzero` | `Rdot(0)/Uc` | 0 | 1 (dimensionless) |
| `p0star` | – | `P0/P8`; computed from equilibrium at `Req` if not given | computed | 1 (dimensionless) |
| `stress0` | `Szero` | Initial stress state: empty for stress 0–2, one value for stress 3–4, `2*Nv` values for stress 5 with `spectral = 1` and 2 values with `spectral = 0`; computed if not given (from the collapse for stress 3–4 with `collapse = 1`, zero otherwise; stress 5 with `collapse = 1` must give it) | computed | 1 (dimensionless) |

### Acoustic options

| Input | Variable | Definition | Default | Units |
|-------|----------|------------|------------|---|
| `cstar` | `Cstar` | `C8/Uc` | 152.071 | 1 (dimensionless) |
| `gama` | `GAMa` | Tait parameter `GAM/P8` | 3009.26 | 1 (dimensionless) |
| `nstate` | `nstate` | Tait exponent | 7.15 | 1 (dimensionless) |
| `hugoniot_s` | `hugoniot_s` | Hugoniot slope (Mie–Grüneisen) | 1.65 | 1 (dimensionless) |

### Pressure waveform

| Input | Variable | Definition | Default | Units |
|-------|----------|------------|------------|---|
| `wave_type` | `wave_type` | -1/-2/-3 measured waveforms (`d_hn`, `d_mn`, `d_ml`, already dimensionless, from [Hersey et al. 2023](https://doi.org/10.1121/10.0017342)), 0 impulse, 1 Gaussian, 2 histotripsy, 3 Heaviside | 0 | 1 (dimensionless) |
| `ee` | `ee` | Amplitude `pA/P8` | 0 | 1 (dimensionless) |
| `om` | `om` | Frequency `omega*t0` | 0 | 1 (dimensionless) |
| `tw` | `tw` | Width `TW/t0` | 0 | 1 (dimensionless) |
| `dt` | `dt` | Delay `DT/t0` | 0 | 1 (dimensionless) |
| `mn` | `mn` | Histotripsy power | 0 | 1 (dimensionless) |

### Stress and viscosity

| Input | Variable | Definition | Default | Units |
|-------|----------|------------|------------|---|
| `we` | `We` | Weber `P8*R0/(2*S)` (Inf: no surface tension) | 36.1875 | 1 (dimensionless) |
| `re` | `Re8` | Reynolds `P8*R0/(mu8*Uc)` (Inf: inviscid) | 623.366 | 1 (dimensionless) |
| `ca` | `Ca` | Cauchy `P8/G` (Inf: no elasticity) | 101.325 | 1 (dimensionless) |
| `de` | `De` | Deborah `lambda1*Uc/R0`; at least `1e-4` for stress 3–5 (below that, run the `De -> 0` limit: stress 1, 2, or 1 with `Ca = Inf`) | 0.019517 | 1 (dimensionless) |
| `lam` | `LAM` | `lambda2/lambda1` | 0.1 | 1 (dimensionless) |
| `alphax` | `alphax` | Strain-stiffening parameter | 0.25 | 1 (dimensionless) |
| `dre` | `DRe` | `P8*R0/((mu0-mu8)*Uc)`, 0 for no non-Newtonian contribution | 0 | 1 (dimensionless) |
| `v_a`, `v_nc` | `v_a`, `v_nc` | Exponents of the non-Newtonian models (see `f_viscosity`) | 0, 0 | 1 (dimensionless) |
| `v_lambda` | `v_lambda_star` | Viscosity time scale `v_lambda/t0` | 0 | 1 (dimensionless) |
| `nu_model` | `nu_model` | 0 Newtonian, 1 Carreau, 2 Carreau–Yasuda, 3 Powell–Eyring, 4 modified Powell–Eyring, 5 Cross, 6 simplified Cross, 7 modified Cross (see `f_viscosity`) | 0 | 1 (dimensionless) |

### Thermal transport

| Input | Variable | Definition | Default | Units |
|-------|----------|------------|------------|---|
| `kappa` | `kappa` | Ratio of specific heats | 1.47 | 1 (dimensionless) |
| `foh` | `Foh` | Thermal Fourier `Dm/(Uc*R0)` | 2.5339e-4 | 1 (dimensionless) |
| `br` | `Br` | Brinkman `Uc^2/(Cp*T8)` | 7.6394e-5 | 1 (dimensionless) |
| `alpha_g`, `beta_g` | | Gas conductivity `K = K8*(alpha_g*T + beta_g)`. `K8` is the mean of the gas and vapor conductivities at `T8`, so `alpha_g + beta_g` need not be 1; the solvers use both as given | 0.57610, 0.42634 | 1 (dimensionless) |
| `alpha_v`, `beta_v` | | Vapor conductivity coefficients | 0.36006, 0.63750 | 1 (dimensionless) |
| `chi` | `chi` | `T8*K8/(P8*R0*Uc)` | 0.16479 | 1 (dimensionless) |
| `kmstar` | `Km_star` | Medium to gas conductivity ratio `Km/K8` | 20.1276 | 1 (dimensionless) |

### Mass transfer and vapor pressure

| Input | Variable | Definition | Default | Units |
|-------|----------|------------|------------|---|
| `fom` | `Fom` | Mass Fourier `D0/(Uc*R0)` | 0.049597 | 1 (dimensionless) |
| `lheat` | `L_heat_star` | Latent heat `L_heat/Uc^2` | 23781.9 | 1 (dimensionless) |
| `rv`, `ra` | `Rv_star`, `Ra_star` | Gas constants `R/(P8/(rho8*T8))` | 1444.95, 898.681 | 1 (dimensionless) |
| `pv` | `Pv_sat` | Saturated vapor pressure at `T8`, `Pv/P8` | 0.030760 | 1 (dimensionless) |
| `thetav` | `thetav` | Vapor pressure exponent: `Pv(T)/P8 = pv*exp(thetav*(1 - 1/T))` | 17.4409 | 1 (dimensionless) |


## Input validation

- All inputs are validated by [`f_checkers`](../reference/functions/f_checkers.md); inputs that look dimensional raise the `IMR:dimensional` warning.
- [Configuration](../running/configuration.md) explains structs, cell arrays and precedence. There are no case files: a case is the list of inputs that differ from these defaults, kept as a struct from [`f_imr_inputs`](../api/f_imr_inputs.md) or as a cell array `varin` passed as `f_imr_nodal(varin{:})`. The former `casefile` input now raises an error that says so. To run many cases, see [`f_imr_sweep`](../api/f_imr_sweep.md).
- See [`compatibility`](../reference/compatibility.md) for which features run together in each solver.
