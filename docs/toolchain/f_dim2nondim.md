# `f_dim2nondim.m` Documentation

## Overview
`f_dim2nondim` is a side function that converts dimensional (SI) physical properties into the dimensionless inputs taken by `f_imr_nodal` and `f_imr_spectral`. The solvers never call it: they accept only dimensionless inputs.

```matlab
addpath('src'); addpath('toolchain');
[varin, ref] = f_dim2nondim('r0',50e-6,'req',5e-6,'mu',1e-2,'g',1e3, ...
    'tvector',linspace(0,12e-6,100));
[t, R] = f_imr_nodal(varin{:}, 'radial', 2, 'stress', 1);
plot(t*ref.t0, R*ref.R0);   % restore dimensions yourself
```

## Inputs (SI)
`r0`, `req`, `u0`, `tfin`, `tvector`, `rho8`, `p8`, `c8`, `gam`, `surft`, `g`, `mu`, `mu0`, `du`, `lambda1`, `lambda2`, `v_lambda`, `t8`, `atg`, `btg`, `atv`, `btv`, `km`, `cp`, `dm`, `dmass`, `lheat`, `rv`, `ra`, `pv`, `p0`, `pa`, `omega`, `tw`, `dt`. Defaults match water at 1 atm and 298.15 K with a 50 µm bubble. Any other name/value pair (model flags, grid sizes, `kappa`, `alphax`, `nstate`, `hugoniot_s`, ...) is passed through unchanged.

## Outputs
- `varin`: dimensionless name/value pairs (see [`default_case`](../src/default_case.md) for their definitions).
- `ref`: reference scales `R0`, `Uc`, `t0`, `P8`, `T8`, `rho8`.

The saturated vapor pressure comes from [`f_pvsat`](f_pvsat.md), which also supplies the vapor-pressure exponent `thetav = 5200/T8`.
