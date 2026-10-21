# Spatial discretization

`f_imr_nodal` carries field values at nodes; `f_imr_spectral` carries Chebyshev expansion coefficients. Select resolution through [parameters](../model/parameters.md) and check convergence for the chosen case.

## Finite differences (`scheme = 0`)

`f_diff_mat` builds second-order central differences. With `N = nodes - 1`, the interior nodes are `y = (0:N)/N` and exterior nodes are `xi = 1 - 2*(0:N)/N`. The medium mapping uses the stretching length `Lt`.

| Operator row | Stencil or treatment |
|---|---|
| Interior first derivative | Centered difference away from boundaries |
| Spherical Laplacian | Neighbor coefficients `(1 ± dy/y)/dy^2` |
| Bubble center Laplacian | `(-6, 6)/dy^2`, enforcing spherical symmetry |
| Bubble-wall first derivative | One-sided `(0.5, -2, 1.5)/dy` |
| Other boundary rows | Zero where the solver imposes boundary values separately |

## Chebyshev (`scheme = 1`)

The bubble operators start with a doubled Chebyshev grid of `2*nodes - 1` points on `[-1,1]`. Folding even functions onto `y >= 0` orders them from center to wall. The first-derivative center row is exactly zero. The Laplacian is `d2/dy2 + (2/y)d/dy` away from the center and `3*d2/dy2` at the center, formed before folding.

The exterior operators are the standard Chebyshev matrix and its square on `cos(pi*j/N)`, ordered from 1 to -1. Global operators couple all nodes. Wall heat and vapor balances therefore use every field node; finite-difference wall rows use three points. The physics and time integration are shared between nodal schemes.

## Spectral expansions

`f_imr_spectral` uses an even Chebyshev expansion for bubble temperature and a mapped Chebyshev expansion for the medium. Operators `gDA`, `gDDA`, `mDA` and `mDDA` map coefficients to derivative values and are formed once.

The thermal boundary system has a single time-dependent row. The solver factors the constant matrix once and applies a Sherman–Morrison rank-one correction. Stress fields use the mapped basis described in [stress heating theory](../theory/stress-heating.md#spectral-stress-field).

Both solvers precompute constant medium-equation and thermal-map factors. For `wave_type = 0` they use the constant forcing directly instead of calling `f_pinfinity` at every right-hand-side evaluation.

The [operator suite](../validation/operators.md) checks analytic derivatives and stencil definitions. [Reported convergence](../validation/convergence.md) distinguishes historical observed errors from those executable thresholds. Strong collapse, sharp forcing and high relaxation time scales can require finer grids.
