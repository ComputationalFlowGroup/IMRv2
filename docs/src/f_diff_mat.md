# `f_diff_mat.m` Documentation

## Overview
`[y, D] = f_diff_mat(nodes, order, tm_check, scheme)` returns the node positions and the differentiation matrix used by the nodal thermal solver, for either spatial discretization. Swapping `scheme` is the only difference between the finite difference and Chebyshev versions of the solver.

| Argument | Meaning |
|---|---|
| `nodes` | number of grid points |
| `order` | 1 for `d/dy`, 2 for the second-derivative operator |
| `tm_check` | 0: interior (bubble) grid on `y = r/R` in [0, 1], where the order-2 operator is the spherical Laplacian `d2/dy2 + (2/y) d/dy` with the symmetry condition at `y = 0`; 1: exterior (medium) grid on `xi` in [1, -1] |
| `scheme` | 0 finite differences (uniform grid), 1 Chebyshev collocation |

## Finite differences (`scheme = 0`)
Identical to `f_finite_diff_mat` in `f_imr_fd`: second-order central differences in the interior, `(1 ± dy/y)/dy^2` for the Laplacian, `(-6, 6)/dy^2` at the center, the one-sided `(0.5, -2, 1.5)/dy` wall row for the interior first derivative, and zero boundary rows otherwise. Node positions are `y = (0:N)/N` (interior) and `xi = 1 - 2*(0:N)/N` (exterior).

## Chebyshev (`scheme = 1`)
- Interior: the operators are built on the doubled Chebyshev grid of `2*nodes - 1` points on [-1, 1] and folded onto the half `y >= 0` for even functions, ordered from the center `y = 0` to the wall `y = 1`. The first-derivative row at the center is zero exactly (`f'(0) = 0`). The order-2 operator is the spherical Laplacian, with `(2/y) d/dy` added at every node except the center and `3 d2/dy2` at the center, built on the doubled grid before folding.
- Exterior: the standard Chebyshev differentiation matrix on the points `cos(pi*j/N)` of [-1, 1] (ordered from 1 to -1) and its square.

## Verification (`toolchain/s_test_suite_diffmat.m`)
`scheme = 0` reproduces the solver's operator bitwise for `nodes` in {10, 100, 200}, both orders and both grids. The Chebyshev interior operators differentiate `cos(pi*y/2)` and `y^2` to 1e-11 or better for 8–32 nodes and annihilate constants to 1e-13; the exterior operators converge spectrally (`exp(xi)`: 6e-6 at 8 nodes, 4e-14 at 16 nodes).
