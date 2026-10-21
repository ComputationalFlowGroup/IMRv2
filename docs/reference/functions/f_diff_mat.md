# `f_diff_mat`

`[y, D] = f_diff_mat(nodes, order, tm_check, scheme)` returns the node positions and the differentiation matrix for either spatial discretization. It is the operator module of the nodal solver [`f_imr_nodal`](../../api/f_imr_nodal.md), whose `scheme` input is passed straight through: swapping it is the only difference between the finite difference and Chebyshev versions of the solver.

| Argument | Meaning |
|---|---|
| `nodes` | number of grid points |
| `order` | 1 for `d/dy`, 2 for the second-derivative operator |
| `tm_check` | 0: interior (bubble) grid on `y = r/R` in [0, 1], where the order-2 operator is the spherical Laplacian `d2/dy2 + (2/y) d/dy` with the symmetry condition at `y = 0`; 1: exterior (medium) grid on `xi` in [1, -1] |
| `scheme` | 0 finite differences (uniform grid), 1 Chebyshev collocation |

For the algorithm and assumptions, see the [topic guide](../../method/discretization.md). Numerical measurements are retained as [historical reports](../../validation/convergence.md).
