# `imr.numerics.jacobian_pattern`

`Sp = imr.numerics.jacobian_pattern(n, ith, itm, ikv, iz, scheme)` builds the Jacobian sparsity pattern of the nodal solver from its state indices. `imr.solve_nodal` passes it to the stiff integrators when `jacobian = 1`, and when `jacobian = 2` without bubble heat transfer, where no analytic Jacobian exists.

## Inputs

- `n`: length of the state vector `[R, Rdot, P, theta, Tm, kv, Z]`.
- `ith`, `itm`, `ikv`, `iz`: indices of the bubble temperature, medium temperature, vapor concentration and stress states.
- `scheme`: the differentiation scheme. For Chebyshev (`scheme = 1`) every node of a field couples to every other, so the pattern is dense.

## Notes

- The pattern is a conservative superset of the true Jacobian. The bubble pressure equation and the wall boundary condition feed every field row, so all field rows couple to `R`, `Rdot`, `P` and the wall stencil columns, while the finite difference operators themselves are tridiagonal.
