# `imr.numerics.dual_operations`

`d = imr.numerics.dual_operations(n)` returns the operations used by the analytic Jacobians [`imr.numerics.jacobian_nodal`](jacobian_nodal.md) and [`imr.numerics.jacobian_spectral`](jacobian_spectral.md). Each quantity is a struct with its value `v` and its sparse derivative `J = dv/dX` with respect to the `n`-dimensional state `X`. Composing the operations therefore applies the chain rule exactly.

## Operations

| Field | Result |
|---|---|
| `var(X, idx)` | the state entries `X(idx)`, with identity derivative rows |
| `con(v)` | a constant, with zero derivative |
| `el(b, k)`, `set(a, k, b)` | read or write entries `k` |
| `lin(M, b)` | `M*b` for a matrix `M` |
| `fun(b, val, dval)` | an elementwise function with value `val` and derivative `dval` |
| `scale(b, s)` | `s.*b` for a constant `s` |
| `add`, `sub`, `mul`, `div` | elementwise arithmetic of two quantities |
| `diag(v)` | a sparse diagonal matrix |

You can mix scalars and vectors. The operations broadcast a scalar where needed.
