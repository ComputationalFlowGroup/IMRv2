# `f_dual_ops.m` Documentation

## Overview
`d = f_dual_ops(n)` returns the operations used by the analytic Jacobians [`f_jacobian_nodal`](f_jacobian_nodal.md) and [`f_jacobian_spectral`](f_jacobian_spectral.md). A quantity is carried as a struct with its value `v` and its sparse derivative `J = dv/dX` with respect to the `n`-dimensional state `X`, so composing the operations applies the chain rule exactly.

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

Scalars and vectors may be mixed; a scalar is broadcast where needed.
