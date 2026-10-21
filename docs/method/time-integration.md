# Time integration and Jacobians

`f_odesolve` integrates the coupled state with a stiff solver. `jacobian` chooses derivative evaluation, while `scheme` chooses nodal spatial operators.

## Time stepping

The supported integrators, codes and MATLAB requirements live in the [parameter table](../model/parameters.md#solver-options). Thermal, vapor and stress fields can make the state stiff, so the wrapper does not offer a nonstiff integrator.

| Setting | Implementation |
|---|---|
| `reltol`, `abstol` | `RelTol` and `AbsTol` in `f_odesolve` |
| `divisions = 0` | Automatic time-step selection |
| Positive `divisions` | `MaxStep = tspan(end)/divisions` |
| Two-element time span | Return the integrator's output steps |
| Longer time vector | Return values at requested times |

CVODES uses MATLAB's `ode` object with `Solver = 'cvodesstiff'`. `f_odesolve` assigns object properties separately so older MATLAB releases can parse the file even when that object is absent. It constructs `odeJacobian` from a derivative function or sparsity pattern. CVODES compiles the integrator's step control and Newton iterations, while the right-hand side and Jacobian still execute MATLAB code. The [benchmark](integrator-benchmarks.md) measures the combination.

## Derivative choices

| `jacobian` | Nodal solver | Spectral solver |
|---|---|---|
| 0 | Dense numerical | Dense numerical |
| 1 | Sparsity pattern, dense for Chebyshev | Dense numerical |
| 2 | Hybrid analytic with bubble heat; sparsity fallback otherwise | Hybrid analytic |
| 3 | Batched finite differences | Rejected |

The hybrid Jacobians use `f_dual_ops` to carry each field quantity as its value and sparse state derivative. The chain rule through transport and thermal transforms is exact in those operations. Radial, constitutive and other low-dimensional closure partials use central differences; “analytic” does not mean every operation is differentiated symbolically.

`f_jacobian_nodal` differentiates wall temperature implicitly through `f_wall_mass_bc`, or through the exact quadratic thermal root. Its derivative includes every node used by the wall row. Zener-field rows are linear in their own field values; those entries are exact, while their radius and velocity partials use central differences. The heating's field derivative is `-2*Br*a` on the diagonal.

`f_jacobian_spectral` differentiates the boundary solve
`A(s)*q = rhs`, with `A(s) = A0 - (s-1)*e1*v'`, as
`dq = A(s)\(d(rhs) + e1*(v'*q)*ds)`.
It reuses the prefactored `A0` and Sherman–Morrison correction. The local polymer heating contributes `-2*Br*a*ZZT` in field coefficients. Radius and velocity partials of heating and Zener transport use central differences at fixed field.

`f_jpattern` builds a conservative coupling superset. With finite differences, interior operators are tridiagonal but pressure and wall stencils couple fields to radius, velocity and pressure. Chebyshev operators couple all nodes of each field, so their pattern is dense.

## Batched right-hand side

The nodal `SVBDODE` accepts a matrix with one state per column. Its helpers evaluate each column elementwise in the same arithmetic order. The time integrators call it with one column; `jacobian = 3` uses a batch.

It forms `F = SVBDODE(t, [X, X + diag(del)])` and
`J = (F(:,2:end) - F(:,1))./del.'`.
The perturbation is `sqrt(eps)*max(abs(X),1)`, rounded to an exactly representable state increment. It zeros entries outside `f_jpattern` and returns a dense Jacobian for dense LU.

Batching amortizes interpreted operations across columns. The earlier implementation description counted 55 sparse dual-number operations for the hybrid analytic path; actual timing depends on configuration and MATLAB. [Reported measurements](../validation/convergence.md#reported-batched-jacobian-measurements) preserve historical comparisons without claiming a new speedup.

## Thermal variables

Both solvers carry a Kirchhoff variable for linear conductivity:

$$\theta=\alpha(T^2-1)/2+\beta(T-1),\qquad
T=\frac{\sqrt{(\alpha+\beta)^2+2\alpha\theta}-\beta}{\alpha}.$$

The spectral solver uses gas coefficients. The nodal solver uses the initial mixture without mass transfer and local vapor-dependent mixture with mass transfer. Both use the supplied `alpha` and `beta`; their sum need not equal 1. The wall-condition derivative must use the same transform. See [wall boundaries](wall-boundaries.md) and the [reported correction](../validation/convergence.md#reported-kirchhoff-transform-correction).
