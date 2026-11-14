# Differentiation operators

`s_test_suite_diffmat.m` checks the differentiation operators returned by [`imr.numerics.differentiation`](../reference/functions/differentiation.md). The `tests.yml` workflow runs this suite.

Nobody ran MATLAB for this documentation change. The author checked the descriptions below against the suite implementation. You must verify the numerical outcomes on your own MATLAB installation.


## Checks

- **Finite differences (`scheme = 0`):** the interior and exterior stencils against their defining values, including the symmetry row at the center, the one-sided wall row and the zero boundary rows. Tolerance `1e-12`.
- **Chebyshev interior (`scheme = 1`, 12 and 24 nodes):** the first derivative and the spherical Laplacian of `cos(pi*y/2)` and `y^2`, and that both operators annihilate constants. Tolerance `1e-8`.
- **Chebyshev exterior (12 and 24 nodes):** the first and second derivatives of `exp(xi)`, and that the first derivative annihilates constants. Tolerance `1e-8`.

## Outputs

- Console output with the maximum error of each check.
- An error if any check fails, which exits a batch run (`matlab -batch`, as in CI) with a nonzero code.

Run this suite from `tests/`, after you call `imr_test_setup`. `imr_test_setup` exposes the source and suites, without recursive path additions.

These assertions use analytic functions and stencil definitions. These assertions do not use stored simulations. The [reported comparisons](convergence.md#reported-measurements-imrnumericsdifferentiation) page keeps historical measurements with other node counts and tighter observed errors.
