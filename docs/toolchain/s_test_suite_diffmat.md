# `s_test_suite_diffmat.m` Documentation

## Overview
`s_test_suite_diffmat.m` checks the differentiation operators returned by [`f_diff_mat`](../src/f_diff_mat.md). The CI workflow `test_operators.yml` runs it.

## Checks
- **Finite differences (`scheme = 0`):** the interior and exterior stencils against their defining values, including the symmetry row at the center, the one-sided wall row and the zero boundary rows. Tolerance `1e-12`.
- **Chebyshev interior (`scheme = 1`, 12 and 24 nodes):** the first derivative and the spherical Laplacian of `cos(pi*y/2)` and `y^2`, and that both operators annihilate constants. Tolerance `1e-8`.
- **Chebyshev exterior (12 and 24 nodes):** the first and second derivatives of `exp(xi)`, and that the first derivative annihilates constants. Tolerance `1e-8`.

## Outputs
- Console output with the maximum error of each check.
- Exit code `0` if all checks pass, `1` otherwise.

Adds `../src/` to the path, so it runs from `toolchain/`.
