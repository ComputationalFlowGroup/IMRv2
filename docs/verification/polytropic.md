# Constant-vapor polytropic verification

`tests/suites/s_test_suite_polytropic.m` checks the common gas/vapor pressure law in both solvers. From `toolchain/`, run this command: `matlab -batch "imr_test_setup; s_test_suite_polytropic"`.

For fixed ambient vapor pressure, total bubble pressure is

$$
P=P_{g0}R^{-3\kappa}+P_v,\qquad
\dot P=-3\kappa P_{g0}R^{-3\kappa}\frac{\dot R}{R}.
$$

$P_{g0}$ is the initial noncondensable gas pressure. $P_v$ is the constant vapor pressure. $\kappa$ belongs to the gas. The derivative differentiates only the gas pressure. The old expression wrongly used total pressure in the derivative. So the old expression differentiated a constant vapor contribution, which it should not have done.

The suite compares the derivative with finite differences. The suite verifies invariance to a constant vapor offset. The suite compares integrated pressure against the algebraic pressure law during motion, in nodal FD, nodal Chebyshev and spectral modes. This suite needs no stored baseline files. You must still run MATLAB to verify these checks locally. This page does not claim any new numerical result.
