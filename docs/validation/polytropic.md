# Constant-vapor polytropic validation

`tests/suites/s_test_suite_polytropic.m` checks the common gas/vapor pressure law in both solvers. From `toolchain/`, run `matlab -batch "imr_test_setup; s_test_suite_polytropic"`.

For fixed ambient vapor pressure, total bubble pressure is

$$
P=P_{g0}R^{-3\kappa}+P_v,\qquad
\dot P=-3\kappa P_{g0}R^{-3\kappa}\frac{\dot R}{R}.
$$

$P_{g0}$ is initial noncondensable gas pressure, $P_v$ is constant vapor pressure and $\kappa$ belongs to the gas. Only gas pressure is differentiated. The old expression incorrectly used total pressure in the derivative and therefore differentiated a constant vapor contribution.

The suite compares the derivative with finite differences, verifies invariance to a constant vapor offset, and compares integrated pressure against the algebraic pressure law during motion in nodal FD, nodal Chebyshev and spectral modes. It needs no stored baseline files. MATLAB execution remains for local verification; no new numerical result is claimed here.
