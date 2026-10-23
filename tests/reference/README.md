# Approved golden simulation histories

The regression suite expects `nodal_fd.mat`, `nodal_cheb.mat` and `spectral.mat` here. Each file must contain matching `keys`, radius histories `R` and generation metadata `info` for the cases in `f_test_cases('regression')`.

These approved `.mat` files are absent from the current checkout. This README identifies their destination without supplying fabricated or newly generated baselines. Obtain approved files before expecting the regression suite to pass.

The separate [reference generator](../../toolchain/s_generate_references.m) can generate new histories when that action has been explicitly agreed. Routine suites only read expected results. See [generation safeguards and baseline availability](../../docs/validation/regression.md).
