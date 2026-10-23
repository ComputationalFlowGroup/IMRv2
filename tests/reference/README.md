# Approved golden simulation histories

The regression suite expects `nodal_fd.mat`, `nodal_cheb.mat` and `spectral.mat` here. Each file must contain matching `keys`, radius histories `R` and generation metadata `info` for the cases in `f_test_cases('regression')`.

These `.mat` files are absent from the current checkout and will be recreated locally in MATLAB. The unused older histories have been removed.

Use `run_validation('generate')` from `tests/`, or the existing [reference generator script](../../toolchain/s_generate_references.m), to create the keyed histories through one shared [helper](../helpers/generate_references.m). Review the results before accepting them as baselines. Ordinary test runs only read these files. See the [runner modes](../README.md#one-command-for-all-tests) and [generation command, safeguards and baseline availability](../../docs/validation/regression.md#generating-the-stored-results).
