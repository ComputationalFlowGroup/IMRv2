# `s_test_suite_compatibility.m` Documentation

## Overview
`s_test_suite_compatibility.m` verifies that the finite difference (`f_imr_nodal`) and spectral (`f_imr_spectral`) solvers produce compatible results under various model configurations, including thermal effects.

## Major Sections

1. **Initialization**
   - Clears workspace (`clc`, `clear`, `close`).
   - Adds the paths `../toolchain/`, `../src/` and `../tests/`, so it runs from `toolchain/`.
   - Loads test identifiers from `file_ids.mat`.

2. **Preallocation**
   - Computes `num_tests` as the product of parameter combinations.
   - Initializes `errors_fd`, `errors_sp`, and `failed_tests` arrays.

3. **Parameter Setup**
   - Defines parameter vectors:  
     - `muvec`: viscosity values  
     - `Gvec`: shear modulus values  
     - `alphaxvec`: strain-softening exponents  
     - `lambda1vec`: relaxation times  
   - Sets the error threshold `threshold = 1e-4`, `vapor = 1`, `masstrans = collapse = 0`, `R0 = 50e-6` and `Req = R0/10`, in SI units converted with `f_dim2nondim`.

4. **Nested Test Loop**
   - Loops over combinations of:  
     `radial` (1–4), `bubtherm` (0/1), `medtherm` (0/1), `stress` (0–5), skipping `medtherm = 1` with radial 3–4 (their stored results are kept but not compared),  
     `mu`, `G`, `alphax`, `lambda1` indices.
   - For each combination:
     1. Constructs `varin` parameter-value pairs.
     2. Loads reference results: `Rf` and `Rs` from test MAT files.
     3. Runs `f_imr_nodal` and compares output `Rf_test`; computes L2 error.
     4. Runs `f_imr_spectral` and compares output `Rs_test`; computes L2 error.
     5. Records failures where error exceeds `threshold`.

5. **Summary and Exit**
   - Truncates and cleans `failed_tests`.
   - If all tests passed, prints success and calls `exit(0)`.
   - Otherwise, prints failed test indices and raises an error.

## Inputs

| Name             | Description                                    |
|------------------|------------------------------------------------|
| `file_ids.mat`   | MAT-file in `toolchain/` defining `ids` (list of test case IDs)|
| `f_imr_nodal`       | Finite difference solver function              |
| `f_imr_spectral` | Spectral solver function                       |

## Outputs

- **Console output**:
  - L2 norm error for each finite-difference and spectral test.
  - Summary line: pass or fail with indices.
- **Exit code**:
  - `0` if all tests pass.
  - Error thrown (non-zero) if any test fails.

## Notes

- `num_tests` covers the valid combinations of radial, thermal, stress, and material parameters; `medtherm = 1` without `bubtherm` is skipped.
- Tests use reduced mesh sizes (`Nt=70, Mt=70` for FD; `Nt=12, Mt=12` for spectral) for speed.
