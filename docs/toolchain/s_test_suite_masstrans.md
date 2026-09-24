# `s_test_suite_masstrans.m` Documentation

## Overview
`s_test_suite_masstrans.m` validates the finite difference solver’s mass transfer implementation under various model settings.

## Major Sections

1. **Initialization**
   - Clears workspace.
   - Adds the paths `../toolchain/`, `../src/` and `../tests/`, so it runs from `toolchain/`.
   - Loads `file_ids.mat`.

2. **Preallocation**
   - Computes `num_tests`.
   - Initializes `errors_fd` and `failed_tests`.
   - Sets `threshold = 5e-5`.

3. **Test Parameter Setup**
   - Defines `tvector` for simulation.
   - Fixes `masstrans = 1`, `vapor = 1`, `collapse = 0`.
   - Sets the bubble geometry `R0 = 200e-6`, `Req = 35e-6` (SI units, converted with `f_dim2nondim`).
   - Loops over `radial` 1–4 and `stress` 0–5 with `bubtherm = medtherm = 1`.

4. **Loop Over Test Cases**
   - Iterates over all combinations.
   - Builds `varin` parameters.
   - Loads reference `Rm` from test files.
   - Executes `f_imr_fd(...,'Nt',70,'Mt',70)`, computes L2 error.
   - Records failures above `threshold`.

5. **Summary**
   - Trims `failed_tests`.
   - Prints pass/fail summary (no automatic exit).

## Inputs

| Name            | Description                          |
|-----------------|--------------------------------------|
| `file_ids.mat`  | Test case identifiers                |
| `f_imr_fd`      | Finite difference solver function    |

## Outputs

- **Console output**: L2 norm errors per test and summary.
- **Error flags**: `failed_tests` indices where mass transfer failed.

## Notes

- Script does not exit on failure to allow further analysis.
- Uses `Nt = Mt = 70`.
