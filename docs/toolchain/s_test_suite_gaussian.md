# `s_test_suite_gaussian.m` Documentation

## Overview
`s_test_suite_gaussian.m` checks that a bubble driven by a Gaussian pressure pulse returns to its equilibrium radius, in both solvers. The CI workflow `test_equilibrium.yml` runs it after the equilibrium suite.

## Setup
- Adds the paths `../toolchain/` and `../src/`, so it runs from `toolchain/`.
- SI inputs converted with `f_dim2nondim`: `Req = R0 = 10e-6`, `T8 = 298.15`, `mu = 0.5`, a Gaussian pulse (`wave_type = 1`) with amplitude `pA = 1e5`, width `tw = 1e-4` and delay `dt = 3e-4`, over `tfin = 2e-3` with 1000 output times.
- `masstrans = collapse = 0`, integrator `ode23tb` (`method = 23`).

## Test loop
- Loops over `radial` 1–4, `bubtherm` 0/1, `medtherm` 0/1, `vapor` 0/1 and `stress` 1–5, skipping `medtherm = 1` without `bubtherm`.
- Runs `f_imr_nodal` with `Nt = Mt = 50` and `f_imr_spectral` with `Nt = Mt = 12`.
- The error is the departure of the final radius from equilibrium: `|1 - R(end)|` for the nodal solver (finite difference scheme) and `|R(1) - R(end)|` for the spectral solver.
- A test fails when the error exceeds `threshold = 1e-2`.

## Outputs
- Console output with the error of each test and a pass or fail summary.
- Exit code `0` if all tests pass, `1` otherwise.
