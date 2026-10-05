# `s_test_suite_regression.m` and `s_generate_references.m`

## Overview

`s_test_suite_regression.m` runs the regression cases in the three solver modes and compares `R(t)` with the stored results in `tests/reference/<mode>.mat`. It catches unintended changes of the results.

## Cases

From [`f_test_cases('regression')`](f_test_cases.md). Only the combinations the compatibility matrix allows are kept: 230 cases each for `nodal_fd` and `nodal_cheb`, and 155 for `spectral`.

- **Radial, thermal and stress models, for two materials.** The cases are:

    - Radial models 1–6.
    - Five thermal configurations: polytropic, bubble, bubble and medium, bubble with mass transfer, and bubble and medium with mass transfer.
    - Stress models 0–5.

    The two materials are a soft, weakly viscous one (`mu = 1e-4 Pa s`, `G = 100 Pa`, `alphax = 1e-3`, `lambda1 = 1e-7 s`) and a stiff, viscous one (`mu = 0.1 Pa s`, `G = 2.5 kPa`, `alphax = 1`, `lambda1 = 1 ms`).

    The runs are:

    - Without mass transfer: a 50 µm bubble released from four times its equilibrium radius and followed for 12 µs.
    - The polytropic bubble: released from twice its equilibrium radius. It loses no heat, so its collapse is far sharper.
    - With mass transfer: a 200 µm bubble released from four times its equilibrium radius and followed for 30 µs.

    These runs exercise every term without being so violent that the result mostly tests the integrator. From ten times the equilibrium radius, as the earlier suites used, about half of the cases stop at the first collapse.

- **Other features**, on Keller–Miksis with bubble and medium heat transfer:

    - The seven non-Newtonian viscosity models.
    - The Gaussian, histotripsy and Heaviside waveforms.
    - The UCM / Oldroyd-B stress field (`spectral = 1`, spectral solver only).
    - Collapse with stress 0, 1, 3 and 4, with the soft material, up to the first collapse and rebound (20 µs). This is nodal solver only, since collapse needs mass transfer.

## Check

A case fails if:

- its run fails;
- it has no stored result;
- its output size differs from the stored one; or
- `norm(R./R_ref - 1, 2) > 1e-4` over its 100 output times.

The suite also notes stored results that no case uses any more.

## Generating the stored results

`s_generate_references.m` runs every regression case and writes `tests/reference/nodal_fd.mat`, `nodal_cheb.mat` and `spectral.mat`. Each file holds:

- `keys`: the case keys.
- `R`: the results.
- `info`: the date and program version.

Regenerate them only when a change of results is intended and agreed (see Contributing).

- **Overwrite guard.** It refuses to overwrite existing files unless `overwrite = true` is set in the script.
- **Parallel runs.** It runs the cases with `parfor`, which uses the default parallel pool when the Parallel Computing Toolbox is available and runs in series otherwise.
- **All or nothing.** If any case fails, it lists the failures and writes nothing. Otherwise it writes all files to a temporary folder and then moves them into place.
- **Cross-mode report.** For the cases the modes share, it reports the largest relative difference between them, and lists any case above `1e-2`. This is a sanity check on the new results, not a test: the grids differ, and violent runs may differ more.

A new case added to `f_test_cases` fails the suite until its result is stored.

## Agreement between the modes

In the cross-mode report, Chebyshev and the spectral solver agree with finite differences to about 2% on the cases without mass transfer. The mass-transfer cases differ by up to 11% between the two nodal schemes. At `Nt = Mt = 20` (Chebyshev) and 70 (finite differences), neither nodal scheme is resolved for mass transfer: the convergence study in the `f_imr_nodal` documentation needs about 40 Chebyshev points. The suite checks that results are reproduced, not that they are accurate, so these grids are adequate.
