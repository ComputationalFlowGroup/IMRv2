# Contributing to IMR

This guide walks you through a contribution from start to finish: fork the repository, make your change, run the same checks that the continuous integration (CI) runs, and open a pull request.

IMR is released under the [GNU General Public License](https://github.com/mrodrig6/imrv2/blob/main/LICENSE), version 3 or (at your option) any later version (`GPL-3.0-or-later`). By submitting a contribution, you agree that it is licensed under the same terms, and you confirm that you have the right to contribute it. For example, do not copy code from sources whose license is incompatible with the GPL.

If you are planning a large change, such as a new physical model, a new solver option or a change to stored results, open an issue first or ask on the IMR Slack workspace. Agreeing on the approach early saves rework later.

## Contents

1. [Set up your fork](#1-set-up-your-fork)
2. [Make your change](#2-make-your-change)
3. [Run the checks](#3-run-the-checks)
4. [Commit and push](#4-commit-and-push)
5. [Open a pull request](#5-open-a-pull-request)
6. [Review and updates](#6-review-and-updates)
7. [Reporting bugs](#reporting-bugs)

## 1. Set up your fork

1. Fork [mrodrig6/imrv2](https://github.com/mrodrig6/imrv2) with the **Fork** button on GitHub.
2. Clone your fork and add this repository as `upstream`:

    ```bash
    git clone https://github.com/<your-username>/IMRv2.git
    cd IMRv2
    git remote add upstream https://github.com/mrodrig6/imrv2.git
    ```

3. Create a branch for your change from an up-to-date `main`:

    ```bash
    git fetch upstream
    git checkout -b my-change upstream/main
    ```

    Use one branch per change, and give it a descriptive name, such as `fix-wall-temperature` or `add-cross-viscosity`.

### Tools you need

| Tool | Used for |
|---|---|
| MATLAB (R2021a or newer; CI uses R2023a) | Running the solvers, the formatter, the linter and the test suites |
| Python with [`codespell`](https://github.com/codespell-project/codespell) (`pip install codespell`) | The spelling check |
| Git | Version control |

## 2. Make your change

### Where things live

| Folder | Contents |
|---|---|
| `src/` | Everything a user calls: the solvers (`imr.solve_nodal`, `imr.solve_spectral`), their inputs (`imr.default_inputs`), sweeps (`imr.sweep`), unit conversion (`imr.nondimensionalize`), and the physics functions the solvers call |
| `toolchain/` | Developer tooling only: the test suites, the formatter, the linter and the test generators |
| `tests/` | Stored results: `tests/reference/` for the regression suite. The older `tests/*.mat` files are no longer used by any suite |
| `docs/api/`, `docs/reference/functions/`, `docs/validation/` | One Markdown page per function or script |
| `docs/reference/compatibility.md` | Which features run together in each solver |
| `examples/`, `benchmarks/` | Example scripts (including the `f_multirun` study and its cluster job), scaling benchmarks, and the integrator comparison (`benchmarks/integrators`) |

### Code conventions

Use lowercase `snake_case` without an `f_` prefix for new functions. Place public calls in `+imr` and implementations in their responsibility package. The [namespace guide](namespaces.md) covers file/function naming, compatibility and runtime checks.

- **Base MATLAB only.** Code in `src/` must not need any toolbox, so do not add calls to toolbox functions.
- **Dimensionless inputs.** The solvers take and return dimensionless quantities only. Dimensional values are converted beforehand with `src/+imr/nondimensionalize.m`.
- **Adding an input.** Follow the [extension guide](../api/extending.md): add the default in `imr.configuration.defaults`, its public/internal mapping and discrete rule in `imr.configuration.schema`, and its named prepared-case field in `imr.configuration.prepare_case`. Keep coupled physical restrictions in `imr.configuration.validate_case`, document the input, and add an independent numerical check.
- **Named internal contracts.** Both solvers consume `imr.configuration.prepare_case` and `imr.solvers.build`. Use `imr.solvers.state_layout` blocks for initialization, RHS, Jacobians and reconstruction. Do not add new positional unpacking to the solvers. `imr.configuration.legacy_options` is a compatibility serializer; preserve the order of its existing vector outputs for legacy callers.
- **File headers.** Start every `.m` file with the same header as the existing files: the file name, a one-line summary, the license block, then the description. The license block is a separate comment block, so it stays out of MATLAB's `help` text, which shows only the first block.

    ```matlab
    % file example_quantity.m
    % brief contains function example_quantity

    % Copyright (C) The IMR developers (see README.md)
    % SPDX-License-Identifier: GPL-3.0-or-later

    % brief This function computes ... (one or two sentences on what it does)
    function out = example_quantity(in)
    ```

    Shell scripts carry the same two license lines as `#` comments after the `#!` line.

- **Avoid semicolons inside strings and comments.** To keep one statement per line, `s_format` splits every line that contains more than one `;` at the first one. It does not know about strings or comments, so a line such as `error(['... not available; use radial = 1']);` is split inside the string, which breaks the code. Use a comma or a full stop instead.

- **Comments.** Explain why the code does something, especially for numerical choices such as tolerances, brackets and quadrature rules.
- **Keep changes focused.** One topic per pull request. Unrelated clean-ups go in their own pull request.

### Documentation

- When you change a function's behavior, inputs or outputs, update its page in the appropriate topic under `docs/`. A new function gets a new page, modeled on the existing ones.
- If your change affects which features run together, update [`docs/reference/compatibility.md`](compatibility.md) and the rules in `imr.configuration.validate_case` together.
- User-facing changes, such as a new model or option, also belong in the feature list of the [README](https://github.com/mrodrig6/imrv2/blob/main/README.md).

### Documentation checks

Source-checking and spelling commands live in [CI](ci.md). Read and review the Markdown on GitHub. The Material site configuration in `mkdocs.yml` is dormant; CI does not build or deploy a website.

- Add new pages to the matching section folder and navigation.
- Keep `docs/index.md` as a normal Markdown topic index and link to the README for setup. Do not add custom Python source or duplicate README prose in other pages.
- Use relative links between Markdown pages so they also work on GitHub. Link repository code files explicitly when useful.
- Start a page with a short scope statement. Keep tables for parameters and contracts; link canonical facts rather than repeating them.
- Give figures descriptive alt text and keep their source/rebuild instructions beside them.
- Lists need a blank line before them. Indent nested content by four spaces for the Markdown engine.
- Use `$...$` and `$$...$$` for equations.

### Stored test results

The files in `tests/reference/` are the stored results that `s_test_suite_regression` compares against, one file per solver mode with each result under its case key (such as `r2_bubble+medium_s3_mat1`). Do not edit or regenerate them as part of an ordinary change. If your change deliberately changes the results, for example by fixing a physics bug:

1. Explain in the pull request which results change, why, and by how much.
2. Regenerate them only after the maintainers agree, with `toolchain/s_generate_references.m` in MATLAB (set `overwrite = true` in the script). It writes nothing unless every case succeeds.

A new case added to `toolchain/f_test_cases.m` fails the regression suite until its result is stored.

## 3. Run the checks

Follow [CI and local equivalents](ci.md) for spelling, MATLAB lint, formatting and source-documentation checks. Run the [numerical suites](../validation/tests.md#run-the-suites) from `toolchain/` in separate MATLAB processes and preserve failure status.

The regression suite requires approved references that are absent from the audited checkout. See [baseline availability](../validation/regression.md#baseline-availability). Do not fabricate or silently replace them.

`s_format` rewrites MATLAB files, including files under `docs/` and `benchmarks/`, and skips itself. Inspect every resulting diff before committing. It does not format `examples/`. `s_lint` calls `exit`, so run it in a separate MATLAB process.

The user runs MATLAB validation locally for this documentation change; no numerical pass is claimed. A successful local run does not guarantee CI licensing or a different MATLAB release will succeed.

## 4. Commit and push

- Write commit messages with a short summary line (about 70 characters) in the imperative, such as "Fix the wall temperature output". Explain what changed and why in the body.
- Commit the formatter's changes together with your own.
- Do not commit generated files, local data or editor files (`.gitignore` covers the common ones).

```bash
git add <files>
git commit
git push -u origin my-change
```

## 5. Open a pull request

1. On GitHub, open a pull request from `my-change` in your fork to `main` in `mrodrig6/imrv2`.
2. In the description, cover:
    - **What** the change does and **why** it is needed. Link the issue if there is one.
    - **How you tested it:** the checks and suites you ran, and any extra runs, such as convergence tests or comparisons with the previous version.
    - **Changes in results:** whether any solver output changes, and by how much. A change that leaves `R(t)` bit-identical is worth stating as such.
    - **Documentation** you updated.
3. Check each [CI workflow](ci.md), including logs for jobs marked `continue-on-error`. Documentation checks must pass without publishing a site.

## 6. Review and updates

- A maintainer will review the pull request. Answer each comment, either with a change or with the reason you are keeping the code as it is.
- Push fixes to the same branch; the pull request updates automatically. Rerun the checks for anything you changed.
- If `main` has moved on and your branch conflicts, bring it up to date:

    ```bash
    git fetch upstream
    git merge upstream/main
    ```

    Resolve the conflicts, rerun the checks and push.

- Once the pull request is approved and CI passes, a maintainer merges it.

## Reporting bugs

Open an issue with:

- the solver call (all name/value inputs) that shows the problem,
- what you expected and what happened, including the full error message,
- your MATLAB version and operating system.

A short script that reproduces the problem is the most helpful thing you can include.

For questions, contact [Mauro Rodriguez](mailto:mrodri97@brown.edu) or ask on the IMR Slack workspace.
