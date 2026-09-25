# Contributing to IMR

Thank you for helping improve IMR. This guide walks you through a contribution from start to finish: fork the repository, make your change, run the same checks that the continuous integration (CI) runs, and open a pull request.

IMR is released under the [GNU General Public License](https://github.com/InertialMicrocavitationRheometry/IMRv2/blob/main/LICENSE), version 3 or (at your option) any later version (`GPL-3.0-or-later`). By submitting a contribution, you agree that it is licensed under the same terms, and you confirm that you have the right to contribute it. For example, do not copy code from sources whose license is incompatible with the GPL.

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

1. Fork [InertialMicrocavitationRheometry/IMRv2](https://github.com/InertialMicrocavitationRheometry/IMRv2) with the **Fork** button on GitHub.
2. Clone your fork and add the main repository as `upstream`:

    ```bash
    git clone https://github.com/<your-username>/IMRv2.git
    cd IMRv2
    git remote add upstream https://github.com/InertialMicrocavitationRheometry/IMRv2.git
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

GNU Octave runs most of the solvers, but the test suites need MATLAB. Octave has no `ode23tb`, the default integrator, and differs from MATLAB in a few other places; see the notes at the end of [`docs/compatibility.md`](docs/compatibility.md).

## 2. Make your change

### Where things live

| Folder | Contents |
|---|---|
| `src/` | The solvers (`f_imr_nodal`, `f_imr_spectral`) and the physics functions they call |
| `toolchain/` | The test suites, the formatter, the linter, the test generators and `f_dim2nondim` |
| `tests/` | Stored reference results (`.mat`) that the test suites compare against |
| `docs/src/`, `docs/toolchain/` | One Markdown page per function or script |
| `docs/compatibility.md` | Which features run together in each solver |
| `examples/`, `benchmarks/` | Example scripts and scaling benchmarks |

### Code conventions

- **Base MATLAB only.** Code in `src/` must not need any toolbox, so do not add calls to toolbox functions.
- **Dimensionless inputs.** The solvers take and return dimensionless quantities only. Dimensional values are converted beforehand with `toolchain/f_dim2nondim.m`.
- **Adding an input.** Follow the path every input takes:
    1. Give it a default in `src/default_case.m`, and add it to the struct `d` at the end of that function.
    2. In `src/f_call_params.m`, copy it out of the defaults (`x = d.x;`) and read it in the name/value parser. If its input name differs from its name in `default_case`, add the pair to the `renamed` list in `src/f_imr_inputs.m`.
    3. Validate it in `src/f_checkers.m`, raising an `INPUT ERROR` for invalid values.
    4. Document it in `docs/src/default_case.md`.
- **Option vectors: append, never insert.** `f_call_params` passes settings to the solvers in positional vectors, such as `eqns_opts`, `solve_opts`, `init_opts` and `sigma_opts`. Add a new entry at the end of a vector. Inserting one in the middle shifts every later position and silently breaks the code that unpacks it.
- **File headers.** Start every `.m` file with the same header as the existing files: the file name, a one-line summary, the license block, then the description. The license block is a separate comment block, so it stays out of MATLAB's `help` text, which shows only the first block.

    ```matlab
    % file f_example.m
    % brief contains function f_example

    % Copyright (C) The IMR developers (see README.md)
    % SPDX-License-Identifier: GPL-3.0-or-later

    % brief This function computes ... (one or two sentences on what it does)
    function out = f_example(in)
    ```

    Shell scripts carry the same two license lines as `#` comments after the `#!` line.

- **Avoid semicolons inside strings and comments.** To keep one statement per line, `s_format` splits every line that contains more than one `;` at the first one. It does not know about strings or comments, so a line such as `error(['... not available; use radial = 1']);` is split inside the string, which breaks the code. Use a comma or a full stop instead.

- **Comments.** Explain why the code does something, especially for numerical choices such as tolerances, brackets and quadrature rules.
- **Keep changes focused.** One topic per pull request. Unrelated clean-ups go in their own pull request.

### Documentation

- When you change a function's behavior, inputs or outputs, update its page in `docs/src/` or `docs/toolchain/`. A new function gets a new page, modeled on the existing ones.
- If your change affects which features run together, update [`docs/compatibility.md`](docs/compatibility.md) and the rules in `f_checkers.m` together.
- User-facing changes, such as a new model or option, also belong in the feature list of the [README](https://github.com/InertialMicrocavitationRheometry/IMRv2/blob/main/README.md).

### Documentation website

The pages in `docs/` are published as a website with [MkDocs](https://www.mkdocs.org) and the [Material](https://squidfunk.github.io/mkdocs-material/) theme, configured in `mkdocs.yml`. To preview your changes, install the packages once and serve the site locally:

```bash
pip install -r docs/requirements.txt
mkdocs serve
```

Then open <http://127.0.0.1:8000>; the page reloads as you edit.

- **New pages** must be added to `nav` in `mkdocs.yml`.
- **Home, Getting started and About** are built from the sections of `README.md` between the `<!-- site:... -->` markers, and **Contributing** from this file, so each text lives in one place. Links in those parts must work both on GitHub and on the website: link to pages in `docs/` with relative links (such as `docs/compatibility.md`), and to other repository files, such as `LICENSE`, with full GitHub URLs.
- **Lists** need a blank line before them, and content nested in a list item (sub-items, code blocks, paragraphs) is indented by 4 spaces. GitHub accepts less, but the website's Markdown engine does not.
- **Equations** can be written as `$...$` inline and `$$...$$` for display.

### Stored test results

The files in `tests/*.mat` are the reference that the test suites compare against. Do not edit or regenerate them as part of an ordinary change. If your change deliberately changes the results, for example by fixing a physics bug:

1. Explain in the pull request which results change, why, and by how much.
2. Regenerate them only after the maintainers agree, with `toolchain/s_generate_compatibility.m` or `toolchain/s_generate_masstransfer.m`.

## 3. Run the checks

Run these before you push. They are the same checks CI runs on every pull request, so a clean local run means a clean CI run.

| Check | What it does | Command |
|---|---|---|
| Spelling | `codespell` on `src/` and `toolchain/` | Shell, from the repository root (see below) |
| Formatting | Re-indents every `.m` file in `src/`, `toolchain/`, `docs/` and `benchmarks/` | `cd toolchain; s_format` in MATLAB |
| Lint | MATLAB's `checkcode` on `src/` and `toolchain/` | `matlab -batch "addpath('toolchain'); s_lint"` |
| Tests | The five test suites | See [Test suites](#test-suites) |
| Documentation website | Builds the site and fails on broken links or anchors, and on pages missing from the navigation | `mkdocs build --strict` (see [Documentation website](#documentation-website)) |

### Spelling

From the repository root:

```bash
codespell toolchain src --ignore-words-list="te,nd,Fom,Foh" --skip="*.git,*.md,*.yml,*.json" --quiet-level=2
```

It prints nothing when there are no misspellings. To also check the documentation you changed:

```bash
codespell docs README.md CONTRIBUTING.md --ignore-words-list="te,nd,Fom,Foh,fom" --quiet-level=2
```

`te`, `nd`, `Fom` and `Foh` are variable names, such as `Fom` for the mass Fourier number, not misspellings. If `codespell` flags a new variable name, rename the variable if you can. Otherwise add the name to `--ignore-words-list` in `.github/workflows/spellcheck.yml` and in this guide.

### Formatting

In MATLAB, from the repository root:

```matlab
cd toolchain
s_format
```

`s_format` rewrites the files in place. Then check what it changed and commit the result:

```bash
git diff
```

The CI format check runs `s_format` and fails if that leaves any change. So "formatting passes" means running `s_format` produces no diff.

### Lint

`s_lint` runs MATLAB's `checkcode` on `src/` and `toolchain/`, skipping `f_multirun.m`. It ends with `exit`, which closes MATLAB, so run it from a shell:

```bash
matlab -batch "addpath('toolchain'); s_lint"
```

Fix the warnings in the files you changed. If a warning is intentional, suppress it on that line with a `%#ok<ID>` comment and say why in a comment.

### Test suites

Each suite runs from `toolchain/`:

| Suite | Checks |
|---|---|
| `s_test_suite_compatibility` | Both solvers against stored results, over radial, thermal and stress models |
| `s_test_suite_equilibrium` | A bubble at equilibrium stays there, in both solvers |
| `s_test_suite_gaussian` | A bubble returns to equilibrium after a Gaussian pressure pulse |
| `s_test_suite_masstrans` | The nodal solver with mass transfer against stored results |
| `s_test_suite_diffmat` | The differentiation operators of `f_diff_mat` |

Most suites end with `exit`, so run each one in its own MATLAB process:

```bash
cd toolchain
for s in s_test_suite_diffmat s_test_suite_equilibrium s_test_suite_gaussian \
         s_test_suite_compatibility s_test_suite_masstrans; do
    matlab -batch "$s" || echo "FAILED: $s"
done
```

Each suite prints `✅ All tests PASSED.` or lists the failing cases. `s_test_suite_masstrans` does not set an exit code, so read its last line rather than relying on the exit status.

If you add a model or an option, add a case that exercises it, or say in the pull request how you tested it.

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

1. On GitHub, open a pull request from `my-change` in your fork to `main` in `InertialMicrocavitationRheometry/IMRv2`.
2. In the description, cover:
    - **What** the change does and **why** it is needed. Link the issue if there is one.
    - **How you tested it:** the checks and suites you ran, and any extra runs, such as convergence tests or comparisons with the previous version.
    - **Changes in results:** whether any solver output changes, and by how much. A change that leaves `R(t)` bit-identical is worth stating as such.
    - **Documentation** you updated.
3. Check the CI results on the pull request. The jobs are:

    | Workflow | Runs |
    |---|---|
    | Code Spellcheck | `codespell` on `src/` and `toolchain/` |
    | Format check | `s_format`, then fails on any diff |
    | Neatness Check | `s_lint` |
    | IMR test suites | The five test suites, one job each |
    | Documentation website | `mkdocs build --strict`; on `main`, also publishes the site |

    The spelling, format and lint workflows are marked `continue-on-error`, so a failure there may not turn the pull request red; open their logs and confirm they passed. The test suites workflow fails if any suite fails, except `s_test_suite_masstrans`, which does not set an exit code: read the last line of its log.

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
