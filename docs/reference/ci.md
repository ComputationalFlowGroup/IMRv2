# Continuous integration

The five workflows below run on pushes and pull requests to `main`, and support manual dispatch. Documentation builds are local artifacts only; no GitHub Pages upload or deployment is configured.

| Workflow | File | What it checks | Local equivalent |
|---|---|---|---|
| IMR test suites | `.github/workflows/tests.yml` | Seven separate MATLAB jobs on R2023a; any failure fails the workflow | [Test commands](../validation/tests.md#run-the-suites) |
| Neatness Check | `.github/workflows/lint.yml` | `checkcode` on `src/` and `toolchain/`, R2023a | `matlab -batch "addpath('toolchain'); s_lint"` from root |
| Format check | `.github/workflows/format.yml` | Runs `s_format` and checks for a Git diff; default MATLAB release from setup action | `matlab -batch "cd('toolchain'); s_format"`, then `git diff --exit-code` |
| Code Spellcheck | `.github/workflows/spellcheck.yml` | Code and documentation spelling | Commands below, from root |
| Documentation build | `.github/workflows/docs.yml` | Python 3.12, installs docs requirements, `mkdocs build --strict` | Same build command from root |

The lint, format and spelling jobs currently use `continue-on-error: true`. Read their logs even if the overall pull request is green. The tests workflow uses a matrix with `fail-fast: false` so every suite gets a result.

## Spelling

```bash
codespell toolchain src --ignore-words-list="te,nd,Fom,Foh" --skip="*.git,*.md,*.yml,*.json" --quiet-level=2
codespell docs README.md CONTRIBUTING.md benchmarks/integrators/README.md --ignore-words-list="te,nd,Fom,Foh,fom" --skip="*.mat,*.jfif,*.svg" --quiet-level=2
```

`te`, `nd`, `Fom` and `Foh` are accepted variable names. A clean run prints nothing. Both commands passed during this documentation change.

## MATLAB licensing and scope

The existing MATLAB jobs rely on MathWorks GitHub Actions activation for public repositories. A private repository needs its `MLM_LICENSE_TOKEN` GitHub Actions secret; the workflows fail explicitly if it is missing. Enter credentials in GitHub settings, never documentation or chat.

R2023a numerical CI cannot exercise CVODES, which needs R2024a+. The permutation case generator excludes it on older releases. No MATLAB version matrix is added by this documentation change, and a passing R2023a job is not evidence of CVODES validation.

No MATLAB or Octave was run in the cloud documentation task. The user will run MATLAB checks locally. Regression cannot pass without the [approved baseline files](../validation/regression.md#baseline-availability).

## Strict local build

```bash
python -m pip install -r docs/requirements.txt
mkdocs build --strict
```

The build checks navigation, missing pages, links and anchors. The configured include-markdown plugin reads README as the home page and rewrites relative links. No custom Python source is bundled with IMRv2; simulation and numerical tests run in MATLAB alone. Retain the Material light/dark theme, copy buttons and MathJax configuration when adding pages. Add every new page to navigation.
