# Continuous integration

The five workflows below run on pushes and pull requests to `main`, and support manual dispatch. Documentation builds are local artifacts only; no GitHub Pages upload or deployment is configured.

| Workflow | File | What it checks | Local equivalent |
|---|---|---|---|
| IMR test suites | `.github/workflows/tests.yml` | Separate MATLAB jobs on R2023a; any failure fails the workflow | [Test commands](../validation/tests.md#run-the-suites) |
| Neatness Check | `.github/workflows/lint.yml` | `checkcode` on `src/` and `toolchain/`, R2023a | `matlab -batch "addpath('toolchain'); s_lint"` from root |
| Format check | `.github/workflows/format.yml` | Runs `s_format` and checks for a Git diff; default MATLAB release from setup action | `matlab -batch "cd('toolchain'); s_format"`, then `git diff --exit-code` |
| Code Spellcheck | `.github/workflows/spellcheck.yml` | Code and documentation spelling | Commands below, from root |
| Documentation build | `.github/workflows/docs.yml` | Strict MkDocs build, MathJax compilation and Chromium rendering/navigation | Commands below |

The lint, format and spelling jobs currently use `continue-on-error: true`. Read their logs even if the overall pull request is green. The tests workflow uses a matrix with `fail-fast: false` so every suite gets a result.

## Spelling

```bash
codespell toolchain src tests/fixtures imr_setup.m --ignore-words-list="te,nd,Fom,Foh,fom" --skip="*.git,*.md,*.yml,*.json" --quiet-level=2
codespell docs README.md CONTRIBUTING.md benchmarks/integrators/README.md --ignore-words-list="te,nd,Fom,Foh,fom" --skip="*.mat,*.jfif,*.svg,*.pdf" --quiet-level=2
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

## Equation rendering

From the checkout root, after the strict build:

```bash
npm ci --ignore-scripts --prefix toolchain/docs
cd toolchain/docs
npx playwright install chromium
cd ../..
node toolchain/docs/check.cjs site
```

The pinned documentation-only JavaScript dependencies compile every HTML `arithmatex` expression with MathJax 3 and fail on undefined commands. A deliberate invalid-command fixture checks that failure detection works. Chromium then checks radial, thermal and wall equations on initial load and navigation. The browser uses the pinned local MathJax bundle, so the check does not depend on CDN availability. A local system Chromium may be selected with `IMR_CHROMIUM=/path/to/chromium`. No JavaScript tooling is required to run IMRv2 in MATLAB.

These checks validate the configured MkDocs renderer. GitHub Markdown uses its own renderer; inspect the source Markdown preview when changing delimiters or advanced TeX. No GitHub Pages build or deployment is configured.
