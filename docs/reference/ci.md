# Continuous integration

The five workflows below run on pushes and pull requests to `main` and `water`. They also support manual dispatch. Documentation checks read Markdown sources directly. Website building and deployment are dormant.

| Workflow | File | What it checks | Local equivalent |
|---|---|---|---|
| IMR test suites | `.github/workflows/tests.yml` | Separate MATLAB jobs on R2023a; any failure fails the workflow | [Test commands](../verification/tests.md#run-the-suites) |
| Neatness Check | `.github/workflows/lint.yml` | `checkcode` on `src/`, `tests/` and `toolchain/`, R2023a | `matlab -batch "addpath('toolchain'); s_lint"` from root |
| Format check | `.github/workflows/format.yml` | Runs `s_format` and checks for a Git diff; default MATLAB release from setup action | `matlab -batch "cd('toolchain'); s_format"`, then `git diff --exit-code` |
| Code Spellcheck | `.github/workflows/spellcheck.yml` | Code and documentation spelling | Commands below, from root |
| Documentation checks | `.github/workflows/docs.yml` | Relative Markdown destinations, math delimiters and LaTeX compilation | `bash toolchain/docs/check.sh` |

The lint, format and spelling jobs currently use `continue-on-error: true`. Read their logs even if the overall pull request is green. The tests workflow uses a matrix with `fail-fast: false` so every suite gets a result.

## Spelling

```bash
codespell toolchain src tests imr_setup.m --ignore-words-list="te,nd,Fom,Foh,fom" --skip="*.git,*.md,*.yml,*.json,*.mat" --quiet-level=2
codespell docs README.md CONTRIBUTING.md benchmarks/integrators/README.md --ignore-words-list="te,nd,Fom,Foh,fom" --skip="*.mat,*.jfif,*.svg,*.pdf" --quiet-level=2
```

`te`, `nd`, `Fom` and `Foh` are accepted variable names. A clean run prints nothing.

## MATLAB licensing and scope

The existing MATLAB jobs rely on MathWorks GitHub Actions activation for public repositories. A private repository needs its `MLM_LICENSE_TOKEN` GitHub Actions secret. The workflows fail explicitly if this secret is missing. Enter credentials in GitHub settings. Never enter credentials in documentation or chat.

R2023a numerical CI cannot exercise CVODES. CVODES needs R2024a or a later release. The permutation case generator excludes it on older releases. A passing R2023a job is not proof that CVODES is correct. Regression also needs its [approved baseline files](../verification/regression.md#baseline-availability).

Documentation checking does not require MATLAB. This documentation cleanup does not claim any numerical MATLAB or Octave execution.

## Source documentation checks

Read the [documentation index](../index.md) on GitHub. GitHub renders its Markdown, LaTeX equations and Mermaid diagrams. IMRv2 bundles no JavaScript source, Node.js package manifest or browser-testing dependency.

The checker requires Bash 4.4+, Git, Pandoc 3+, jq and pdfLaTeX with the AMS math packages. From the checkout root on Ubuntu:

```bash
sudo apt-get install --no-install-recommends -y pandoc jq texlive-latex-base texlive-latex-recommended texlive-fonts-recommended
bash toolchain/docs/test.sh
bash toolchain/docs/check.sh
```

`check.sh` reads repository Markdown files, including newly added nonignored files. Pandoc's GitHub-flavored Markdown parser extracts math and link/image destinations. The parser excludes inline and fenced code examples. Relative file destinations must exist. The checker does not check external URLs or anchor fragments.

pdfLaTeX compiles every parsed equation with `amsmath` and `amssymb`. Undefined commands and malformed TeX fail the check. An unparsed dollar sign also fails the check. This catches math that a broken delimiter would otherwise leave as prose. Put literal dollar signs in inline or fenced code. The error message identifies the source page. For a TeX error, the message also identifies the equation number.

The checker removes the generated TeX, PDF and log files on exit. A run with no equations fails. `test.sh` verifies code exclusion and encoded relative links. It also proves that each of these cases fails the check: an undefined command, malformed braces, an unmatched delimiter, a missing file, and an empty equation check.

These checks validate file destinations and equation syntax. They do not validate GitHub anchor slugs, browser appearance or scientific correctness. Preview changed equations, links and diagrams on GitHub before you merge. Some MathJax-specific commands might not work in LaTeX. Use the syntax both renderers support.

## Reactivating MkDocs

`mkdocs.yml` and `docs/requirements.txt` retain the Material theme, navigation and optional website dependencies. No CI job installs them or builds `site/`. The Markdown index is also a usable future home page.

When someone requests the website:

1. Review the retained navigation.
2. Install the optional requirements in a virtual environment.
3. Choose and configure an equation renderer. The retained `arithmatex` markup configuration alone does not render equations. This cleanup removed the previous custom JavaScript and CDN entries.
4. Run `mkdocs build --strict`.
5. Inspect equations, Mermaid diagrams and links with `mkdocs serve`. Check the navigation and the light/dark themes.
6. Add a website-build CI job only after those checks pass. Deployment is a separate decision. This cleanup does not configure GitHub Pages.

Keep source-documentation checks active independently of future website tooling. Simulation code and numerical tests remain MATLAB-only.
