# Full and paper distributions

The main development repository keeps all supported physics and optional material add-ons. The `water` branch shares the numerical advancements, documentation and artwork, but omits the mercury-specific profile, example, documentation, test suite and frozen observations. It retains generic EoS families, the six radial models, transport, rheology, Jacobians and validation of common physics.

The add-on boundary consists of:

- `src/+imr/+fluids/mercury.m`
- `examples/s_run_mercury.m`
- `tests/suites/s_test_suite_mercury.m`
- `tests/+imrtest/mercury_reference.m`
- `docs/model/mercury.md` and `docs/validation/mercury.md`
- Small optional-material links and the material-specific CI entry.

No generic source function or water test depends on those files. Water default and polytropic checks remain separate suites. CI also runs on `water` pushes and pull requests targeting that branch.

## Branch versus new repository

A branch without those files still shares earlier Git history with the full repository. Do not push that history to a paper repository intended to contain no mercury information. Instead, seed the new repository from a clean snapshot of the `water` tree, retaining the license and existing attribution.

For example, from the original checkout, with a new empty destination directory:

```bash
git archive water | tar -x -C /path/to/empty-paper-directory
cd /path/to/empty-paper-directory
git init -b main
git add .
git commit -m "Initial paper code snapshot"
git remote add origin https://github.com/ORGANIZATION/REPOSITORY.git
git push -u origin main
```

Create the destination with the intended visibility first. Do not run these initialization commands inside an existing repository containing unrelated work. This exports the current files without bringing source-repository history into the new repository; it is not a history rewrite of the full repository.

## Maintaining both versions

Keep common solver changes in commits that do not mix in material-specific coefficients or documentation. Cherry-pick those commits onto `water`, run the shared tests and strict docs build, and publish the reviewed snapshot to the paper repository. Merge all full-physics changes through the usual review process rather than repeatedly deleting the add-on after unrestricted merges.

The paper branch does not add GitHub Pages or a Python simulation runtime. Numerical suites and examples still require licensed MATLAB; local static checks alone do not establish numerical success.
