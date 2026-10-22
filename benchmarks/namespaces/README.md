# Namespace timing and numerical comparison

Run this benchmark before merging the namespace reorganization. It needs base MATLAB R2021a+ and an unchanged checkout of the previous flat layout. It does not download or install MATLAB or change a repository branch.

From a shell, create the comparison checkout at the last pre-namespace commit:

```bash
git worktree add --detach /tmp/imrv2-before-namespace 9d77559
```

Then, from the namespace checkout root in MATLAB:

```matlab
addpath('benchmarks/namespaces');
report = compare_layouts('/tmp/imrv2-before-namespace',pwd);
```

Equivalent shell invocation from the namespace checkout root:

```bash
matlab -batch "addpath('benchmarks/namespaces'); compare_layouts('/tmp/imrv2-before-namespace',pwd)"
```

The paths may be changed. The baseline must contain the original flat-layout implementations, not the namespace compatibility adapters. The benchmark compares the seven numerical outputs and the input record before timing. It measures nodal finite differences, nodal Chebyshev, spectral thermal fields and nodal heat/mass transfer, using the same physical inputs and tolerances in both layouts.

Each case uses five paired `timeit` measurements with alternating execution order. MATLAB's own path is restored between layouts; only the requested source root is added, and function resolution is checked. Both sides request eight outputs and capture console output identically. The original MATLAB path and working directory are restored on exit, including errors. Run in a fresh MATLAB session with a stable machine load and repeat suspicious results.

The default gate flags a median namespace/baseline time ratio above 1.03. This is an explicit measurement-noise allowance, not proof that a smaller slowdown is absent. To use a tighter budget and more measurements:

```matlab
report = compare_layouts('/tmp/imrv2-before-namespace',pwd,9,0.01);
```

Profile any reproducible regression before accepting the change. Numerical output mismatch also fails the comparison. No speed result is claimed here: this environment has no MATLAB runtime, so the comparison remains unexecuted.
