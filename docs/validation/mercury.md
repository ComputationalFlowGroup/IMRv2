# Mercury profile validation

The mercury MATLAB suite checks the optional material profile, properties and KM preparation. The separate [polytropic suite](polytropic.md) checks the shared constant-vapor derivative. They need base MATLAB and bundled `.m` reference values, not Python, another checkout or stored regression `.mat` files.

From the checkout root in a shell with licensed MATLAB:

```bash
cd toolchain
matlab -batch "s_test_suite_mercury"
```

Or change to `toolchain/` in MATLAB and run the script there. It runs as a separate job in the existing MATLAB CI workflow.

## Claims checked

- MATLAB properties reproduce 24 frozen source-reference states covering pressure and temperature, with relative tolerance `5e-12`.
- All 115 sound-speed observations retain the original canonical model predictions and its 0.712489418% RMS residual.
- Scalar and array calls agree; zero-pressure limits, temperature-specific moduli and invalid input handling behave as documented.
- Saturation interpolation matches measured nodes, is monotone and interpolates in reciprocal temperature without extrapolation.
- Input preparation uses Hg density, acoustics, gas constant and saturation pressure, derives consistent scales, records replay metadata and rejects missing transport values or conflicting overrides.
- The nodal FD, nodal Chebyshev and spectral KM calls match equivalent explicitly supplied SI properties and preserve an equilibrium bubble.

`tests/fixtures/+imrtest/mercury_reference.m` stores source calculations and observations as MATLAB literals with the source commit and original data citation. They were prepared from the original evaluator independently of the MATLAB port, not regenerated from the implementation under test. The reference-state property values were checked while preparing this change; they are not MATLAB execution results.

## Execution status and scientific limits

MATLAB and Octave were not installed or run in the cloud task. These numerical suites, the example, MATLAB lint and the repository formatter remain for local execution. Static MATLAB syntax checks do not establish a solver pass.

The source sound-speed RMS is agreement with fit observations, not independent mercury-bubble validation. Saturation acoustics have the source's larger uncertainty. Dynamic mass transfer, variable liquid properties during a run and enthalpy/Gilmore coupling are outside this implementation. The separate [regression suite](regression.md#baseline-availability) still lacks its approved baseline files.
