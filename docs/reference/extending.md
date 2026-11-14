# Extending IMR

You extend IMR by editing its MATLAB functions directly. There is no plugin loader. Each input must pass from defaults through parsing and validation to the equations that use it.

For a new fluid that uses an existing EoS family, start with a [MATLAB fluid profile](../model/fluids.md). Material coefficients and saturation data belong there. The shared evaluator and SI converter handle selection. The trace below applies when you add a numerical option or a mathematical model that the code does not already support.

## Worked trace: the existing `scheme` input

| Step | Existing implementation | What an extension must preserve |
|---|---|---|
| Default | `src/+imr/+configuration/defaults.m`: `d.scheme = 0` | One authoritative default |
| Public struct | `imr.default_inputs` uses `imr.configuration.schema` | Schema public name and internal field match |
| Parse | `imr.configuration.prepare_case` uses that schema | Override and struct calls prepare the same case |
| Validate | Schema discrete rule and `imr.configuration.validate_case` | Invalid values and incompatible combinations fail |
| Prepare | Named group on the case | Both factories read the intended field |
| Use | Private nodal/spectral factory | The setting changes the intended operator or closure |
| Verify | Focused MATLAB suite | Independent numerical check covers the new behavior |

Follow that trace for a new input. Add the parameter to `imr.configuration.defaults` and `imr.configuration.schema`, then to its named group in `imr.configuration.prepare_case`. If the input affects a physical model, update its state layout, right-hand side, Jacobian and output reconstruction together. A parser change alone is not enough.

A new stress or viscosity law also needs checker ranges, labels in `imr.solvers.describe_case`, compatibility decisions and tests of limiting behavior. Keep solver paths compatible with base MATLAB. Discuss intended numerical-result changes before you update stored references. See [contributing](contributing.md).

This worked example explains existing code. It does not add an input or modify solver behavior.
