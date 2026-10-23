# Defaults, material data and units

`imr.configuration.defaults` is still used. Both `imr.configuration.prepare_case` and `imr.default_inputs` call it; removing it would break the default solver API. Its physical groups are now derived through `imr.nondimensionalize` instead of repeating a second list of physical numbers.

## One source per kind of value

| Kind of value | Source | Consumer |
|---|---|---|
| Liquid/vapor coefficients and reference transport | `src/+imr/+fluids/water.m` or an explicitly selected material profile | `imr.fluid_properties` and SI preparation |
| Reference geometry, ambient conditions, gas conductivity and dimensional stress/forcing choices | `src/+imr/nondimensionalize.m` | Conversion to dimensionless groups |
| Numerical/model choices, grid sizes, tolerances and default output times | `src/+imr/+configuration/defaults.m` | Parser and public input struct |
| Physical dimensionless groups | Calculated by `imr.nondimensionalize`, used by `imr.configuration.defaults` | Both solvers |
| Frozen old values used to detect regressions | `tests/+imrtest/default_physics.m` | Tests only |

The snapshot in `tests/+imrtest` is intentionally independent of the current implementation. Replacing it with a call to the converter would make a default-change regression test unable to detect a mistake.

## Why the old `p.defaults` group was confusing

That group mixed SI transport values and conversion defaults, including artificial scales for inactive terms. It was not a set of enthalpy coefficients. It has been removed. Profiles now name reference transport fields individually under `p.transport`, with a unit and meaning beside each assignment. Missing data remain absent; an active model must obtain required inputs from the caller.

The converter still needs finite scales for some inactive terms because the solver input interface contains thermal and transport groups even in polytropic runs. It supplies neutral numerical scales internally when material data are absent. They are documented implementation scales, not measured transport properties. They never satisfy a required active-model input.

Material enthalpy coefficients live under `liquid` and `vapor`. A caloric slice such as $h_l=C_{p,l}T+q_l$ is distinct from the pressure integral used by an enthalpy radial equation. See [radial theory](../theory/radial-equations.md#enthalpy-forms).

## Units and provenance

The unused `p.units` string has been removed. SI is the fixed contract of the property and conversion APIs, with units documented beside material fields and in the [API tables](../api/f_fluid_properties.md). The solver interface remains dimensionless; reference scales restore SI outputs.

`name`, `source` and `commit` remain useful provenance metadata: they are preserved in `ref.fluid.profile` when a material is selected. They do not affect numerical equations and are not a unit-conversion system.

## Three ways to keep the code streamlined

1. Derive physical groups through the shared converter instead of maintaining parallel SI and dimensionless constants. This is implemented.
2. Keep one public path setup and use namespaced data profiles; do not add `genpath` calls throughout examples. This is implemented.
3. Keep optional material profiles, documentation, tests and fixtures separate from generic physics. The paper branch can omit an add-on without editing the solvers or replacing generic tests.

Changes in the computed defaults are limited to floating-point rounding of the previous literals. The fluid suite compares against frozen pre-consolidation values with a relative/absolute tolerance of `5e-13`. Numerical MATLAB execution remains for local verification.
