# Initial stress

`f_init_stress` supplies the stress history at maximum radius for stress 3 and 4 collapse runs.
The [collapse guide](../running/collapse.md) explains when it runs and when `stress0` is required.

## Shooting algorithm

1. **Bracket**: the maximum radius grows with the initial wall velocity and equals `Req` at rest. Starting from `[0, 1]`, the upper end doubles until the maximum radius exceeds 1. The growth equation is singular at `Rdot = CL`, so the bracket stays below the sound speed; if even a nearly sonic start stays below `R = 1`, the function stops with an error that asks for `stress0`.
2. **Shooting**: `fzero` on the bracket finds the initial wall velocity for which `f_maxwell_growth_iter` reaches a maximum radius of 1.
3. **Time Integration**: `ode23tb` on `f_RP_bub_collapse_comp`, with an event at `Rdot = 0` to locate the maximum radius.
4. **Stress Extraction**: the stress state at the maximum radius.

### Internal Functions

- `f_maxwell_growth_iter`: returns the maximum radius for a trial initial wall velocity.
- `f_RP_bub_collapse_comp`: Nondimensional Keller–Miksis ODE for the growth with the stress evolution. For the Maxwell branch the stress relaxes toward the elastic stress `Se` of the (quadratic) neo-Hookean model, the same `Ze/R^3` as stress 3 and 4 in [`f_stress`](../reference/functions/f_stress.md): `Sdot = -4*Rdot/(Re*R*De) + dSe/dt - (S - Se)/De`.

## Standalone approximation

1. **Constant Definitions**: Compute `trc` (time constant) and coefficient `B`.
2. **Barotropic Correction** (`fbarbc`): Accounts for vapor and compressibility effects.
3. **Strain Correction** (`fbarst`): Originating from surface tension scaling.
4. **Compressibility Correction** (`fbarc`): Based on Grüneisen parameter.
5. **Viscous Correction** (`fbarv`): Derived from Reynolds number.
6. **Elastic Correction** (`fbare`): Analytical term for elastic effects.
7. **Stress Summation** (`fsum`): Aggregates all correction terms, with input validation.
8. **Time Constant Calculation** (`tg`): Computes relaxation time for stress growth.
9. **Stress Prediction** (`Smaxpred`): Final predicted maximum stress.

The standalone `f_approx_init_stress` is not called by either solver. It calls `hypergeom`, which is not bundled and normally requires Symbolic Math Toolbox. Confirm availability with `which hypergeom` in your MATLAB installation. Its approximate result is not a substitute for an approved collapse stress history.
