# Water solver distribution

This distribution supports water as the ambient fluid with all six radial models (Rayleigh–Plesset, Keller–Miksis, Gilmore variants), transport properties, viscoelastic stress models, heat and mass transfer, and both nodal and spectral discretization schemes. The water profile uses a legacy reference set with `rho8 = 1064 kg/m³` and constant latent heat / caloric offset model for vapor; it is not a temperature-dependent water equation of state.

All source functions, solvers, Jacobians and tests are generic and do not depend on a specific fluid. New fluid profiles can be added by copying `src/+imr/+fluids/water.m`, editing the material data and registering the fluid with the configuration system. See [Fluid profiles](../model/fluids.md) for the schema and conventions.
