# Measured pressure waveforms

The three measured driving-pressure waveforms (`wave_type = -1, -2, -3`) are from:

E. Hersey, M. Rodriguez and E. Johnsen, "Dynamics of an oscillating microbubble in a blood-like Carreau fluid," *J. Acoust. Soc. Am.* **153**(3), 1836–1845 (2023). <https://doi.org/10.1121/10.0017342>

Please cite this paper if you use these waveforms.

## Files

| File | Contents |
|---|---|
| `patterson_waveform_data.mat` | The three waveforms, stored as dimensionless time and pressure columns: mostly linear (ML), moderately nonlinear (MN) and highly nonlinear (HN) |
| `r_patterson_data.m` | Fits a cubic spline to each waveform and its derivative, plots them, and saves them as `ml.mat`, `mn.mat` and `hn.mat`. It uses `fnder` from the Curve Fitting Toolbox. |

The solvers load the spline fits from `src/`, through `f_call_params`:

| `wave_type` | File | Waveform |
|---|---|---|
| -1 | `src/d_hn.mat` | Highly nonlinear |
| -2 | `src/d_mn.mat` | Moderately nonlinear |
| -3 | `src/d_ml.mat` | Mostly linear |

Each file holds the spline `poly` of the dimensionless pressure and its derivative `dpoly`, over dimensionless times up to 15. [`f_pinfinity`](../src/f_pinfinity.md) evaluates them.
