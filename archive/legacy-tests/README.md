# Legacy test histories

This directory preserves the 2,328 opaque-name MATLAB files formerly stored directly under `tests/`. They were moved without changing their contents. Their original filenames and MATLAB variable names are retained.

The current suites do not read these files. They contain older radius histories rather than the case-keyed baseline structure expected by `tests/suites/s_test_suite_regression.m`. Do not use them as substitutes for missing approved baselines or regenerate them during routine development.

The original case-to-filename mapping and generation conditions are not documented here. Retain these files for historical investigation until their provenance can be established. The files were already present in repository history before this reorganization; no new simulation was run to create them.

Current test organization and the distinction between reference predictions and generated golden histories are explained in [the test guide](../../tests/README.md).
