# Validation — 2026-09-21

Ran `run_measurements_figures` in MATLAB R2023b from an independent checkout
in `/tmp`, with no access through package-relative paths to the original
working project's data. The run completed with exit code 0 and produced
all eight expected manuscript JPEGs. No fits were rerun.

- FRAGMENT: 500 paired draws in each panel; bin counts 103, 218, 255, 86.
  Existing vectorized-model check: maximum relative differences
  2.44e-14 (Zhang) and 1.39e-14 (FRAGMENT).
- OCOTILLO: 26 meteorological states x 500 paired Zhang draws.
- Residual plot: regenerated from the supplied saved tables.
- Map: regenerated from the packaged raster without the external terrain
  dataset or personal MATLAB colormap file.
- Enhancement comparison: exported from the saved editable two-panel FIG.
- Schematic and particle/elemental graphics: copied from reference artwork.

Visually inspected the regenerated residual, map, and combined outputs.
FRAGMENT and OCOTILLO had already been visually checked in this session.
The optional expensive enhancement-panel recalculations were not run.
The separate original manuscript images remain the reference for manual
layout changes and export-resolution differences.
