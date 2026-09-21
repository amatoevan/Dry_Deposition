# Measurements-paper figure archive

Private working archive for the dry-deposition measurements manuscript,
assembled on 2026-09-21. This folder is self-contained: copy it out of the
parent repository and the figure runner still works. It intentionally
contains only figure code, its selected inputs, and figure artwork—not the
full field-processing pipeline or manuscript drafts.

## Reproduce the figures

Use MATLAB R2023b with Statistics and Machine Learning Toolbox (`prctile`).
Start MATLAB in this folder and run:

```matlab
run_measurements_figures
```

The eight manuscript JPEGs are written to `output/`, using the exact filenames
referenced by the current Overleaf manuscript. Computational figures also
produce editable `.fig` and PDF files. The original manuscript exports are
preserved separately in `reference_figures/`; regeneration never overwrites them.
No Monte Carlo fitting is performed. FRAGMENT and OCOTILLO envelopes use up to
500 saved paired parameter draws with reproducible random seed 2.

Individual figures can be regenerated, for example:

```matlab
run_measurements_figures('fragment')
run_measurements_figures('anza')
run_measurements_figures('residuals')
run_measurements_figures('map')
run_measurements_figures('combined')
run_measurements_figures('artwork')
```

## Figure inventory

| Manuscript figure / label | Overleaf filename | Reproduction source |
| --- | --- | --- |
| Measurement schematic | `SchematicMeasure.jpg` | Author-created artwork; exact image copied. Related editable PowerPoint in `assets/schematic/`. |
| Field location map | `figLocation.jpg` | `src/location/figLocation.m`, using the original NOAA-21 raster stored losslessly in `assets/location/NOAA-21.mat`. |
| FRAGMENT deposition | `Fragment_wd_default_fit.jpg` | `src/fragment/sensitivity/plot_fragment_zhang_parameter_ensemble.m`; paired Zhang and FRAGMENT fits. |
| OCOTILLO deposition | `Anza_wd_default.jpg` | `src/ocotillo/plot_default_vs_fragment_lambda0.m`, `ManuscriptDefault=true`. |
| Deposition residuals | `Deltawd.jpg` | `src/enhancement/step03_plot_residuals.m`; saved residual tables. |
| Constant / bin-specific enhancement | `combined.jpg` | Final manually assembled MATLAB layout in `assets/enhancement/combined.fig`; numerical summaries and original panel-generating scripts are also included. |
| Single-particle microscopy | `SingleParticle.jpg` | Author-created artwork; exact image copied. Editable PowerPoint in `assets/particles/`. |
| Elemental composition | `elemental.jpg` | Exact author-supplied export copied. PDF and accompanying `density.xlsx` workbook in `assets/particles/`; no source plotting script was found. |

The artwork files are preserved as supplied. The available schematic PowerPoint
is supporting editable source; the exact manuscript JPEG remains authoritative.
The elemental workbook is accompanying data, not a claim that it independently
reproduces the complete elemental graphic.

## Data and uncertainty

- `data/fragment/`: saved meteorology, accepted measurement Monte Carlo
  realizations, plotting summaries, and paired alpha/gamma fits.
- `data/ocotillo/`: the selected 26-period CDP dataset, cached observation
  summaries, and saved fitted parameters from `eddypro_10min_fdiff_pc10_fc`.
- `data/zhang2014/`: the same 2,000-realization laboratory-fit result used in
  the submitted theory paper. Accepted alpha/gamma medians are approximately
  3.291 and 0.34897. This archive reuses that fit rather than refitting it.
- `data/enhancement/`: residual tables, saved lambda ensembles, and numerical
  summaries for the two enhancement panels.

In the updated FRAGMENT and OCOTILLO default-parameter figures, shading is
pointwise pooled 25th–75th percentiles across meteorological states and paired
alpha/gamma draws. Lines are pooled medians. Titles show full-ensemble
parameter medians and 2.5th–97.5th percentile intervals. Alpha and gamma are
never sampled independently. OCOTILLO uses 26 states; FRAGMENT uses four
friction-velocity bins with 103, 218, 255, and 86 states.

The enhancement figure is preserved as previously prepared. Its panel scripts
specify 2.5th–97.5th percentile shading: meteorological variation with fixed
lambda on the left, meteorological and lambda variation on the right. These
are not the newly added alpha/gamma IQR envelopes. This packaging task does
not revise that figure's scientific calculation or the manuscript captions.

## Recompute enhancement panels (optional)

The default runner exports the final editable two-panel figure, including its
stored plotted data. To rebuild the constituent panels from the included
inputs, run these original scripts in order:

```matlab
run('src/enhancement/step06_plot_wd_bin_specific_lambda.m')
run('src/enhancement/step07_plot_wd_constant_lambda_four_bins.m')
```

These scripts recompute theory envelopes, can be substantially slower, write
individual panels to `figures/enhancement/`, and do not recreate the manual
combined-panel arrangement. No refitting is needed. Their cached output
statistics are already supplied. The physical-model functions are the actual
workspace versions used when this archive was assembled; older cached figures
remain archived rather than being silently replaced by new calculations.

## Provenance and scope

`provenance.json` maps included files to their source paths in the working
project. `SHA256SUMS` records the packaged inputs, source files, and reference
artwork. Paths stored inside legacy MAT metadata describe the original run;
executable scripts use paths relative to this folder.

The existing scripts are retained. Packaging changes are limited to portable
input/output paths, removal of unused map terrain/colormap dependencies, a
losslessly compressed map input, and the small figure runner. No fitting
algorithms were rewritten and no simulation was rerun during packaging.

The large selected FRAGMENT data file and schematic PowerPoint are included
because they support the figure sources; all files are below GitHub's 100 MiB
per-file limit. Raw sensor archives, unrelated analyses, review comments,
manuscript text, temporary outputs, and credentials are excluded.

Keep this archive private until the authors choose to release it. No new
license or public-release commitment is made by this archive.
