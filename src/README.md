# Source code

`run_paper_figures.m` is the entry point. All figure scripts live directly in `src/`; the `zhang2014/` subdirectory holds the data ingest and fitting workflow.

| File | Role |
| --- | --- |
| `deposition_velocity.m` | Mass-conserving and electric-analogue deposition calculations used in the paper |
| `settling_velocity.m` | Settling schemes compared in Figures 4, 5, and 9 |
| `zhang2014/step01_ingest_zhang2014.m` | Convert the included CSV to a standardized MATLAB table |
| `zhang2014/step02_fit_deposition_default_n2000.m` | Run the paired 2,000-realization fits for Figures 3, 6-9 |
| `zhang2014/step02_fit_deposition_default.m` | Fit one deposition formulation |
| `zhang2014/utils/fit_zhang2014_alpha_gamma_mc.m` | Parameter optimization and uncertainty sampling |
| `zhang2014/step03_plot_deposition_default.m` | Plot the calibrated curves for Figures 3 and 6 |
| `plot_rebound_probability.m` | Figure 2 |
| `plot_physics_comparisons.m` | Figures 4 and 5, plus comparison plots used as bases for Figures 7-9 |
| `plot_parameter_uncertainty.m` | Uncertainty ribbons for Figures 5 and 7-9 |

The paired result files are `zhang14_fit_default_cliftgauvin_n2000.mat` and `zhang14_fit_default_cliftgauvin_zhang01_n2000.mat` under `data/zhang2014/results/`. They are reproducible outputs and are not committed. `run_paper_figures(false)` needs them to be present.

