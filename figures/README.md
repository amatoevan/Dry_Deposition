# Submitted paper figures

`paper/` contains the nine JPG files supplied with the ACP submission. They are retained as the submitted assets; reproduction runs write to `generated/` so they do not overwrite those files.

| Figure | Submitted asset | Source |
| --- | --- | --- |
| 1 | `paper/Figure_01_Schematic.jpg` | Author-created schematic; no MATLAB generator |
| 2 | `paper/Figure_02_R1.jpg` | `src/plot_rebound_probability.m` |
| 3 | `paper/Figure_03_ZhangFit.jpg` | `src/zhang2014/step03_plot_deposition_default.m`, mass-conserving fit |
| 4 | `paper/Figure_04_ws.jpg` | `src/plot_physics_comparisons.m` |
| 5 | `paper/Figure_05_wdws.jpg` | `src/plot_physics_comparisons.m`; uncertainty version from `src/plot_parameter_uncertainty.m` |
| 6 | `paper/Figure_06_ZhangFit_wde.jpg` | `src/zhang2014/step03_plot_deposition_default.m`, electric-analogue fit |
| 7 | `paper/Figure_07_wdewd0_parameter_IQR.jpg` | `src/plot_parameter_uncertainty.m` |
| 8 | `paper/Figure_08_wdewdz_parameter_IQR.jpg` | `src/plot_parameter_uncertainty.m` |
| 9 | `paper/Figure_09_TotalDiff_parameter_IQR.jpg` | `src/plot_parameter_uncertainty.m` |

Generated `.fig`, image, and PDF versions are local outputs. The uncertainty script also produces a version of Figure 5 with parameter ribbons (`wdws_parameter_IQR`); the submitted Figure 5 image in `paper/` is the `wdws` file.
