# Calibration data

`zhang2014/raw/Zhang14.csv` contains the measurements transcribed for this study from Zhang, J., Shao, Y., and Huang, N. (2014), “Measurements of dust deposition velocity in a wind-tunnel experiment,” *Atmospheric Chemistry and Physics* 14, 8869-8882, [doi:10.5194/acp-14-8869-2014](https://doi.org/10.5194/acp-14-8869-2014).

The table includes wood and water surfaces as well as the three bare surfaces. The fitting code selects the bare surfaces (sand, sandy loam, and sand and gravel). Diameters are in micrometres; deposition velocities and their reported uncertainty are in millimetres per second. The header rows record friction velocity and height or roughness in metres. `src/zhang2014/step01_ingest_zhang2014.m` converts these to SI units and writes `data/zhang2014/derived/Zhang14.mat`.

The 2,000-realization fitting step writes two result files to `data/zhang2014/results/`. These derived files are generated on demand and ignored by Git. The included CSV is the only input data needed by the released workflow.
