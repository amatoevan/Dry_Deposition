function run_paper_figures(refit)
%RUN_PAPER_FIGURES Rebuild the analysis figures for the submitted paper.
%
%   run_paper_figures() regenerates the Zhang et al. fits and then
%   all computational figures. run_paper_figures(false) reuses the
%   two existing n=2000 result files.

if nargin<1
    refit = true;
end

thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(thisDir);
srcDir = thisDir;
figureDir = fullfile(projectRoot,'figures','generated');
if ~exist(figureDir,'dir'), mkdir(figureDir); end
zhangDir = fullfile(srcDir,'zhang2014');
addpath(srcDir,zhangDir,fullfile(zhangDir,'utils'))

step01_ingest_zhang2014();
if refit
    step02_fit_deposition_default_n2000();
end

step03_plot_deposition_default( ...
    'fitTag','default_cliftgauvin_n2000', ...
    'outputBase','ZhangFit','figDir',figureDir);
step03_plot_deposition_default( ...
    'fitTag','default_cliftgauvin_zhang01_n2000', ...
    'outputBase','ZhangFit_wde','figDir',figureDir);

run(fullfile(thisDir,'plot_rebound_probability.m'))
run(fullfile(thisDir,'plot_physics_comparisons.m'))
run(fullfile(thisDir,'plot_parameter_uncertainty.m'))

thisDir = fileparts(mfilename('fullpath'));
fprintf('Generated figures written to %s\n',figureDir)
end
