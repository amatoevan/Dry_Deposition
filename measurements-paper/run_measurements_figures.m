function run_measurements_figures(selection)
%RUN_MEASUREMENTS_FIGURES Rebuild manuscript figures without refitting.
% run_measurements_figures() rebuilds all computational figures and copies
% the author-created schematic/microscopy artwork into output/.
% Select 'fragment', 'anza', 'residuals', 'map', 'combined', or 'artwork'
% to rebuild only that item. Original manuscript images: reference_figures/.
if nargin<1, selection = 'all'; end
selection = validatestring(selection, ...
    {'all','fragment','anza','residuals','map','combined','artwork'});
root = fileparts(mfilename('fullpath'));
outDir = fullfile(root,'output');
if ~isfolder(outDir), mkdir(outDir), end
oldPath = path;
restorePath = onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(fullfile(root,'src'),fullfile(root,'src','ocotillo'))
if any(strcmp(selection,{'all','artwork'}))
    for name = {'SchematicMeasure.jpg','SingleParticle.jpg','elemental.jpg'}
        copyfile(fullfile(root,'reference_figures',name{1}),outDir);
    end
end
if any(strcmp(selection,{'all','fragment'}))
    run_isolated(fullfile(root,'src','fragment','sensitivity', ...
        'plot_fragment_zhang_parameter_ensemble.m'));
end
if any(strcmp(selection,{'all','anza'}))
    plot_default_vs_fragment_lambda0(ManuscriptDefault=true);
end
if any(strcmp(selection,{'all','residuals'}))
    run_isolated(fullfile(root,'src','enhancement','step03_plot_residuals.m'));
end
if any(strcmp(selection,{'all','map'}))
    run_isolated(fullfile(root,'src','location','figLocation.m'));
end
if any(strcmp(selection,{'all','combined'}))
    % This final two-panel layout was assembled manually in MATLAB.
    % Its editable FIG contains the plotted X/Y data and uncertainty bands.
    % Preserve that layout rather than repeat the costly lambda calculations.
    fig = openfig(fullfile(root,'assets','enhancement','combined.fig'),'invisible');
    savefig(fig,fullfile(outDir,'combined.fig'));
    exportgraphics(fig,fullfile(outDir,'combined.jpg'),'Resolution',300);
    exportgraphics(fig,fullfile(outDir,'combined.pdf'),'ContentType','vector');
    close(fig)
end
fprintf('Manuscript figure outputs: %s\n',outDir)
end

function run_isolated(script)
% Existing scripts use clearvars; keep each in its own workspace.
run(script)
end
