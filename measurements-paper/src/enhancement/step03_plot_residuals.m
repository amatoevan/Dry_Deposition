%% STEP03_PLOT_RESIDUALS
% Plot theoretical-minus-observed deposition velocity for four FRAGMENT
% ustar bins and all retained OCOTILLO CDP observations. Particle
% diameter is displayed on a linear horizontal axis.

clearvars
close all
clc

fprintf('Running step03_plot_residuals\n')

%% Paths and inputs
thisDir = fileparts(mfilename('fullpath'));
srcDir = fileparts(thisDir);
projectRoot = fileparts(srcDir);
addpath(srcDir)
inputFile = fullfile(projectRoot,'data','enhancement','results', ...
    'wd_difference_alpha0p06_gamma0p4.mat');
anzaInputFile = fullfile(projectRoot,'data','ocotillo','derived', ...
    'eddypro_10min_fdiff_pc10_fc', ...
    'step01b_cdp_mc_selected_10min_7bins_paired_reference_heights.mat');
figureDir = fullfile(projectRoot,'output');
baseName = 'Deltawd';
pngFile = fullfile(figureDir,[baseName '.png']);
pdfFile = fullfile(figureDir,[baseName '.pdf']);
figFile = fullfile(figureDir,[baseName '.fig']);

S = load(inputFile,'fragment','anzaAll');
F = S.fragment;
A = S.anzaAll;
anzaInput = load(anzaInputFile,'nominal','referenceHeight_m');

colors = [ ...
    0.0000 0.4470 0.7410; ...
    0.8500 0.3250 0.0980; ...
    0.9290 0.6940 0.1250; ...
    0.4660 0.6740 0.1880];
anzaColor = [0.4940 0.1840 0.5560];
symlogScale = 0.1; % cm s^-1

%% Figure
fig = figure(96);
clf(fig)
set(fig,'Color','w')
figure(fig)
ax = subplot(2,2,1);
hold(ax,'on')
box(ax,'on')
disableDefaultInteractivity(ax)
ax.Toolbar.Visible = 'off';
yline(ax,0,'k-','LineWidth',0.8,'HandleVisibility','off');

allDiameter = [F.Dp_um;A.Dp_um];
xMaximum = 10*ceil(max(allDiameter)/10);
dpSettling = linspace(0,xMaximum,400)';
dpSettling(1) = max(0.01,0.5*min(allDiameter));

% One settling-velocity reference across the complete diameter range.
met = anzaInput.nominal.meteorology;
LReference = finite_median_obukhov(met.L(:),1e10);
zReference = median(expand_height( ...
    anzaInput.referenceHeight_m.fitReference,numel(met.ustar)), ...
    'omitnan');
outSettling = deposition_velocity( ...
    median(met.ustar,'omitnan'),median(met.T,'omitnan'), ...
    median(met.p,'omitnan'),LReference,zReference, ...
    dpSettling*1e-6,'vdMethod','MCzhang01','wsMethod','bb', ...
    'z0',median(met.z0,'omitnan'),'rhop',2680,'ARp',1.5, ...
    'alpha',0.06,'gamma',0.4);
hSettling = plot(ax,dpSettling, ...
    signed_log(-100*outSettling.ws(:),symlogScale),'--k', ...
    'LineWidth',1,'DisplayName','$-w_s$');

windowLow = unique(F.ustarBinLow,'stable');
hFragment = gobjects(numel(windowLow),1);
for iw = 1:numel(windowLow)
    use = F.ustarBinLow==windowLow(iw);
    B = sortrows(F(use,:),'Dp_um');
    hFragment(iw) = plot(ax,B.Dp_um, ...
        signed_log(B.deltaWd_cm_s,symlogScale),'-o', ...
        'Color',colors(iw,:),'MarkerFaceColor','w', ...
        'MarkerSize',4.5,'LineWidth',1, ...
        'DisplayName',sprintf('$u_* = %.2f$--$%.2f$', ...
        B.ustarBinLow(1),B.ustarBinHigh(1)));
end

A = sortrows(A,'Dp_um');
hAnza = plot(ax,A.Dp_um, ...
    signed_log(A.deltaWd_cm_s,symlogScale),'-s', ...
    'Color',anzaColor,'MarkerFaceColor','w', ...
    'MarkerSize',5.5,'LineWidth',1, ...
    'DisplayName',sprintf('OCOTILLO: median $u_*=%.2f$',A.ustarMedian(1)));

set(ax,'XScale','linear','FontSize',12,'TickLabelInterpreter','latex')
xlim(ax,[0 xMaximum])
ax.XTick = 0:10:xMaximum;
yTickValues = [-100 -30 -10 -3 -1 -0.3 0 0.3 1];
ax.YTick = signed_log(yTickValues,symlogScale);
ax.YTickLabel = compose('%g',yTickValues);
ylim(ax,signed_log([-100 1],symlogScale))
xlabel(ax,'$D_p$ ($\mu$m)','Interpreter','latex','FontSize',12)
ylabel(ax,'$\Delta w_d$ (cm s$^{-1}$)', ...
    'Interpreter','latex','FontSize',12)
title(ax,'$\alpha=0.06$, $\gamma=0.4$', ...
    'Interpreter','latex','FontWeight','normal','FontSize',12)

panelPosition = ax.Position;
legend(ax,[hFragment;hAnza;hSettling], ...
    'Location','eastoutside','Interpreter','latex','FontSize',8);
ax.Position = panelPosition;

if ~isfolder(figureDir)
    mkdir(figureDir)
end
exportgraphics(fig,pngFile,'Resolution',300)
exportgraphics(fig,fullfile(figureDir,[baseName '.jpg']),'Resolution',300)
exportgraphics(fig,pdfFile,'ContentType','vector')
savefig(fig,figFile)

fprintf('Saved four-bin linear-Dp figures:\n  %s\n  %s\n  %s\n', ...
    pngFile,pdfFile,figFile)


function transformed = signed_log(value,scale)
transformed = sign(value).*log10(1+abs(value)./scale);
end


function value = finite_median_obukhov(L,threshold)
L(isinf(L) & L>0) = threshold;
L(isinf(L) & L<0) = -threshold;
value = median(L,'omitnan');
end


function z = expand_height(z,nTime)
if isscalar(z)
    z = repmat(z,nTime,1);
else
    z = z(:);
end
if numel(z)~=nTime
    error('Reference-height length does not match the time dimension.')
end
end
