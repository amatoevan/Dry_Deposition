%% STEP06_PLOT_WD_BIN_SPECIFIC_LAMBDA
% Plot four FRAGMENT ustar bins and the selected OCOTILLO CDP data.
% Each group uses its own fitted lambda distribution. The theoretical bands
% contain both retained-time meteorological variation and lambda uncertainty.

clearvars
close all
clc

fprintf('Running step06_plot_wd_bin_specific_lambda\n')

%% Settings and paths
saveFig = true;
fragmentUstarEdges = [0 0.05 0.10 0.15 0.20];
DpCurve_um = logspace(log10(0.1),log10(50),200);
xLimits = [0.3 50];
xTicks = [0.3 1 3 10 20 50];
yLimits = [3e-3 1e2];

thisDir = fileparts(mfilename('fullpath'));
srcDir = fileparts(thisDir);
projectRoot = fileparts(srcDir);
addpath(srcDir)

fragmentFile = fullfile(projectRoot,'data','fragment','derived', ...
    'Fragment_fit_ready_MC.mat');
lambdaFile = fullfile(projectRoot,'data','enhancement','results', ...
    'enhancement_lambda_four_fragment_bins_ocotillo_mc.mat');
anzaSelectedFile = fullfile(projectRoot,'data','ocotillo','derived', ...
    'eddypro_10min_fdiff_pc10_fc', ...
    'step01b_cdp_mc_selected_10min_7bins_paired_reference_heights.mat');
anzaStatsFile = fullfile(projectRoot,'data','ocotillo','results', ...
    'eddypro_10min_fdiff_pc10_fc', ...
    'ocotillo_wd_10min_7bins_default_fragment_lambda_at_2p1m_statistics.mat');
anzaFitFile = fullfile(projectRoot,'data','ocotillo','results', ...
    'eddypro_10min_fdiff_pc10_fc', ...
    'ocotillo_CDP_fit_10min_7bins_fourCase_alpha_gamma_lambda_MC.mat');

figDir = fullfile(projectRoot,'figures','enhancement');
dataFile = fullfile(projectRoot,'data','enhancement','results', ...
    'enhance_wd_bin_specific_lambda_statistics.mat');
baseName = 'Figure_wd_bin_specific_lambda_four_fragment_bins_cdp';
pngFile = fullfile(figDir,[baseName '.png']);
pdfFile = fullfile(figDir,[baseName '.pdf']);
figFile = fullfile(figDir,[baseName '.fig']);

inputFiles = {fragmentFile;lambdaFile;anzaSelectedFile; ...
    anzaStatsFile;anzaFitFile};
for i = 1:numel(inputFiles)
    if ~isfile(inputFiles{i})
        error('Required input not found:\n%s',inputFiles{i})
    end
end

F = load(fragmentFile);
E = load(lambdaFile,'R');
D = load(anzaSelectedFile);
A = load(anzaStatsFile,'R');
Q = load(anzaFitFile,'fits');

alphaFragment = A.R.settings.alphaFragment;
gammaFragment = A.R.settings.gammaFragment;
fragmentLambdaMC = E.R.fragment.lambdaMonteCarlo;
cdpFit = Q.fits.fragmentAlphaGamma_fitLambda;
cdpLambdaSamples = cdpFit.monteCarlo.lambda(cdpFit.monteCarlo.valid);

validate_inputs(F,E.R,D,A.R,cdpFit,fragmentUstarEdges)

%% Build observed summaries and lambda-aware theoretical envelopes
fragment = summarize_fragment(F,fragmentUstarEdges,DpCurve_um, ...
    alphaFragment,gammaFragment,fragmentLambdaMC);
cdp = summarize_cdp(D,A.R,DpCurve_um,alphaFragment,gammaFragment, ...
    cdpLambdaSamples);

%% Plot
colors = [ ...
    0.0000 0.4470 0.7410
    0.8500 0.3250 0.0980
    0.9290 0.6940 0.1250
    0.4660 0.6740 0.1880
    0.4940 0.1840 0.5560];

fig = figure(90);
clf(fig)
set(fig,'Color','w')
ax = subplot(2,2,1,'Parent',fig);
ax.Position(2) = 0.25;
ax.Position(4) = 0.67;
hold(ax,'on')
box(ax,'on')

h = gobjects(6,1);
labels = strings(6,1);
for ib = 1:4
    plot_envelope(ax,DpCurve_um,fragment.theory,ib,colors(ib,:))
    plot_observations(ax,F.diam_mean_um, ...
        fragment.obs.med(ib,:),fragment.obs.q25(ib,:), ...
        fragment.obs.q75(ib,:),'s',colors(ib,:))
    labels(ib) = sprintf( ...
        ['FRAGMENT: $u_*=(%.2f$--$%.2f)$, ' ...
         '$\\lambda=%.2f$ [%.2f, %.2f], $N=%d$'], ...
        fragmentUstarEdges(ib),fragmentUstarEdges(ib+1), ...
        fragment.lambda.median(ib),fragment.lambda.ci95(ib,1), ...
        fragment.lambda.ci95(ib,2),fragment.Ntime(ib));
    h(ib) = plot(ax,NaN,NaN,'s-','Color',colors(ib,:), ...
        'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1.6, ...
        'DisplayName',labels(ib));
end

plot_envelope(ax,DpCurve_um,cdp.theory,1,colors(5,:))
plot_observations(ax,cdp.Dp_um,cdp.obs.med,cdp.obs.q25, ...
    cdp.obs.q75,'o',colors(5,:))
labels(5) = sprintf( ...
    'OCOTILLO CDP: $\\lambda=%.2f$ [%.2f, %.2f], $N=%d$', ...
    cdp.lambda.median,cdp.lambda.ci95(1),cdp.lambda.ci95(2), ...
    cdp.Ntime);
h(5) = plot(ax,NaN,NaN,'o-','Color',colors(5,:), ...
    'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1.6, ...
    'DisplayName',labels(5));

labels(6) = '$w_s$';
h(6) = plot(ax,DpCurve_um,100*cdp.wsReference,'--k', ...
    'LineWidth',1.5,'DisplayName',labels(6));

set(ax,'XScale','log','YScale','log', ...
    'XLim',xLimits,'XTick',xTicks,'YLim',yLimits, ...
    'FontSize',12,'TickLabelInterpreter','latex')
xlabel(ax,'$D_p$ ($\mu$m)','Interpreter','latex')
ylabel(ax,'$w_d$ (cm s$^{-1}$)','Interpreter','latex')
title(ax,sprintf( ...
    '$\\alpha=%.2f$, $\\gamma=%.1f$; bin-specific $\\lambda$', ...
    alphaFragment,gammaFragment), ...
    'Interpreter','latex','FontWeight','normal','FontSize',12)
legend(ax,h,labels,'Location','northwest', ...
    'Interpreter','latex','FontSize',7.2)
text(ax,0.98,0.025, ...
    'Shading: joint 95\% range from meteorology and $\lambda$ MC', ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','bottom','Interpreter','latex','FontSize',7.2)
grid(ax,'on')

%% Save figure and plotting statistics
S = struct;
S.workflow = 'enhancement/bin-specific-lambda';
S.created = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss Z'));
S.inputs = struct('fragmentData',fragmentFile, ...
    'fragmentLambdaFits',lambdaFile,'anzaSelectedData',anzaSelectedFile, ...
    'anzaPlotStatistics',anzaStatsFile,'anzaLambdaFits',anzaFitFile);
S.parameters = struct('alpha',alphaFragment,'gamma',gammaFragment, ...
    'DpCurve_um',DpCurve_um,'fragmentUstarEdges',fragmentUstarEdges);
S.uncertainty = ['central line: median retained-time prediction at the ' ...
    'median fitted lambda; shading: 2.5th-97.5th percentiles of all ' ...
    'predictions across retained meteorological times and all finite ' ...
    'lambda Monte Carlo realizations'];
S.fragment = fragment;
S.cdp = cdp;
S.figureFiles = struct('png',pngFile,'pdf',pdfFile,'fig',figFile);

if saveFig
    if ~isfolder(figDir)
        mkdir(figDir)
    end
    exportgraphics(ax,pngFile,'Resolution',300)
    exportgraphics(ax,pdfFile,'ContentType','vector')
    savefig(fig,figFile)
    save(dataFile,'S','-v7.3')
    fprintf('Saved bin-specific-lambda figure:\n  %s\n  %s\n  %s\n', ...
        pngFile,pdfFile,figFile)
    fprintf('Saved plotting statistics:\n  %s\n',dataFile)
end

fprintf('FRAGMENT retained-time N: %s\n',mat2str(fragment.Ntime(:)'))
fprintf('FRAGMENT lambda medians: %s\n', ...
    mat2str(fragment.lambda.median(:)',6))
fprintf('CDP retained-time N: %d\n',cdp.Ntime)
fprintf('CDP lambda median [95%% CI]: %.6g [%.6g, %.6g]\n', ...
    cdp.lambda.median,cdp.lambda.ci95)


%% Local functions
function S = summarize_fragment(F,edges,DpCurve_um,alpha,gamma,M)
u = F.meteorology.ustar(:);
T = F.meteorology.T_K(:);
p = F.meteorology.p_Pa(:);
L = F.meteorology.L_m(:);
z0 = F.meteorology.z0_m(:);
zRef = F.configuration.zReference_m;
nBin = numel(F.diam_mean_um);
[nTime,~,nMC] = size(F.mc.wd_obs);

timeEligible = isnan(F.q_15avg(:)) & ...
    isfinite(u) & u>0 & u<edges(end);
nominalMask = repmat(timeEligible,1,nBin) & ...
    isfinite(F.nominal.Fdiff) & F.nominal.Fdiff<0 & ...
    isfinite(F.nominal.wd_obs) & F.nominal.wd_obs>0;
timeRetained = any(nominalMask,2);

wdObs = F.mc.wd_obs;
for k = 1:nMC
    W = wdObs(:,:,k);
    use = nominalMask & isfinite(F.mc.Fdiff(:,:,k)) & ...
        F.mc.Fdiff(:,:,k)<=0 & isfinite(W) & W>0;
    W(~use) = NaN;
    wdObs(:,:,k) = W;
end

binID = discretize(u,edges);
nGroup = numel(edges)-1;
S.Ntime = zeros(nGroup,1);
S.obs.med = NaN(nGroup,nBin);
S.obs.q25 = NaN(nGroup,nBin);
S.obs.q75 = NaN(nGroup,nBin);
S.theory.med = NaN(nGroup,numel(DpCurve_um));
S.theory.q025 = NaN(nGroup,numel(DpCurve_um));
S.theory.q975 = NaN(nGroup,numel(DpCurve_um));
S.theory.nMeteorologicalTimes = zeros(nGroup,1);
S.theory.nJointPredictions = zeros(nGroup,1);
S.lambda.samples = M.lambda;
S.lambda.median = M.median(:);
S.lambda.ci95 = M.ci95;
S.lambda.nSuccessful = M.nSuccessful(:);
S.lambda.nRequested = M.nRequested;

validMet = timeRetained & isfinite(T) & T>0 & ...
    isfinite(p) & p>0 & ~isnan(L) & isfinite(z0) & z0>0;

for ib = 1:nGroup
    inBin = binID==ib & timeRetained;
    S.Ntime(ib) = nnz(inBin);
    [S.obs.med(ib,:),S.obs.q25(ib,:),S.obs.q75(ib,:)] = ...
        summarize_observations(wdObs,inBin);

    idx = find(binID==ib & validMet);
    lambdaSamples = M.lambda(:,ib);
    lambdaSamples = lambdaSamples(isfinite(lambdaSamples));
    E = theory_envelope(u,T,p,L,z0,repmat(zRef,nTime,1),idx, ...
        DpCurve_um,alpha,gamma,2500,1.46,lambdaSamples);
    S.theory.med(ib,:) = E.med;
    S.theory.q025(ib,:) = E.q025;
    S.theory.q975(ib,:) = E.q975;
    S.theory.nMeteorologicalTimes(ib) = E.nMeteorologicalTimes;
    S.theory.nJointPredictions(ib) = E.nJointPredictions;
end

if nTime~=numel(u)
    error('FRAGMENT observation and meteorology dimensions disagree.')
end
end


function S = summarize_cdp(D,R,DpCurve_um,alpha,gamma,lambdaSamples)
ib = R.validPlotBins(1);
met = D.nominal.meteorology;
timeRetained = any(D.selection.nominalUseMask,2);
validMet = timeRetained & isfinite(met.ustar) & met.ustar>0 & ...
    isfinite(met.T) & met.T>0 & isfinite(met.p) & met.p>0 & ...
    ~isnan(met.L) & isfinite(met.z0) & met.z0>0 & ...
    isfinite(D.referenceHeight_m.displayReference) & ...
    D.referenceHeight_m.displayReference>0;
idx = find(validMet);

S.Ntime = nnz(timeRetained);
S.Dp_um = R.DpCDP_um;
S.obs.med = R.obs.med(ib,:);
S.obs.q25 = R.obs.q25(ib,:);
S.obs.q75 = R.obs.q75(ib,:);
S.lambda.samples = lambdaSamples(:);
S.lambda.median = median(lambdaSamples);
S.lambda.ci95 = prctile(lambdaSamples,[2.5 97.5]);
S.lambda.nSuccessful = numel(lambdaSamples);
S.lambda.nRequested = 500;
S.theory = theory_envelope(met.ustar,met.T,met.p,met.L,met.z0, ...
    D.referenceHeight_m.displayReference,idx,DpCurve_um,alpha,gamma, ...
    R.settings.rhop_kg_m3,R.settings.ARp,lambdaSamples);
S.wsReference = settling_velocity( ...
    median(met.T(timeRetained),'omitnan'), ...
    median(met.p(timeRetained),'omitnan'),DpCurve_um*1e-6,'bb', ...
    'rhop',R.settings.rhop_kg_m3,'ARp',R.settings.ARp);
S.wsReference = S.wsReference(:)';
S.wsReferenceConditions = struct( ...
    'T_K',median(met.T(timeRetained),'omitnan'), ...
    'p_Pa',median(met.p(timeRetained),'omitnan'), ...
    'method','bb','rhop_kg_m3',R.settings.rhop_kg_m3, ...
    'ARp',R.settings.ARp);
end


function E = theory_envelope(u,T,p,L,z0,zRef,idx,DpCurve_um, ...
        alpha,gamma,rhop,ARp,lambdaSamples)
lambdaSamples = lambdaSamples(isfinite(lambdaSamples));
if isempty(idx) || isempty(lambdaSamples)
    error('Cannot calculate theory envelope without times and lambda draws.')
end

DpCurve_m = DpCurve_um*1e-6;
nTime = numel(idx);
nLambda = numel(lambdaSamples);
nDp = numel(DpCurve_um);
base = NaN(nTime,nDp);
gain = NaN(nTime,nDp);
for ii = 1:nTime
    it = idx(ii);
    out = deposition_velocity(u(it),T(it),p(it),L(it), ...
        zRef(it),DpCurve_m,'vdMethod','MCzhang01','wsMethod','bb', ...
        'z0',z0(it),'rhop',rhop,'ARp',ARp, ...
        'alpha',alpha,'gamma',gamma);
    base(ii,:) = out.wd(:)';
    gain(ii,:) = u(it)*out.ws(:)';
end

lambdaMedian = median(lambdaSamples);
central = base + lambdaMedian*gain;
E.med = median(central,1,'omitnan');

joint = NaN(nTime*nLambda,nDp);
for ii = 1:nTime
    rows = (ii-1)*nLambda+(1:nLambda);
    joint(rows,:) = base(ii,:) + lambdaSamples(:)*gain(ii,:);
end
q = prctile(joint,[2.5 97.5],1);
E.q025 = q(1,:);
E.q975 = q(2,:);
E.nMeteorologicalTimes = nTime;
E.nLambdaRealizations = nLambda;
E.nJointPredictions = nTime*nLambda;
end


function [med,q25,q75] = summarize_observations(wdObs,inBin)
nBin = size(wdObs,2);
med = NaN(1,nBin);
q25 = NaN(1,nBin);
q75 = NaN(1,nBin);
for j = 1:nBin
    x = wdObs(inBin,j,:);
    x = x(isfinite(x) & x>0);
    if ~isempty(x)
        q = prctile(x,[25 50 75]);
        q25(j) = q(1);
        med(j) = q(2);
        q75(j) = q(3);
    end
end
end


function plot_envelope(ax,x,E,ib,color)
if nargin==4
    color = ib;
    ib = 1;
end
if size(E.q025,1)==1
    y1 = 100*E.q025;
    y2 = 100*E.q975;
    y = 100*E.med;
else
    y1 = 100*E.q025(ib,:);
    y2 = 100*E.q975(ib,:);
    y = 100*E.med(ib,:);
end
good = isfinite(x) & isfinite(y1) & isfinite(y2) & ...
    y1>0 & y2>0 & y2>=y1;
if nnz(good)>=2
    fill(ax,[x(good),fliplr(x(good))], ...
        [y1(good),fliplr(y2(good))],color, ...
        'FaceAlpha',0.16,'EdgeColor','none','HandleVisibility','off')
end
good = isfinite(x) & isfinite(y) & y>0;
plot(ax,x(good),y(good),'-','Color',color,'LineWidth',1.6, ...
    'HandleVisibility','off')
end


function plot_observations(ax,Dp,med,q25,q75,marker,color)
x = Dp(:);
y = med(:);
lo = q25(:);
hi = q75(:);
good = isfinite(x) & isfinite(y) & isfinite(lo) & isfinite(hi) & ...
    y>0 & lo>0 & hi>0;
x = x(good);
y = y(good);
lo = lo(good);
hi = hi(good);
elow = max(y-lo,1e-5);
ehigh = max(hi-y,1e-5);
elow = min(elow,0.99*y);
errorbar(ax,x,100*y,100*elow,100*ehigh,marker, ...
    'Color',color,'MarkerFaceColor','w','MarkerSize',6, ...
    'LineWidth',1,'CapSize',6,'HandleVisibility','off')
end


function validate_inputs(F,E,D,R,C,edges)
requiredF = {'nominal','mc','meteorology','q_15avg', ...
    'diam_mean_um','configuration'};
for i = 1:numel(requiredF)
    if ~isfield(F,requiredF{i})
        error('FRAGMENT input is missing "%s".',requiredF{i})
    end
end
if ~isfield(E,'fragment') || ~isfield(E.fragment,'lambdaMonteCarlo')
    error('FRAGMENT lambda Monte Carlo product is incomplete.')
end
M = E.fragment.lambdaMonteCarlo;
if size(M.lambda,2)~=numel(edges)-1 || M.nRequested~=500 || ...
        any(M.nSuccessful~=500)
    error('Expected four complete 500-member FRAGMENT lambda fits.')
end
if any(abs(M.windowStart(:)-edges(1:end-1)')>1e-10) || ...
        any(abs(M.windowEnd(:)-edges(2:end)')>1e-10)
    error('FRAGMENT lambda-fit windows do not match the plotted bins.')
end
if D.selection.timeCount~=26 || D.selection.cellCount~=56 || ...
        nnz(any(D.selection.nominalUseMask,2))~=26
    error('Updated CDP selection must contain 26 times and 56 cells.')
end
if R.Ntime(R.validPlotBins(1))~=26
    error('CDP plotting statistics do not report retained N=26.')
end
if C.monteCarlo.nSuccessful~=500 || nnz(C.monteCarlo.valid)~=500
    error('Expected 500/500 successful updated CDP lambda fits.')
end
end
