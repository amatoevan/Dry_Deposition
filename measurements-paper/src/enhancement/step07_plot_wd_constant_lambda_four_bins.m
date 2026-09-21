%% STEP07_PLOT_WD_CONSTANT_LAMBDA_FOUR_BINS
% Matched four-bin FRAGMENT and OCOTILLO CDP comparison using one fixed
% enhancement coefficient for every theoretical deposition-velocity curve.

clearvars
close all
clc

fprintf('Running step07_plot_wd_constant_lambda_four_bins\n')

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
anzaSelectedFile = fullfile(projectRoot,'data','ocotillo','derived', ...
    'eddypro_10min_fdiff_pc10_fc', ...
    'step01b_cdp_mc_selected_10min_7bins_paired_reference_heights.mat');
anzaStatsFile = fullfile(projectRoot,'data','ocotillo','results', ...
    'eddypro_10min_fdiff_pc10_fc', ...
    'ocotillo_wd_10min_7bins_default_fragment_lambda_at_2p1m_statistics.mat');
binSpecificStatsFile = fullfile(projectRoot,'data','enhancement','results', ...
    'enhance_wd_bin_specific_lambda_statistics.mat');

figDir = fullfile(projectRoot,'figures','enhancement');
dataFile = fullfile(projectRoot,'data','enhancement','results', ...
    'enhance_wd_constant_lambda_four_fragment_bins_cdp_statistics.mat');
baseName = 'Figure_wd_constant_lambda_four_fragment_bins_cdp';
pngFile = fullfile(figDir,[baseName '.png']);
pdfFile = fullfile(figDir,[baseName '.pdf']);
figFile = fullfile(figDir,[baseName '.fig']);

inputFiles = {fragmentFile;anzaSelectedFile;anzaStatsFile; ...
    binSpecificStatsFile};
for i = 1:numel(inputFiles)
    if ~isfile(inputFiles{i})
        error('Required input not found:\n%s',inputFiles{i})
    end
end

F = load(fragmentFile);
D = load(anzaSelectedFile);
A = load(anzaStatsFile,'R');
P = load(binSpecificStatsFile,'S');

alphaFragment = A.R.settings.alphaFragment;
gammaFragment = A.R.settings.gammaFragment;
lambdaConstant = A.R.lambda.fragmentReference;
validate_inputs(F,D,A.R,P.S,fragmentUstarEdges,lambdaConstant)

%% Constant-lambda theoretical envelopes
fragmentTheory = fragment_constant_theory(F,fragmentUstarEdges, ...
    DpCurve_um,alphaFragment,gammaFragment,lambdaConstant);
cdpTheory = cdp_constant_theory(D,A.R,DpCurve_um, ...
    alphaFragment,gammaFragment,lambdaConstant);

%% Plot
colors = [ ...
    0.0000 0.4470 0.7410
    0.8500 0.3250 0.0980
    0.9290 0.6940 0.1250
    0.4660 0.6740 0.1880
    0.4940 0.1840 0.5560];

fig = figure(91);
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
    plot_envelope(ax,DpCurve_um,fragmentTheory,ib,colors(ib,:))
    plot_observations(ax,F.diam_mean_um, ...
        P.S.fragment.obs.med(ib,:),P.S.fragment.obs.q25(ib,:), ...
        P.S.fragment.obs.q75(ib,:),'s',colors(ib,:))
    labels(ib) = sprintf( ...
        'FRAGMENT: $u_*=(%.2f$--$%.2f)$, $N=%d$', ...
        fragmentUstarEdges(ib),fragmentUstarEdges(ib+1), ...
        P.S.fragment.Ntime(ib));
    h(ib) = plot(ax,NaN,NaN,'s-','Color',colors(ib,:), ...
        'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1.6, ...
        'DisplayName',labels(ib));
end

plot_envelope(ax,DpCurve_um,cdpTheory,1,colors(5,:))
plot_observations(ax,P.S.cdp.Dp_um,P.S.cdp.obs.med, ...
    P.S.cdp.obs.q25,P.S.cdp.obs.q75,'o',colors(5,:))
labels(5) = sprintf('OCOTILLO CDP: $N=%d$',P.S.cdp.Ntime);
h(5) = plot(ax,NaN,NaN,'o-','Color',colors(5,:), ...
    'MarkerFaceColor','w','MarkerSize',5,'LineWidth',1.6, ...
    'DisplayName',labels(5));

labels(6) = '$w_s$';
h(6) = plot(ax,DpCurve_um,100*P.S.cdp.wsReference,'--k', ...
    'LineWidth',1.5,'DisplayName',labels(6));

set(ax,'XScale','log','YScale','log', ...
    'XLim',xLimits,'XTick',xTicks,'YLim',yLimits, ...
    'FontSize',12,'TickLabelInterpreter','latex')
xlabel(ax,'$D_p$ ($\mu$m)','Interpreter','latex')
ylabel(ax,'$w_d$ (cm s$^{-1}$)','Interpreter','latex')
title(ax,sprintf( ...
    '$\\alpha=%.2f$, $\\gamma=%.1f$, $\\lambda=%.1f$ s m$^{-1}$', ...
    alphaFragment,gammaFragment,lambdaConstant), ...
    'Interpreter','latex','FontWeight','normal','FontSize',12)
legend(ax,h,labels,'Location','northwest', ...
    'Interpreter','latex','FontSize',7.2)
text(ax,0.98,0.025, ...
    'Shading: 95\% range across retained meteorological times', ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','bottom','Interpreter','latex','FontSize',7.2)
grid(ax,'on')

%% Save figure and plotting statistics
S = struct;
S.workflow = 'enhancement/constant-lambda-four-fragment-bins-cdp';
S.created = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss Z'));
S.inputs = struct('fragmentData',fragmentFile, ...
    'anzaSelectedData',anzaSelectedFile, ...
    'anzaPlotStatistics',anzaStatsFile, ...
    'observationSummaries',binSpecificStatsFile);
S.parameters = struct('alpha',alphaFragment,'gamma',gammaFragment, ...
    'lambda_s_m',lambdaConstant,'DpCurve_um',DpCurve_um, ...
    'fragmentUstarEdges',fragmentUstarEdges);
S.uncertainty = ['central line: median prediction across retained ' ...
    'meteorological times; shading: 2.5th-97.5th percentiles across ' ...
    'retained meteorological times; lambda is fixed without uncertainty'];
S.fragment = struct('Ntime',P.S.fragment.Ntime, ...
    'obs',P.S.fragment.obs,'theory',fragmentTheory);
S.cdp = struct('Ntime',P.S.cdp.Ntime,'Dp_um',P.S.cdp.Dp_um, ...
    'obs',P.S.cdp.obs,'theory',cdpTheory, ...
    'wsReference',P.S.cdp.wsReference, ...
    'wsReferenceConditions',P.S.cdp.wsReferenceConditions);
S.figureFiles = struct('png',pngFile,'pdf',pdfFile,'fig',figFile);

if saveFig
    if ~isfolder(figDir)
        mkdir(figDir)
    end
    exportgraphics(ax,pngFile,'Resolution',300)
    exportgraphics(ax,pdfFile,'ContentType','vector')
    savefig(fig,figFile)
    save(dataFile,'S','-v7.3')
    fprintf('Saved constant-lambda figure:\n  %s\n  %s\n  %s\n', ...
        pngFile,pdfFile,figFile)
    fprintf('Saved plotting statistics:\n  %s\n',dataFile)
end

fprintf('Constant lambda: %.6g s m^-1\n',lambdaConstant)
fprintf('FRAGMENT retained-time N: %s\n', ...
    mat2str(P.S.fragment.Ntime(:)'))
fprintf('CDP retained-time N: %d\n',P.S.cdp.Ntime)


%% Local functions
function E = fragment_constant_theory(F,edges,DpCurve_um,alpha,gamma,lambda)
u = F.meteorology.ustar(:);
T = F.meteorology.T_K(:);
p = F.meteorology.p_Pa(:);
L = F.meteorology.L_m(:);
z0 = F.meteorology.z0_m(:);
zRef = repmat(F.configuration.zReference_m,numel(u),1);
nDiameter = numel(F.diam_mean_um);
timeEligible = isnan(F.q_15avg(:)) & ...
    isfinite(u) & u>0 & u<edges(end);
nominalMask = repmat(timeEligible,1,nDiameter) & ...
    isfinite(F.nominal.Fdiff) & F.nominal.Fdiff<0 & ...
    isfinite(F.nominal.wd_obs) & F.nominal.wd_obs>0;
timeRetained = any(nominalMask,2);
validMet = timeRetained & isfinite(T) & T>0 & isfinite(p) & p>0 & ...
    ~isnan(L) & isfinite(z0) & z0>0;
binID = discretize(u,edges);
nGroup = numel(edges)-1;
E.med = NaN(nGroup,numel(DpCurve_um));
E.q025 = NaN(nGroup,numel(DpCurve_um));
E.q975 = NaN(nGroup,numel(DpCurve_um));
E.nMeteorologicalTimes = zeros(nGroup,1);
for ib = 1:nGroup
    idx = find(binID==ib & validMet);
    G = constant_theory(u,T,p,L,z0,zRef,idx,DpCurve_um, ...
        alpha,gamma,2500,1.46,lambda);
    E.med(ib,:) = G.med;
    E.q025(ib,:) = G.q025;
    E.q975(ib,:) = G.q975;
    E.nMeteorologicalTimes(ib) = G.nMeteorologicalTimes;
end
end


function E = cdp_constant_theory(D,R,DpCurve_um,alpha,gamma,lambda)
met = D.nominal.meteorology;
timeRetained = any(D.selection.nominalUseMask,2);
validMet = timeRetained & isfinite(met.ustar) & met.ustar>0 & ...
    isfinite(met.T) & met.T>0 & isfinite(met.p) & met.p>0 & ...
    ~isnan(met.L) & isfinite(met.z0) & met.z0>0 & ...
    isfinite(D.referenceHeight_m.displayReference) & ...
    D.referenceHeight_m.displayReference>0;
idx = find(validMet);
E = constant_theory(met.ustar,met.T,met.p,met.L,met.z0, ...
    D.referenceHeight_m.displayReference,idx,DpCurve_um,alpha,gamma, ...
    R.settings.rhop_kg_m3,R.settings.ARp,lambda);
end


function E = constant_theory(u,T,p,L,z0,zRef,idx,DpCurve_um, ...
        alpha,gamma,rhop,ARp,lambda)
if isempty(idx)
    error('Cannot calculate a theory envelope without retained times.')
end
DpCurve_m = DpCurve_um*1e-6;
pool = NaN(numel(idx),numel(DpCurve_um));
for ii = 1:numel(idx)
    it = idx(ii);
    out = deposition_velocity(u(it),T(it),p(it),L(it), ...
        zRef(it),DpCurve_m,'vdMethod','MCzhang01','wsMethod','bb', ...
        'z0',z0(it),'rhop',rhop,'ARp',ARp, ...
        'alpha',alpha,'gamma',gamma);
    pool(ii,:) = out.wd(:)' + lambda*u(it)*out.ws(:)';
end
q = prctile(pool,[2.5 50 97.5],1);
E.q025 = q(1,:);
E.med = q(2,:);
E.q975 = q(3,:);
E.nMeteorologicalTimes = numel(idx);
end


function plot_envelope(ax,x,E,ib,color)
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


function validate_inputs(F,D,R,P,edges,lambda)
requiredF = {'nominal','meteorology','q_15avg', ...
    'diam_mean_um','configuration'};
for i = 1:numel(requiredF)
    if ~isfield(F,requiredF{i})
        error('FRAGMENT input is missing "%s".',requiredF{i})
    end
end
if D.selection.timeCount~=26 || D.selection.cellCount~=56 || ...
        nnz(any(D.selection.nominalUseMask,2))~=26
    error('Updated CDP selection must contain 26 times and 56 cells.')
end
if ~isequal(P.fragment.Ntime(:)',[103 218 255 86]) || ...
        P.cdp.Ntime~=26
    error('Observation summaries do not contain the expected retained N.')
end
if numel(edges)~=5 || abs(lambda-11.7)>1e-10 || ...
        abs(R.settings.lambdaFragmentReference_s_m-lambda)>1e-10
    error('Expected four bins and constant lambda = 11.7 s m^-1.')
end
end
