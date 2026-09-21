%% PLOT_FRAGMENT_ZHANG_PARAMETER_ENSEMBLE
% Rebuild both panels of the FRAGMENT deposition-velocity comparison.
% Uses the paired alpha/gamma Monte Carlo ensemble fitted to the Zhang et
% al. (2014) laboratory measurements. Lambda is zero and the impaction
% exponent is fixed at two. Lines and shading are the pooled median and
% IQR across FRAGMENT meteorological states and paired parameter draws:
% Zhang fits in (a), FRAGMENT alpha/gamma fits in (b). No refitting is done.
% Titles use the full accepted ensembles for medians and 95% intervals.

clearvars
close all
clc

fprintf('Running plot_fragment_zhang_parameter_ensemble\n')

%% Paths and fixed physical settings
thisDir = fileparts(mfilename('fullpath'));
fragmentDir = fileparts(thisDir);
srcDir = fileparts(fragmentDir);
projectRoot = fileparts(srcDir);
addpath(srcDir)

plotDataFile = fullfile(projectRoot,'data','fragment','derived', ...
    'Fragment_wd_plot_data_MC.mat');
fitReadyFile = fullfile(projectRoot,'data','fragment','derived', ...
    'Fragment_fit_ready_MC.mat');
zhangFitFile = fullfile(projectRoot,'data','zhang2014','results', ...
    'zhang14_fit_default_cliftgauvin_n2000.mat');
fragmentFitFile = fullfile(projectRoot,'data','fragment','results', ...
    'Fragment_fit_alpha_gamma_MC.mat');
outputDir = fullfile(projectRoot,'output');
baseName = 'FragmentFit';
pngFile = fullfile(outputDir,[baseName '.png']);
pdfFile = fullfile(outputDir,[baseName '.pdf']);
figFile = fullfile(outputDir,[baseName '.fig']);
dataFile = fullfile(outputDir,[baseName '_data.mat']);

vdMethod = 'MCzhang01';
wsMethod = 'bb';
rhop = 2500;
ARp = 1.46;
lambda = 0;
beta = 2;
NsampMC = 500; % Same cap as ZhangFit; preserve alpha/gamma pairing.
rngSeed = 2;
rng(rngSeed,'twister')

%% Load observations, meteorology, and both paired parameter ensembles
Q = load(plotDataFile,'plotData');
P = Q.plotData;
S = load(fitReadyFile,'meteorology','q_15avg','configuration');
Z = load(zhangFitFile,'results');

F = load(fragmentFitFile,'fit');
ensembles = {Z.results.nonsticky.mc,F.fit.monteCarlo};
acceptedMasks = {Z.results.nonsticky.accepted,F.fit.monteCarlo.valid};
parameterDraws = cell(1,2);
parameterSummary = repmat(struct('alphaMedian',[], 'alphaIQR',[], ...
    'alpha95',[], 'gammaMedian',[], 'gammaIQR',[], 'gamma95',[], ...
    'nPairs',[]),1,2);
for ip = 1:2
    M = ensembles{ip};
    useParameter = logical(acceptedMasks{ip}(:)) & ...
        isfinite(M.alpha(:)) & M.alpha(:)>0 & ...
        isfinite(M.gamma(:)) & M.gamma(:)>0;
    alpha = M.alpha(useParameter);
    gamma = M.gamma(useParameter);
    if isempty(alpha)
        error('Panel %d has no accepted alpha/gamma pairs.',ip)
    end
    parameterSummary(ip) = struct( ...
        'alphaMedian',median(alpha), ...
        'alphaIQR',prctile(alpha,[25 75]), ...
        'alpha95',prctile(alpha,[2.5 97.5]), ...
        'gammaMedian',median(gamma), ...
        'gammaIQR',prctile(gamma,[25 75]), ...
        'gamma95',prctile(gamma,[2.5 97.5]), ...
        'nPairs',numel(alpha));
    idx = randperm(numel(alpha),min(NsampMC,numel(alpha)));
    parameterDraws{ip} = struct('alpha',alpha(idx),'gamma',gamma(idx));
    fprintf('Panel %d: %d accepted pairs; using %d for envelopes.\n', ...
        ip,numel(alpha),numel(idx))
end

%% Match the current FRAGMENT selection and four ustar bins
u = S.meteorology.ustar(:);
T = S.meteorology.T_K(:);
p = S.meteorology.p_Pa(:);
L = S.meteorology.L_m(:);
z0 = S.meteorology.z0_m(:);
zRef = S.configuration.zReference_m;
edges = P.ustarEdges(:);
binID = discretize(u,edges);
timeEligible = isnan(S.q_15avg(:)) & isfinite(u) & u>0 & u<0.20;
binID(~timeEligible) = NaN;

Ntime = arrayfun(@(ib)nnz(binID==ib),(1:numel(edges)-1)');
if ~isequal(Ntime(:),P.Ntime(:))
    error('Reconstructed ustar-bin counts do not match the cached figure.')
end

Dp_m = P.DpCurve_um(:)*1e-6;
validMet = timeEligible & isfinite(T) & T>0 & isfinite(p) & p>0 & ...
    ~isnan(L) & isfinite(z0) & z0>0;
base = build_base_cache(u,T,p,L,z0,zRef,Dp_m,validMet, ...
    vdMethod,wsMethod,rhop,ARp);


%% Pool each parameter ensemble with all states in each ustar bin
theoryPanels = repmat(struct('med',[], 'q25',[], 'q75',[], ...
    'Npool',[]),1,2);
for ip = 1:2
    alpha = parameterDraws{ip}.alpha;
    gamma = parameterDraws{ip}.gamma;
    validate_evaluator(base,u,T,p,L,z0,zRef,Dp_m,validMet, ...
        alpha,gamma,vdMethod,wsMethod,rhop,ARp);
    nUstarBin = numel(edges)-1;
    nCurve = numel(Dp_m);
    theory.med = NaN(nUstarBin,nCurve);
    theory.q25 = NaN(nUstarBin,nCurve);
    theory.q75 = NaN(nUstarBin,nCurve);
    theory.Npool = zeros(nUstarBin,1);

    for ib = P.validPlotBins(:)'
        idx = find(binID==ib & base.valid);
        nState = numel(idx);
        nParameter = numel(alpha);
        pool = NaN(nState*nParameter,nCurve,'single');
        for ii = 1:nState
            rows = (ii-1)*nParameter+(1:nParameter);
            pool(rows,:) = single(evaluate_ensemble( ...
                base,idx(ii),u,alpha,gamma,beta)');
        end
        q = prctile(pool,[25 50 75],1);
        theory.q25(ib,:) = double(q(1,:));
        theory.med(ib,:) = double(q(2,:));
        theory.q75(ib,:) = double(q(3,:));
        theory.Npool(ib) = size(pool,1);
        fprintf(['  %.2f <= ustar < %.2f: %d states x %d parameter ' ...
            'pairs = %d curves\n'],edges(ib),edges(ib+1), ...
            nState,nParameter,theory.Npool(ib))
        clear pool
    end

    theoryPanels(ip) = theory;
end

%% Plot both panels in the existing FRAGMENT style
colors = lines(numel(P.validPlotBins));
fig = figure(98);
clf(fig)
set(fig,'Color','w')
panelLabels = {'a.','b.'};
for ip = 1:2
    theory = theoryPanels(ip);
    summary = parameterSummary(ip);
    % Original FRAGMENT geometry; leave the figure-window size unchanged.
    ax = subplot(2,2,ip,'Parent',fig);
    ax.Position(2) = 0.25;
    ax.Position(4) = 0.67;
    hold(ax,'on')
    box(ax,'on')

    hBin = gobjects(numel(P.validPlotBins),1);
    binLabels = strings(numel(P.validPlotBins),1);
    for kk = 1:numel(P.validPlotBins)
        ib = P.validPlotBins(kk);
        color = colors(kk,:);

        plot_observations(ax,P.DpFragment_um,P.obs.med(ib,:), ...
            P.obs.q25(ib,:),P.obs.q75(ib,:),color)
        plot_envelope(ax,P.DpCurve_um,theory,ib,color)

        binLabels(kk) = sprintf( ...
            '$u_*=(%.2f$--$%.2f)$ m s$^{-1}$, $N=%d$', ...
            edges(ib),edges(ib+1),P.Ntime(ib));
        hBin(kk) = plot(ax,NaN,NaN,'-','Color',color, ...
            'LineWidth',1.6,'DisplayName',binLabels(kk));
    end

    hWs = plot(ax,P.DpCurve_um,100*P.wsMedian,'--k', ...
        'LineWidth',1.4,'DisplayName','$w_s$');
    set(ax,'XScale','log','YScale','log', ...
        'XLim',[0.3 20],'XTick',[0.3 1 10 20], ...
        'YLim',[3e-3 5e1],'FontSize',12, ...
        'TickLabelInterpreter','latex')
    xlabel(ax,'$D_p$ ($\mu$m)','Interpreter','latex')
    if ip==1
        ylabel(ax,'$w_d$ (cm s$^{-1}$)','Interpreter','latex')
    end
    title(ax,{sprintf('$\\alpha=%.2g\\;[%.2g,\\,%.2g]$', ...
        summary.alphaMedian,summary.alpha95), ...
        sprintf('$\\gamma=%.3f\\;[%.3f,\\,%.3f]$', ...
        summary.gammaMedian,summary.gamma95)}, ...
        'Interpreter','latex','FontWeight','normal','FontSize',11)
    if ip==1
        legend(ax,[hBin;hWs],[binLabels;'$w_s$'], ...
            'Location','northwest','Interpreter','latex','FontSize',9)
    end
    grid(ax,'on')
    text(ax,-0.11,1.01,panelLabels{ip},'Units','normalized', ...
        'HorizontalAlignment','left','VerticalAlignment','bottom', ...
        'FontSize',12,'FontWeight','bold','Interpreter','latex')

end

if ~isfolder(outputDir)
    mkdir(outputDir)
end
settings = struct('lambda_s_m',lambda,'impactionExponentBeta',beta, ...
    'vdMethod',vdMethod,'wsMethod',wsMethod,'rhop_kg_m3',rhop, ...
    'ARp',ARp,'parameterPooling', ...
    ['sampled paired alpha/gamma realizations pooled with ' ...
    'all eligible FRAGMENT meteorological states in each ustar bin'], ...
    'theorySummary','median and IQR of pooled theory curves', ...
    'observationSource',plotDataFile,'fitReadySource',fitReadyFile, ...
    'zhangFitSource',zhangFitFile,'fragmentFitSource',fragmentFitFile, ...
    'NsampMC',NsampMC,'rngSeed',rngSeed, ...
    'panelSources',{{'Zhang','FRAGMENT'}});
save(dataFile,'theoryPanels','parameterSummary','parameterDraws','settings','-v7.3')
exportgraphics(fig,pngFile,'Resolution',300)
% Keep the image name already referenced by the measurements main.tex.
exportgraphics(fig,fullfile(outputDir,'Fragment_wd_default_fit.jpg'), ...
    'Resolution',300)
exportgraphics(fig,pdfFile,'ContentType','vector')
savefig(fig,figFile)
fprintf('Saved FRAGMENT two-panel parameter-ensemble figure:\n  %s\n  %s\n  %s\n', ...
    pngFile,pdfFile,figFile)
fprintf('Saved plotted theory data:\n  %s\n',dataFile)


function B = build_base_cache(u,T,p,L,z0,zRef,Dp,valid, ...
        vdMethod,wsMethod,rhop,ARp)
nTime = numel(u);
nCurve = numel(Dp);
for name = {'ws','ra','Sc','St','R1'}
    B.(name{1}) = NaN(nTime,nCurve);
end
B.valid = false(nTime,1);

for it = find(valid(:))'
    out = deposition_velocity(u(it),T(it),p(it),L(it),zRef,Dp, ...
        'vdMethod',vdMethod,'wsMethod',wsMethod,'z0',z0(it), ...
        'rhop',rhop,'ARp',ARp,'alpha',1,'gamma',1);
    St = out.ws(:)*u(it)^2/(9.81*out.air.nu);
    Sc = 1./out.Eb(:);
    Eim = (St./(1+St)).^2;
    R1 = 1./(3*u(it)*out.rb(:).*(out.Eb(:)+Eim));

    B.ws(it,:) = out.ws(:)';
    B.ra(it,:) = out.ra(:)';
    B.Sc(it,:) = Sc(:)';
    B.St(it,:) = St(:)';
    B.R1(it,:) = R1(:)';
    B.valid(it) = all(isfinite([out.ws(:);out.ra(:);Sc(:);St(:);R1(:)]));
end
end


function wd = evaluate_ensemble(B,it,u,alpha,gamma,beta)
ws = B.ws(it,:)';
ra = B.ra(it,:)';
Sc = B.Sc(it,:)';
St = B.St(it,:)';
R1 = B.R1(it,:)';
alpha = alpha(:)';
gamma = gamma(:)';

Eb = exp(-log(Sc)*gamma);
Eim = (St./(alpha+St)).^beta;
rb = 1./(3*u(it)*(Eb+Eim).*R1);
wd = ws./(1-exp(-ws.*(ra+rb)));
end


function validate_evaluator(B,u,T,p,L,z0,zRef,Dp,valid, ...
        alpha,gamma,vdMethod,wsMethod,rhop,ARp)
it = find(valid & B.valid,1,'first');
ip = unique([1 ceil(numel(alpha)/2) numel(alpha)]);
fast = evaluate_ensemble(B,it,u,alpha(ip),gamma(ip),2);
maxRelativeError = 0;
for j = 1:numel(ip)
    out = deposition_velocity(u(it),T(it),p(it),L(it),zRef,Dp, ...
        'vdMethod',vdMethod,'wsMethod',wsMethod,'z0',z0(it), ...
        'rhop',rhop,'ARp',ARp, ...
        'alpha',alpha(ip(j)),'gamma',gamma(ip(j)));
    relativeError = max(abs(fast(:,j)-out.wd(:))./max(abs(out.wd(:)),eps));
    maxRelativeError = max(maxRelativeError,relativeError);
end
fprintf('Maximum vectorized-model relative error: %.3g\n',maxRelativeError)
if maxRelativeError>1e-10
    error('Vectorized ensemble evaluator does not match deposition_velocity.')
end
end


function plot_observations(ax,Dp,med,q25,q75,color)
use = isfinite(Dp) & isfinite(med) & med>0 & ...
    isfinite(q25) & q25>0 & isfinite(q75) & q75>=med;
x = Dp(use);
y = 100*med(use);
lower = 100*(med(use)-q25(use));
upper = 100*(q75(use)-med(use));
errorbar(ax,x,y,lower,upper,'s','Color',color, ...
    'MarkerFaceColor','w','MarkerSize',5.5,'LineWidth',1, ...
    'CapSize',6,'HandleVisibility','off')
end


function plot_envelope(ax,Dp,E,ib,color)
q25 = 100*E.q25(ib,:);
med = 100*E.med(ib,:);
q75 = 100*E.q75(ib,:);
use = isfinite(Dp) & isfinite(q25) & q25>0 & ...
    isfinite(med) & med>0 & isfinite(q75) & q75>=q25;
x = Dp(use);
fill(ax,[x fliplr(x)],[q25(use) fliplr(q75(use))],color, ...
    'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off')
plot(ax,x,med(use),'-','Color',color,'LineWidth',1.6, ...
    'HandleVisibility','off')
end
