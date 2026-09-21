function figureFiles = plot_default_vs_fragment_lambda0(options)
%PLOT_DEFAULT_VS_FRAGMENT_LAMBDA0 Compare two lambda=0 parameter sets.
% Uses the saved plotting statistics so observations and theory envelopes
% are exactly those from the validated retained-time workflow.
% plot_default_vs_fragment_lambda0(ManuscriptDefault=true) regenerates
% Anza_wd_default with paired Zhang parameter and meteorological-state IQRs.

arguments
    options.DerivedFile (1,:) char = ''
    options.StatsFile (1,:) char = ''
    options.FigureDir (1,:) char = ''
    options.SaveFigures (1,1) logical = true
    options.ManuscriptDefault (1,1) logical = false
end

thisDir = fileparts(mfilename('fullpath'));
srcDir = fileparts(thisDir);
projectRoot = fileparts(srcDir);
addpath(srcDir)

workflowName = 'eddypro_10min_fdiff_pc10_fc';
if isempty(options.DerivedFile)
    options.DerivedFile = fullfile(projectRoot,'data','ocotillo', ...
        'derived',workflowName, ...
        'step01b_cdp_mc_selected_10min_7bins_paired_reference_heights.mat');
end
if isempty(options.StatsFile)
    options.StatsFile = fullfile(projectRoot,'data','ocotillo', ...
        'results',workflowName, ...
        ['ocotillo_wd_10min_7bins_default_fragment_' ...
         'lambda_at_2p1m_statistics.mat']);
end
if isempty(options.FigureDir)
    if options.ManuscriptDefault
        options.FigureDir = fullfile(projectRoot,'output');
    else
        options.FigureDir = fullfile(projectRoot,'figures','ocotillo', ...
            'sensitivity',workflowName);
    end
end

S = load(options.DerivedFile,'nominal','referenceHeight_m');
Q = load(options.StatsFile,'R');
R = Q.R;
validate_inputs(S,R)

alphaDefault = 2.4;
gammaDefault = 0.37;
alphaFragment = 0.06;
gammaFragment = 0.434742;
lambda = 0;
xLimits = [2 50];
xTicks = [2 10 50];
yLimits = [1e-1 1e2];
rhop = 2680;
ARp = 1.5;
wsMethod = 'bb';

u = S.nominal.meteorology.ustar(:);
T = S.nominal.meteorology.T(:);
p = S.nominal.meteorology.p(:);
binID = discretize(u,R.ustarEdges);
wsMedian = settling_velocity(median(T,'omitnan'),median(p,'omitnan'), ...
    R.DpCurve_um*1e-6,wsMethod,'rhop',rhop,'ARp',ARp);
wsMedian = wsMedian(:)';

theory = {R.theory.default.lambda0,R.theory.fragment.lambda0};
panelTitles = { ...
    parameter_title(alphaDefault,gammaDefault,lambda), ...
    parameter_title(alphaFragment,gammaFragment,lambda)};
if options.ManuscriptDefault
    yLimits = [1e-1 2e2]; % Preserve the manuscript range and upper error bar.
    zhangFitFile = fullfile(projectRoot,'data','zhang2014','results', ...
        'zhang14_fit_default_cliftgauvin_n2000.mat');
    Z = load(zhangFitFile,'results');
    M = Z.results.nonsticky.mc;
    use = logical(Z.results.nonsticky.accepted(:)) & ...
        isfinite(M.alpha(:)) & M.alpha(:)>0 & ...
        isfinite(M.gamma(:)) & M.gamma(:)>0;
    alpha = M.alpha(use);
    gamma = M.gamma(use);
    assert(~isempty(alpha),'No accepted Zhang parameter pairs.')
    parameterSummary = struct('alphaMedian',median(alpha), ...
        'alpha95',prctile(alpha,[2.5 97.5]), ...
        'gammaMedian',median(gamma),'gamma95',prctile(gamma,[2.5 97.5]));
    rng(2,'twister')
    idx = randperm(numel(alpha),min(500,numel(alpha)));
    parameterDraws = struct('alpha',alpha(idx),'gamma',gamma(idx));
    theory = {zhang_met_envelope(S,R,parameterDraws,rhop,ARp,wsMethod)};
    panelTitles = {{ ...
        sprintf('$\\alpha=%.2g\\;[%.2g,\\,%.2g]$', ...
            parameterSummary.alphaMedian,parameterSummary.alpha95), ...
        sprintf('$\\gamma=%.3f\\;[%.3f,\\,%.3f]$', ...
            parameterSummary.gammaMedian,parameterSummary.gamma95)}};
end
colors = lines(numel(R.validPlotBins));

fig = figure('Color','w');
for ip = 1:numel(theory)
    ax = subplot(2,2,ip,'Parent',fig);
    hold(ax,'on')
    box(ax,'on')
    hBin = gobjects(numel(R.validPlotBins),1);
    binLabels = strings(numel(R.validPlotBins),1);

    for kk = 1:numel(R.validPlotBins)
        ib = R.validPlotBins(kk);
        color = colors(kk,:);
        plot_observations(ax,R.DpCDP_um,R.obs.med(ib,:), ...
            R.obs.q25(ib,:),R.obs.q75(ib,:),'o',color);
        plot_envelope(ax,R.DpCurve_um,theory{ip},ib,color);

        medianUstar = median(u(binID==ib),'omitnan');
        binLabels(kk) = sprintf( ...
            '$\\widetilde{u_*}=%s$ m s$^{-1}$', ...
            format_parameter_value(medianUstar));
        if options.ManuscriptDefault
            binLabels(kk) = sprintf( ...
                'Median $u_*=%.1f$ m s$^{-1}$, $N=%d$', ...
                medianUstar,R.Ntime(ib));
        end
        hBin(kk) = plot(ax,NaN,NaN,'-','Color',color, ...
            'LineWidth',1.5,'DisplayName',binLabels(kk));
    end

    wsLabel = '$w_s$';
    hWs = plot(ax,R.DpCurve_um,100*wsMedian,'--k','LineWidth',1.5, ...
        'DisplayName',wsLabel);
    set(ax,'XScale','log','YScale','log', ...
        'XLim',xLimits,'XTick',xTicks,'YLim',yLimits, ...
        'FontSize',12,'TickLabelInterpreter','latex')
    xlabel(ax,'$D_p$ ($\mu$m)','Interpreter','latex')
    if ip==1
        ylabel(ax,'$w_d$ (cm s$^{-1}$)','Interpreter','latex')
    end
    title(ax,panelTitles{ip},'Interpreter','latex', ...
        'FontWeight','normal')
    legend(ax,[hBin;hWs],[binLabels;wsLabel], ...
        'Location','southeast','Interpreter','latex','FontSize',9)
    grid(ax,'on')
    if ~options.ManuscriptDefault
    text(ax,0.02,0.97,char('a'+ip-1)+".",'Units','normalized', ...
        'HorizontalAlignment','left','VerticalAlignment','top', ...
        'FontSize',12,'FontWeight','bold','Interpreter','latex')
    end

    panelPosition = ax.Position;
    panelStretch = 1.7;
    ax.Position = [ ...
        panelPosition(1), ...
        panelPosition(2)-(panelStretch-1)*panelPosition(4), ...
        panelPosition(3), ...
        panelStretch*panelPosition(4)];
end

baseName = 'Figure_wd_default_vs_fragment_lambda0_at_2p1m';
if options.ManuscriptDefault, baseName = 'Anza_wd_default'; end
figureFiles = struct( ...
    'png',fullfile(options.FigureDir,[baseName '.png']), ...
    'pdf',fullfile(options.FigureDir,[baseName '.pdf']), ...
    'fig',fullfile(options.FigureDir,[baseName '.fig']));
if options.SaveFigures
    if ~isfolder(options.FigureDir), mkdir(options.FigureDir), end
    if options.ManuscriptDefault
        exportTarget = ax; % Crop to the original single stretched panel.
    else
        exportTarget = fig;
    end
    exportgraphics(exportTarget,figureFiles.png,'Resolution',300)
    exportgraphics(exportTarget,figureFiles.pdf,'ContentType','vector')
    savefig(fig,figureFiles.fig)
end

if options.ManuscriptDefault
    if options.SaveFigures
        exportgraphics(ax,fullfile(options.FigureDir,[baseName '.jpg']), ...
            'Resolution',300)
        settings = struct('zhangFitSource',zhangFitFile, ...
            'derivedSource',options.DerivedFile,'statsSource',options.StatsFile, ...
            'rngSeed',2,'NsampMC',500,'lambda',0, ...
            'pooling','all retained meteorological states x paired Zhang draws');
        save(fullfile(options.FigureDir,[baseName '_data.mat']), ...
            'theory','parameterSummary','parameterDraws','settings')
    end
    fprintf('Saved Anza parameter-ensemble figure: %s\n',figureFiles.fig)
    return
end

R.figureFiles.defaultVsFragmentLambda0 = figureFiles;
R.additionalComparisons.defaultVsFragmentLambda0 = struct( ...
    'left',struct('alpha',alphaDefault,'gamma',gammaDefault, ...
        'lambda_s_m',lambda), ...
    'right',struct('alpha',alphaFragment,'gamma',gammaFragment, ...
        'lambda_s_m',lambda), ...
    'usesRetainedTimeCount',sum(R.Ntime));
save(options.StatsFile,'R','-v7.3')
fprintf('Saved default-versus-FRAGMENT lambda=0 comparison:\n%s\n', ...
    figureFiles.fig)
end


function validate_inputs(S,R)
requiredR = {'ustarEdges','validPlotBins','DpCDP_um','DpCurve_um', ...
    'obs','theory','Ntime','figureFiles'};
for i = 1:numel(requiredR)
    if ~isfield(R,requiredR{i})
        error('Plotting statistics lack R.%s.',requiredR{i})
    end
end
if ~isfield(R.theory,'default') || ...
        ~isfield(R.theory.default,'lambda0') || ...
        ~isfield(R.theory,'fragment') || ...
        ~isfield(R.theory.fragment,'lambda0')
    error('Plotting statistics lack the required lambda=0 envelopes.')
end
if ~isfield(S,'nominal') || ~isfield(S.nominal,'meteorology')
    error('Derived product lacks retained meteorology.')
end
if sum(R.Ntime)~=numel(S.nominal.meteorology.ustar)
    error('Plot N and retained meteorological-row count differ.')
end
end


function value = parameter_title(alpha,gamma,lambda)
value = sprintf(['$\\alpha=%s$, $\\gamma=%s$, ' ...
    '$\\lambda=%s$ s m$^{-1}$'], ...
    format_parameter_value(alpha),format_parameter_value(gamma), ...
    format_parameter_value(lambda));
end


function value = format_parameter_value(value)
if value==0
    value = '0';
    return
end
if abs(value)>=1
    nDecimal = 1;
else
    nDecimal = max(0,ceil(-log10(abs(value))));
end
value = sprintf(['%.' num2str(nDecimal) 'f'],round(value,nDecimal));
end


function plot_envelope(ax,x,E,ib,color)
y1 = 100*E.q25(ib,:);
y2 = 100*E.q75(ib,:);
good = isfinite(x) & isfinite(y1) & isfinite(y2) & ...
    y1>0 & y2>0 & y2>=y1;
if nnz(good)>=2
    fill(ax,[x(good),fliplr(x(good))], ...
        [y1(good),fliplr(y2(good))],color, ...
        'FaceAlpha',0.20,'EdgeColor','none','HandleVisibility','off')
end
y = 100*E.med(ib,:);
good = isfinite(x) & isfinite(y) & y>0;
plot(ax,x(good),y(good),'-','Color',color,'LineWidth',1.5, ...
    'HandleVisibility','off')
end


function plot_observations(ax,Dp_um,med,q25,q75,marker,color)
x = Dp_um(:);
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


function E = zhang_met_envelope(S,R,P,rhop,ARp,wsMethod)
% Reuse the FRAGMENT vectorized evaluator: one physical-model call per
% meteorological state, then evaluate all paired alpha/gamma draws at once.
M = S.nominal.meteorology;
u = M.ustar(:);
binID = discretize(u,R.ustarEdges);
zRef = S.referenceHeight_m.displayReference;
Dp = R.DpCurve_um(:)*1e-6;
alpha = P.alpha(:)';
gamma = P.gamma(:)';
nPair = numel(alpha);
E.med = NaN(numel(R.Ntime),numel(Dp));
E.q25 = E.med;
E.q75 = E.med;
E.Npool = zeros(numel(R.Ntime),1);
for ib = R.validPlotBins(:)'
    states = find(binID==ib);
    pool = NaN(numel(states)*nPair,numel(Dp),'single');
    for ii = 1:numel(states)
        it = states(ii);
        iz = min(it,numel(zRef));
        out = deposition_velocity(u(it),M.T(it),M.p(it),M.L(it), ...
            zRef(iz),Dp,'vdMethod','MCzhang01','wsMethod',wsMethod, ...
            'z0',M.z0(it),'rhop',rhop,'ARp',ARp,'alpha',1,'gamma',1);
        ws = out.ws(:);
        St = ws*u(it)^2/(9.81*out.air.nu);
        Sc = 1./out.Eb(:);
        R1 = 1./(3*u(it)*out.rb(:).*(out.Eb(:)+(St./(1+St)).^2));
        Eb = exp(-log(Sc)*gamma);
        Eim = (St./(alpha+St)).^2;
        rb = 1./(3*u(it)*(Eb+Eim).*R1);
        wd = ws./(1-exp(-ws.*(out.ra(:)+rb)));
        rows = (ii-1)*nPair+(1:nPair);
        pool(rows,:) = single(wd');
    end
    q = prctile(pool,[25 50 75],1);
    E.q25(ib,:) = double(q(1,:));
    E.med(ib,:) = double(q(2,:));
    E.q75(ib,:) = double(q(3,:));
    E.Npool(ib) = size(pool,1);
    fprintf('Anza bin %d: %d states x %d paired draws\n', ...
        ib,numel(states),nPair)
end
end
