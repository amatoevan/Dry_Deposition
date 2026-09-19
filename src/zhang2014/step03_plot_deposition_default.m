function step03_plot_deposition_default(varargin)%(vdMethod,wsMethod,varargin)
%STEP03_PLOT_DEPOSITION_DEFAULT Plot Zhang14 observations and fitted theory.
%
% USAGE
%   step03_plot_deposition_default();
%   step03_plot_deposition_default('MCzhang01','bb');
%
% Optional Name/Value:
%   'ustarBins' default [0.10 0.20; 0.30 0.40; 0.40 0.60]
%   'NsampMC'   default 500
%   'NobsMC'    default 1000
%   'rngSeed'   default 2
%   'fitTag'     default 'default_cliftgauvin_n2000'
%   'outputBase' default 'ZhangFit'
%   'saveFigs'  default true
%   'figDir'    default 'figures/generated'

thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
projectRoot = fileparts(fileparts(fileparts(thisFile)));

addpath(fullfile(thisDir,'utils'));
addpath(fileparts(thisDir));

% if nargin<1 || isempty(vdMethod)
%     vdMethod = 'zhang2001'; % 'MCzhang01';
% end
% if nargin<2 || isempty(wsMethod)
%     wsMethod = 'cliftgauvin';
% end
%
% vdMethod = char(string(vdMethod));
% wsMethod = char(string(wsMethod));

defaultFigDir = fullfile(projectRoot,'figures','generated');

p = inputParser;
p.FunctionName = mfilename;
addParameter(p,'ustarBins',[0.10 0.20; 0.30 0.40; 0.40 0.60],@(x)isnumeric(x)&&size(x,2)==2);
addParameter(p,'NsampMC',500,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'NobsMC',1000,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'rngSeed',2,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'fitTag','default_cliftgauvin_n2000',@(s)ischar(s)||isstring(s));
addParameter(p,'outputBase','ZhangFit',@(s)ischar(s)||isstring(s));
addParameter(p,'saveFigs',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'figDir',defaultFigDir,@(s)ischar(s)||isstring(s));
parse(p,varargin{:});
opt = p.Results;

rng(opt.rngSeed,'twister');

dataFile = fullfile(projectRoot,'data','zhang2014','derived','Zhang14.mat');
if ~isfile(dataFile)
    error('Expected dataset not found: %s',dataFile);
end

S = load(dataFile,'Zhang14');
T = S.Zhang14.T;

isWood = contains(lower(string(T.surface)),'wood');
isWater = contains(lower(string(T.surface)),'water');
T = T(~(isWood|isWater),:);

tag = char(string(opt.fitTag));

resultsDir = fullfile(projectRoot,'data','zhang2014','results');
resFile = fullfile(resultsDir,['zhang14_fit_' tag '.mat']);

if ~isfile(resFile)
    error('Results file not found: %s',resFile);
end

R = load(resFile,'results');
results = R.results;

mc = results.nonsticky.mc;
accepted = results.nonsticky.accepted;
mc = mc(accepted,:);

if isempty(mc)
    error('No accepted MC parameter sets available.');
end

% Parameter medians and 95% ranges.
aMed = median(mc.alpha,'omitnan');
a95  = prctile(mc.alpha,[2.5 97.5]);
gMed = median(mc.gamma,'omitnan');
g95  = prctile(mc.gamma,[2.5 97.5]);

nTake = min(opt.NsampMC,height(mc));
iiMC = randperm(height(mc),nTake);
aMC = mc.alpha(iiMC);
gMC = mc.gamma(iiMC);

Dp_m = logspace(-1,2,20)'*1e-6;
Dp_um = Dp_m*1e6;


% modifying the bins
ustarBins = [0.1 0.2; 0.2 0.4; 0.4 0.6]; % opt.ustarBins;
unique(T.ustar_mps(~isnan(T.wd_mps)))

colors = lines(size(ustarBins,1));

figure(1);
subplot(2,2,1); hold on;
% subplot(2,2,2); hold on;

for ib = 1:size(ustarBins,1)
    u0 = ustarBins(ib,1);
    u1 = ustarBins(ib,2);

    inb = T.ustar_mps>u0 & T.ustar_mps<=u1;
    % if ib<size(ustarBins,1)
    %     inb = T.ustar_mps>=u0 & T.ustar_mps<u1;
    % else
    %     inb = T.ustar_mps>=u0 & T.ustar_mps<=u1;
    % end

    Tb = T(inb,:);
    if isempty(Tb)
        continue
    end

    [xObs,yObs,errLow,errHigh] = observation_errorbars(Tb,opt.NobsMC);
    errorbar(xObs,yObs,errLow,errHigh,'s', ...
        'MarkerSize',6,'LineWidth',1.0,'CapSize',6, ...
        'Color',colors(ib,:),'HandleVisibility','off');

    [~,ia] = unique(string(Tb.case_key),'stable');
    Tc = Tb(ia,:);
    nCase = height(Tc);

    wdPool = nan(nCase*nTake,numel(Dp_m));
    k = 0;

    for ic = 1:nCase
        for imc = 1:nTake
            k = k+1;

            out = deposition_velocity( ...
                Tc.ustar_mps(ic),298,101300,0,Tc.zm_m(ic),Dp_m, ...
                'vdMethod',results.vdMethod,'wsMethod',results.wsMethod, ...
                'z0',Tc.z0_m(ic),'rhop',2200, ...
                'alpha',aMC(imc),'gamma',gMC(imc));

            wdPool(k,:) = out.wd(:).';
        end
    end

    wdMed = median(wdPool,1,'omitnan');
    wdQ25 = prctile(wdPool,25,1);
    wdQ75 = prctile(wdPool,75,1);

    xx = [Dp_um;flipud(Dp_um)];
    yy = [wdQ25(:);flipud(wdQ75(:))]*100;
    fill(xx,yy,colors(ib,:),'FaceAlpha',0.2, ...
        'EdgeColor','none','HandleVisibility','off');

    Nbin = nnz(isfinite(Tb.wd_mps) & Tb.wd_mps>0);
    leglab = sprintf('$u_*=(%.2f-%.2f)$ m s$^{-1}$, N=%d', ...
        u0,u1,Nbin);

    plot(Dp_um,wdMed*100,'-','Color',colors(ib,:), ...
        'LineWidth',1.5,'DisplayName',leglab);
end

% Settling velocity as a dashed black line.
outWs = deposition_velocity( ...
    T.ustar_mps(1),298,101300,0,T.zm_m(1),Dp_m, ...
    'vdMethod',results.vdMethod,'wsMethod',results.wsMethod, ...
    'z0',T.z0_m(1),'rhop',2200, ...
    'alpha',aMed,'gamma',gMed);

plot(Dp_um,outWs.ws(:)*100,'k:', ...
    'LineWidth',1,'DisplayName','$w_s$');

set(gca,'XScale','log','YScale','log','FontSize',12, ...
    'TickLabelInterpreter','latex','XLim',[0.6 60], ...
    'XTick',[0.6 1 10 60],'YLim',[0.2 200]);

xlabel('$D_p$ ($\mu$m)','Interpreter','latex');
ylabel('$w_d$ (cm s$^{-1}$)','Interpreter','latex');
legend('Location','northwest','Interpreter','latex');

title({ ...
    sprintf('$\\alpha=%.2g\\;[%.2g,\\,%.3g],\\,\\gamma=%.2g\\;[%.2g,\\,%.2g]$', ...
        aMed,a95(1),a95(2),gMed,g95(1),g95(2))}, ...
    'Interpreter','latex', ...
    'FontSize',11);

% grid on;
box on;

% ax = gca;
% pp = ax.Position;
% ax.Position = [pp(1),pp(2)+pp(4)-1.5*pp(4),pp(3),1.5*pp(4)];

ax = gca;
pp = ax.Position;

stretch = 1.7;
ax.Position = [ ...
    pp(1), ...
    pp(2) - (stretch-1)*pp(4), ...
    pp(3), ...
    stretch*pp(4)];

if opt.saveFigs
    figDir = char(opt.figDir);
    if ~exist(figDir,'dir')
        mkdir(figDir);
    end
    outputBase = char(string(opt.outputBase));
    savefig(gcf,fullfile(figDir,[outputBase '.fig']));
    exportgraphics(gcf,fullfile(figDir,[outputBase '.jpg']),'Resolution',600);
    exportgraphics(gcf,fullfile(figDir,[outputBase '.pdf']),'ContentType','vector');
    fprintf('Saved %s figure set to: %s\n',outputBase,figDir);
end

end

%%
function [xObs,yObs,errLow,errHigh] = observation_errorbars(Tb,nMC)
DpVals = unique(Tb.Dp_m);
nD = numel(DpVals);

xObs = nan(nD,1);
yObs = nan(nD,1);
q25 = nan(nD,1);
q75 = nan(nD,1);

for id = 1:nD
    idx = Tb.Dp_m==DpVals(id);

    wd = Tb.wd_mps(idx)*100;
    h95 = Tb.sig_mps(idx)*100;

    good = isfinite(wd) & wd>0;
    wd = wd(good);
    h95 = h95(good);

    xObs(id) = DpVals(id)*1e6;

    if isempty(wd)
        continue
    end

    sigma1 = h95/1.96;
    sigma1(~isfinite(sigma1) | sigma1<0) = 0;
    sigmaLog = sqrt(log1p((sigma1./wd).^2));

    Wmc = exp(log(wd(:))+sigmaLog(:).*randn(numel(wd),nMC));
    medMC = median(Wmc,1,'omitnan');

    yObs(id) = median(medMC,'omitnan');
    q25(id) = prctile(medMC,25);
    q75(id) = prctile(medMC,75);
end

errLow = yObs-q25;
errHigh = q75-yObs;

ymin = 1e-3;
bad = yObs-errLow<=ymin;
errLow(bad) = yObs(bad)-ymin;
end

function s = sanitize(x)
s = regexprep(char(string(x)),'[^a-zA-Z0-9_\-]','');
end
