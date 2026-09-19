% WDS_PARAMETER_UNCERTAINTY_IQR
% Uncertainty versions of manuscript Figures 5, 7, 8, and 9 with pointwise IQR
% ribbons propagated from the paired alpha/gamma Monte Carlo fits.
%
% The settling-only first figure is intentionally omitted because alpha and
% gamma do not enter settling_velocity. The original wds.m and its figures
% are not modified.

thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
projectRoot = fileparts(thisDir);
srcDir = fullfile(projectRoot,'src');
resultsDir = fullfile(projectRoot,'data','zhang2014','results');
outputDir = fullfile(projectRoot,'figures','generated');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
addpath(srcDir)

% One-vote-per-diameter parameter ensembles. Both fits used the same
% measurement-error realizations (rngSeed=1), so matching row indices are
% retained when ratios involve both deposition formulations.
mcFile = fullfile(resultsDir,'zhang14_fit_default_cliftgauvin_n2000.mat');
eaFile = fullfile(resultsDir, ...
    'zhang14_fit_default_cliftgauvin_zhang01_n2000.mat');

MC = load_parameter_draws(mcFile,'MCzhang01');
EA = load_parameter_draws(eaFile,'zhang2001');
nDraw = min(height(MC),height(EA));
MC = MC(1:nDraw,:);
EA = EA(1:nDraw,:);

fprintf('Using %d paired parameter realizations.\n',nDraw)
print_parameter_summary('Mass-conserving fit',MC)
print_parameter_summary('Electric-analogue fit',EA)

% Shared settings from wds.m.
T = 299;
p = 1013e2;
L = 0;
zr = 10;
z0 = 1e-3;
ustar0 = 0.2;
Dp_um = [0.1:0.01:1 1.1:0.1:10 11:100]';
Dp_m = Dp_um*1e-6;
ustarVec = 0.1:0.1:0.4;
methodColors = lines(4);
ustarColors = lines(4);
bandAlpha = 0.14;
figureVisibility = 'on';

validate_parameter_velocity(T,p,z0,MC,EA)

%% Figure 2: settling-velocity parameterization effects on mass-conserving wd
wsMethods = {'stokes','cliftgauvin','stokes_asph','bb'};
methodLabels = {'Stokes','Clift \& Gauvin','Huang 2020', ...
    'Triaxial ellipsoid'};
WD = cell(1,4);

for im = 1:4
    WD{im} = parameter_velocity(ustar0,T,p,L,zr,Dp_m, ...
        'MCzhang01',wsMethods{im},z0,MC.alpha,MC.gamma);
end

fig = figure('Visible',figureVisibility,'Color','w'); clf(fig)
set(fig,'Position',[100 100 1200 500])

ax = subplot(1,2,1,'Parent',fig); hold(ax,'on')
for im = 1:4
    plot_iqr_y(ax,Dp_um,WD{im}*100,methodColors(im,:), ...
        bandAlpha,methodLabels{im});
end
set(ax,'TickLabelInterpreter','latex','FontSize',12, ...
    'XLim',[0 100],'YLim',[0 100],'XTick',0:20:100)
xlabel(ax,'$D_p$ ($\mu$m)','Interpreter','latex')
ylabel(ax,'$w_d$ (cm s$^{-1}$)','Interpreter','latex')
legend(ax,'Interpreter','latex','Location','northwest','EdgeColor','none')
text(ax,-0.18,1.02,'a.','Units','normalized','Interpreter','latex', ...
    'FontSize',12)
box(ax,'on')

ax = subplot(1,2,2,'Parent',fig); hold(ax,'on')
for im = 1:4
    ratio = WD{1}./WD{im};
    plot_iqr_y(ax,Dp_um,ratio,methodColors(im,:),bandAlpha,methodLabels{im});
end
set(ax,'TickLabelInterpreter','latex','FontSize',12, ...
    'XLim',[0 100],'YLim',[1 1.8],'XTick',0:20:100)
xlabel(ax,'$D_p$ ($\mu$m)','Interpreter','latex')
ylabel(ax,'$w_{d,\mathrm{Stokes}}/w_d$','Interpreter','latex')
text(ax,-0.18,1.02,'b.','Units','normalized','Interpreter','latex', ...
    'FontSize',12)
box(ax,'on')
save_figure_set(fig,outputDir,'wdws_parameter_IQR')
close(fig)

%% Figure 3: electric-analogue and mass-conserving comparisons versus Dp
WDmc = cell(1,4);
WDea = cell(1,4);
WDmcEaParameters = cell(1,4);
WDr1 = cell(1,4);

for iu = 1:4
    u = ustarVec(iu);
    WDmc{iu} = parameter_velocity(u,T,p,0,zr,Dp_m, ...
        'MCzhang01','bb',z0,MC.alpha,MC.gamma);
    WDea{iu} = parameter_velocity(u,T,p,0,zr,Dp_m, ...
        'zhang01','bb',z0,EA.alpha,EA.gamma);
    WDmcEaParameters{iu} = parameter_velocity(u,T,p,0,zr,Dp_m, ...
        'MCzhang01','bb',z0,EA.alpha,EA.gamma);
    % Use each mass-conserving fit draw in both the MC and R1-saltation
    % formulations so panel d isolates the formulation change.
    WDr1{iu} = parameter_velocity(u,T,p,0,zr,Dp_m, ...
        'zhang01r1','bb',z0,MC.alpha,MC.gamma);
end

% Match the geometry and styling saved in the original wdewd0.fig.
fig = figure('Visible',figureVisibility,'Color',[0.94 0.94 0.94], ...
    'Position',[100 100 560 420]); clf(fig)

ax = subplot(2,2,1,'Parent',fig); hold(ax,'on')
plot_iqr_y_original_style(ax,Dp_um,WDmc{3}*100,methodColors(1,:), ...
    bandAlpha,'$w_d$');
plot_iqr_y_original_style(ax,Dp_um,WDea{3}*100,methodColors(2,:), ...
    bandAlpha,'$w_d^e$');
wsBB = settling_velocity(T,p,Dp_m,'bb');
plot(ax,Dp_um,wsBB*100,':k','LineWidth',0.5,'DisplayName','$w_s$')
set(ax,'TickLabelInterpreter','latex','FontSize',12,'XScale','log', ...
    'YScale','log','YLim',[1e-1 1e2],'YTick',[1e-1 1 10 100])
xlabel(ax,'$D_p$ ($\mu$m)','Interpreter','latex')
ylabel(ax,'Velocity (cm s$^{-1}$)','Interpreter','latex')
legend(ax,'Interpreter','latex','Location','northwest','EdgeColor','none', ...
    'FontSize',12)
title(ax,'$u_*=0.3$','Interpreter','latex','FontSize',12)
text(ax,2.5e-2,1e2,'a.','Interpreter','latex','FontSize',12)
box(ax,'on')

ax = subplot(2,2,2,'Parent',fig); hold(ax,'on')
for iu = 1:4
    plot_iqr_y_original_style(ax,Dp_um,WDea{iu}./WDmc{iu},ustarColors(iu,:), ...
        bandAlpha,sprintf('$u_*=%.1f$',ustarVec(iu)));
end
plot(ax,[0.1 100],[1 1],'-k','LineWidth',0.1,'HandleVisibility','off')
set(ax,'TickLabelInterpreter','latex','FontSize',12,'XScale','log', ...
    'YLim',[0.6 1.4],'YTick',0.6:0.2:1.4)
xlabel(ax,'$D_p$ ($\mu$m)','Interpreter','latex')
ylabel(ax,'$w_d^e/w_d$','Interpreter','latex')
legend(ax,'Interpreter','latex','Location','southwest','EdgeColor','none', ...
    'FontSize',12,'NumColumns',2)
text(ax,2.5e-2,1.4,'b.','Interpreter','latex','FontSize',12)
box(ax,'on')

ax = subplot(2,2,3,'Parent',fig); hold(ax,'on')
for iu = 1:4
    plot_iqr_y_original_style(ax,Dp_um,WDmcEaParameters{iu}./WDmc{iu}, ...
        ustarColors(iu,:),bandAlpha,sprintf('$u_*=%.1f$',ustarVec(iu)));
end
yline(ax,1,':k','HandleVisibility','off')
set(ax,'TickLabelInterpreter','latex','FontSize',12,'XScale','log', ...
    'YLim',[1 1.6],'YTick',0.8:0.2:1.6)
xlabel(ax,'$D_p$ ($\mu$m)','Interpreter','latex')
ylabel(ax,'$w_d(\alpha_e,\gamma_e)/w_d(\alpha_m,\gamma_m)$', ...
    'Interpreter','latex')
text(ax,2.5e-2,1.6,'c.','Interpreter','latex','FontSize',12)
box(ax,'on')

ax = subplot(2,2,4,'Parent',fig); hold(ax,'on')
for iu = 1:4
    plot_iqr_y_original_style(ax,Dp_um,WDr1{iu}./WDmc{iu},ustarColors(iu,:), ...
        bandAlpha,sprintf('$u_*=%.1f$',ustarVec(iu)));
end
set(ax,'TickLabelInterpreter','latex','FontSize',12,'XScale','log', ...
    'YLim',[1 1.3],'YTick',0.8:0.1:1.4)
xlabel(ax,'$D_p$ ($\mu$m)','Interpreter','latex')
ylabel(ax,'$w_d^e(\alpha_m,\gamma_m)/w_d(\alpha_m,\gamma_m)$', ...
    'Interpreter','latex')
text(ax,2.5e-2,1.3,'d.','Interpreter','latex','FontSize',12)
text(ax,2e-1,1.2,{ ...
    'Saltation $R_1$'}, ...
    'Interpreter','latex','FontSize',11)
box(ax,'on')


save_wdewd0_figure_set(fig,outputDir,'wdewd0_parameter_IQR')
close(fig)

%% Figure 4: electric-analogue / mass-conserving ratio versus height
% DpHeight_um = [0.1 1 2.5 10 62.5];
DpHeight_um = [2.5 10 20 62.5];

DpHeight_m = DpHeight_um*1e-6;
z = [0.1:0.01:1 1.1:0.1:10 11:50];
Lvec = [0 10];
ratioHeight = cell(1,2);

for iL = 1:2
    ratioHeight{iL} = nan(nDraw,numel(z),numel(DpHeight_m));
    for iz = 1:numel(z)
        Vmc = parameter_velocity(0.2,T,p,Lvec(iL),z(iz),DpHeight_m, ...
            'MCzhang01','bb',z0,MC.alpha,MC.gamma);
        Vea = parameter_velocity(0.2,T,p,Lvec(iL),z(iz),DpHeight_m, ...
            'zhang01r1','bb',z0,MC.alpha,MC.gamma);
        ratioHeight{iL}(:,iz,:) = Vea./Vmc;
    end
end

fig = figure('Visible',figureVisibility,'Color','w'); clf(fig)
% set(fig,'Position',[100 100 1250 650])
heightColors = lines(numel(DpHeight_m));
heightLabels = arrayfun(@(d)sprintf('$D_p=%.3g$ $\\mu$m',d), ...
    DpHeight_um,'UniformOutput',false);

for iL = 1:2
    ax = subplot(1,2,iL,'Parent',fig); hold(ax,'on')
    for id = 1:numel(DpHeight_m)
        Y = squeeze(ratioHeight{iL}(:,:,id));
        plot_iqr_x(ax,Y,z,heightColors(id,:),bandAlpha,heightLabels{id});
    end
    set(ax,'TickLabelInterpreter','latex','FontSize',12,'XLim',[1 1.3], ...
        'YLim',[0 50])
    xlabel(ax,'$w_d^e/w_d$','Interpreter','latex')
    ylabel(ax,'$z_r$ (m)','Interpreter','latex')
    title(ax,sprintf('$L=%g$, $u_*=0.3$',Lvec(iL)), ...
        'Interpreter','latex')
    text(ax,-0.18,1.02,sprintf('%c.','a'+iL-1),'Units','normalized', ...
        'Interpreter','latex','FontSize',12)
    box(ax,'on')
end
legend(ax,'Interpreter','latex','Location','southwest','EdgeColor','none')
save_figure_set(fig,outputDir,'wdewdz_parameter_IQR')
close(fig)


%% Figure 5: combined electric-analogue/Stokes versus MC/triaxial effects
ratioTotal = cell(1,4);
for iu = 1:4
    u = ustarVec(iu);
    Vmc = parameter_velocity(u,T,p,0,10,Dp_m, ...
        'MCzhang01','bb',z0,MC.alpha,MC.gamma);
    VeaStokes = parameter_velocity(u,T,p,0,10,Dp_m, ...
        'zhang01','stokes',z0,EA.alpha,EA.gamma);
    ratioTotal{iu} = VeaStokes./Vmc;
end

% Match the geometry and styling saved in the original TotalDiff.fig.
fig = figure('Visible',figureVisibility,'Color',[0.94 0.94 0.94], ...
    'Position',[100 100 560 420]); clf(fig)
ax = subplot(2,2,1,'Parent',fig); hold(ax,'on')
for iu = 1:4
    plot_raw_median_iqr(ax,Dp_um,ratioTotal{iu}, ...
        ustarColors(iu,:),bandAlpha, ...
        sprintf('$u_*=%.1f$',ustarVec(iu)));
end
yline(ax,1,':k','HandleVisibility','off')
set(ax,'TickLabelInterpreter','latex','FontSize',12,'XScale','log', ...
    'XLim',[0.1 100],'YLim',[0.8 1.7])
legend(ax,'Interpreter','latex','Location','northwest','EdgeColor','none')
ylabel(ax,'$w_{d,\mathrm{Stokes}}^e/w_{d,\mathrm{triaxial}}$', ...
    'Interpreter','latex')
xlabel(ax,'$D_p$ ($\mu$m)','Interpreter','latex')
box(ax,'on')
save_totaldiff_figure_set(fig,ax,outputDir,'TotalDiff_parameter_IQR')


close(fig)
fprintf('Saved uncertainty-aware figures to %s\n',outputDir)

%% Local functions
function draws = load_parameter_draws(fileName,expectedMethod)
if ~isfile(fileName)
    error('Required fit file not found: %s',fileName)
end
S = load(fileName,'results');
method = validatestring(S.results.vdMethod, ...
    {'MCzhang01','zhang01','zhang2001'});
if strcmpi(expectedMethod,'zhang2001')
    methodOK = any(strcmpi(method,{'zhang01','zhang2001'}));
else
    methodOK = strcmpi(method,expectedMethod);
end
if ~methodOK
    error('Fit %s uses vdMethod=%s, expected %s.', ...
        fileName,method,expectedMethod)
end
fit = S.results.nonsticky;
draws = fit.mc(fit.accepted,{'alpha','gamma'});
if isempty(draws)
    error('No accepted parameter draws in %s.',fileName)
end
end

function print_parameter_summary(label,draws)
fprintf(['%s: alpha %.4g [IQR %.4g, %.4g], ' ...
    'gamma %.4g [IQR %.4g, %.4g]\n'],label, ...
    median(draws.alpha),prctile(draws.alpha,25),prctile(draws.alpha,75), ...
    median(draws.gamma),prctile(draws.gamma,25),prctile(draws.gamma,75))
end

function V = parameter_velocity(ustar,T,p,L,zr,Dp,vdMethod,wsMethod,z0, ...
    alpha,gamma)
% Evaluate all alpha/gamma draws efficiently while reproducing the relevant
% branches of deposition_velocity.m. Meteorology and particle properties
% are fixed within each call, so ra, ws, Sc, and St need only be computed
% once.
base = deposition_velocity(ustar,T,p,L,zr,Dp, ...
    'vdMethod',vdMethod,'wsMethod',wsMethod,'z0',z0, ...
    'alpha',median(alpha),'gamma',median(gamma));

Dp = Dp(:);
ws = base.ws(:).';
ra = base.ra(:).';
nu = base.air.nu;
mu = base.air.mu;
lambdaAir = base.air.lambda;

Cc = 1 + (2*lambdaAir./Dp).*(1.257 + ...
    0.40*exp(-(1.10*Dp)./(2*lambdaAir)));
Db = 1.38065e-23*T.*Cc./(3*pi*mu.*Dp);
Sc = (nu./Db).';
St = (base.ws(:)*ustar^2/(9.81*nu)).';

alpha = alpha(:);
gamma = gamma(:);
Eb = exp(-gamma*log(Sc));
EIM = (St./(alpha + St)).^2;

switch lower(vdMethod)
    case {'mczhang01','mczhang','mczhang01_bare'}
        R1 = max(1 - 0.95*(1-exp(-ws)),1e-6);
        rb = 1./(3*ustar*(Eb+EIM).*R1);
        V = ws./(1-exp(-ws.*(ra+rb)));
    case {'zhang01','zhang2001'}
        R1 = exp(-sqrt(St));
        rb = 1./(3*ustar*(Eb+EIM).*R1);
        V = ws + 1./(ra+rb+ra.*rb.*ws);
    case {'zhang01r1','zhang2001r1'}
        R1 = max(1 - 0.95*(1-exp(-ws)),1e-6);
        rb = 1./(3*ustar*(Eb+EIM).*R1);
        V = ws + 1./(ra+rb+ra.*rb.*ws);
    otherwise
        error('Unsupported fitted deposition method: %s',vdMethod)
end
end

function plot_iqr_y(ax,x,Y,color,faceAlpha,label)
q = prctile(Y,[25 50 75],1);
fill(ax,[x(:);flipud(x(:))], ...
    [q(1,:)';flipud(q(3,:)')],color, ...
    'FaceAlpha',faceAlpha,'EdgeColor','none','HandleVisibility','off');
plot(ax,x,q(2,:)','Color',color,'LineWidth',1.3,'DisplayName',label)
end

function plot_iqr_y_original_style(ax,x,Y,color,faceAlpha,label)
q = prctile(Y,[25 50 75],1);
fill(ax,[x(:);flipud(x(:))], ...
    [q(1,:)';flipud(q(3,:)')],color, ...
    'FaceAlpha',faceAlpha,'EdgeColor','none','HandleVisibility','off');
plot(ax,x,q(2,:)','Color',color,'LineWidth',0.5,'DisplayName',label)
end

function plot_raw_median_iqr(ax,x,Y,color,faceAlpha,label)
q = prctile(Y,[25 50 75],1);
fill(ax,[x(:);flipud(x(:))], ...
    [q(1,:)';flipud(q(3,:)')],color, ...
    'FaceAlpha',faceAlpha,'EdgeColor','none','HandleVisibility','off');
plot(ax,x,q(2,:)','Color',color,'LineWidth',0.5,'DisplayName',label)
end

function plot_iqr_x(ax,X,y,color,faceAlpha,label)
q = prctile(X,[25 50 75],1);
fill(ax,[q(1,:)';flipud(q(3,:)')], ...
    [y(:);flipud(y(:))],color, ...
    'FaceAlpha',faceAlpha,'EdgeColor','none','HandleVisibility','off');
plot(ax,q(2,:)',y,'Color',color,'LineWidth',1.3,'DisplayName',label)
end

function save_figure_set(fig,outputDir,baseName)
savefig(fig,fullfile(outputDir,[baseName '.fig']))
exportgraphics(fig,fullfile(outputDir,[baseName '.png']),'Resolution',300)
exportgraphics(fig,fullfile(outputDir,[baseName '.jpg']),'Resolution',600)
exportgraphics(fig,fullfile(outputDir,[baseName '.pdf']),'ContentType','vector')
end

function save_wdewd0_figure_set(fig,outputDir,baseName)
savefig(fig,fullfile(outputDir,[baseName '.fig']))
exportgraphics(fig,fullfile(outputDir,[baseName '.png']),'Resolution',600)
exportgraphics(fig,fullfile(outputDir,[baseName '.jpg']),'Resolution',600)
exportgraphics(fig,fullfile(outputDir,[baseName '.pdf']),'ContentType','vector')
end

function save_totaldiff_figure_set(fig,ax,outputDir,baseName)
% Save the editable full MATLAB figure, but crop publication exports to the
% axes as was done for the original 3.36-by-2.68-inch TotalDiff output.
savefig(fig,fullfile(outputDir,[baseName '.fig']))
exportgraphics(ax,fullfile(outputDir,[baseName '.png']),'Resolution',600)
exportgraphics(ax,fullfile(outputDir,[baseName '.jpg']),'Resolution',600)
exportgraphics(ax,fullfile(outputDir,[baseName '.pdf']),'ContentType','vector')
end

function validate_parameter_velocity(T,p,z0,MC,EA)
% Check the vectorized evaluator against the project implementation for a
% few representative paired draws before generating any figures.
Dp = [0.2 2 20]'*1e-6;
indices = unique([1 round(height(MC)/2) height(MC)]);
cases = { ...
    'MCzhang01','bb',MC; ...
    'zhang01','bb',EA; ...
    'zhang01r1','bb',MC};

for ic = 1:size(cases,1)
    method = cases{ic,1};
    settling = cases{ic,2};
    draws = cases{ic,3};
    fast = parameter_velocity(0.3,T,p,0,10,Dp,method,settling,z0, ...
        draws.alpha(indices),draws.gamma(indices));
    direct = nan(size(fast));
    for ii = 1:numel(indices)
        out = deposition_velocity(0.3,T,p,0,10,Dp, ...
            'vdMethod',method,'wsMethod',settling,'z0',z0, ...
            'alpha',draws.alpha(indices(ii)), ...
            'gamma',draws.gamma(indices(ii)));
        direct(ii,:) = out.wd(:)';
    end
    relativeError = max(abs(fast-direct)./max(abs(direct),realmin),[],'all');
    assert(relativeError<1e-10, ...
        'Vectorized %s evaluator differs from deposition_velocity (%.3g).', ...
        method,relativeError)
end
fprintf('Vectorized parameter propagation validated against deposition_velocity.\n')
end
