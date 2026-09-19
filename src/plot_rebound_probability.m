% Plot the two rebound/sticking-probability parameterizations in Figure 2.

thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
projectRoot = fileparts(thisDir);
outputDir = fullfile(projectRoot,'figures','generated');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
addpath(fullfile(projectRoot,'src'))

T = 299;          % K
p = 1013e2;       % Pa
ustar = 0.2;      % m s^-1
g = 9.81;         % m s^-2
Rd = 287;         % J kg^-1 K^-1
S = 120;          % K
mu0 = 1.812e-5;   % kg m^-1 s^-1
T0 = 288.15;      % K

Dp_um = [0.1:0.01:1 1.1:0.1:10 11:100]';
Dp_m = Dp_um*1e-6;
rhoa = p/(Rd*T);
mu = mu0*(T0+S)/(T+S)*(T/T0)^(3/2);
nu = mu/rhoa;
ws = settling_velocity(T,p,Dp_m,'stokes');
St = ws*ustar^2/(g*nu);

R1Slinn = exp(-sqrt(St));
% Rebound coefficient used in the manuscript (see PleimPrebA in
% deposition_velocity.m).
R1Saltation = 1 - 0.95*(1-exp(-ws));

fig = figure('Color','w','Position',[100 100 560 420]);
ax = subplot(2,2,1,'Parent',fig);
plot(ax,Dp_um,R1Slinn,Dp_um,R1Saltation,'LineWidth',1.2)
set(ax,'TickLabelInterpreter','latex','FontSize',12,'XScale','log')
xlabel(ax,'$D_p$ ($\mu$m)','Interpreter','latex')
ylabel(ax,'$R_1$','Interpreter','latex')
legend(ax,{'Slinn (1982)','Saltation-consistent'}, ...
    'Interpreter','latex','Location','southwest','EdgeColor','none')
box(ax,'on')

savefig(fig,fullfile(outputDir,'R1.fig'))
exportgraphics(ax,fullfile(outputDir,'R1.jpg'),'Resolution',600)
exportgraphics(ax,fullfile(outputDir,'R1.pdf'),'ContentType','vector')
close(fig)
