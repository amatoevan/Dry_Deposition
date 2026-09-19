% Generate all of the plots needed to quantify the effects of different
% representations of the settling velocity and the depostion velocity,
% relative to the prolate spheroid and mass conserving case
%
% For the prolate spheroid case, the diameter isn't the volume equivalent
% diameter. So this has to be adjusted. The volume of the prolate spheroid
% is
%
% Vp ​= pi*4/3 * (Dp/2)^3 * 1/AR^2
%
% Vs = pi*4/3 * (Ds/2)^3
%
% So the equivalent volume conserving diameter: Ds = Dp * AR^(-2/3)
thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
projectRoot = fileparts(thisDir);
addpath(fullfile(projectRoot,'src'))
outputDir = fullfile(projectRoot,'figures','generated');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

%% ----- Settings
L = 0;      % Obukhov Length [m]
zr = 10;     % reference height

% Constants
T = 299;        % [K]
p = 1013e2;     % [Pa]
ustar = 0.2;    % [m/s]
z0 = 1e-3;      % [m]

% Particle diameters
Dp_um  = [0.1:0.01:1 1.1:.1:10 11:100]';   % [um]
Dp_m = Dp_um * 1e-6;                         % [m]
nDp = numel(Dp_m);

% Height grid
z = [0.1:0.01:1 1.1:0.1:10 11:50];  % [m]
z = z(:)';                           % row vector
[~, iz1] = min(abs(z - 1));          % robust z_r=1 index
nz = numel(z);

% From the 1-vote-per-diameter fitting:
%   alpha = 3.67196
%   gamma = 0.351094

%

%% ---- Plot ws for the three cases -- ws.fig
ws = NaN(numel(Dp_m),3); % rhop default is 2,500 kg/m3
ws(:,1) = settling_velocity(T, p, Dp_m, 'stokes');
ws(:,2) = settling_velocity(T, p, Dp_m, 'cliftgauvin');
ws(:,3) = settling_velocity(T, p, Dp_m, 'stokes_asph');
ws(:,4) = settling_velocity(T, p, Dp_m, 'bb');

figure(1); clf

subplot(2,2,1)
plot(Dp_um,ws*100)
set(gca,'ticklabelinterpreter','latex','fontsize',12,'xtick',0:20:100, ...
    'ylim',[0 80])
xlabel('$D_p$ ($\mu$m)','interpreter','latex')
ylabel('$w_s$ (cm s$^{-1}$)','interpreter','latex')
legend({'Stokes','Clift \& Gauvin','Huang 2020','Triaxial ellipsoid'}, ...
    'interpreter','latex','location','northwest','edgecolor','none')
text(-20,80,'a.','interpreter','latex','fontsize',12)


% subplot(2,2,1)
% plot(Dp_um,ws*100)
% set(gca,'ticklabelinterpreter','latex','fontsize',12,'xtick',0:2:10, ...
%     'xlim',[0 10])
% xlabel('$D_p$ ($\mu$m)','interpreter','latex')
% ylabel('$w_s$ (cm s$^{-1}$)','interpreter','latex')
% legend({'Stokes','Clift \& Gauvin','Triaxial ellipsoid','Huang 2020'}, ...
%     'interpreter','latex','location','northwest','edgecolor','none')
% text(-2,0.8,'a.','interpreter','latex','fontsize',12)

subplot(2,2,2)
y = ws;
y(:,1) = 0;
%plot(Dp_um,y./ws(:,1))
plot(Dp_um,ws(:,1)./y)
set(gca,'ticklabelinterpreter','latex','fontsize',12,'xtick',0:20:100, ...
    'ylim',[1 1.8])
text(-20,1.8,'b.','interpreter','latex','fontsize',12)
xlabel('$D_p$ ($\mu$m)','interpreter','latex')
ylabel('$w_{s,\mathrm{Stokes}}/w_s$','interpreter','latex')

save_figure_set(gcf,outputDir,'ws')
close



%% ---- Plot wd for the three cases - wdws.fig
wd = NaN(numel(Dp_m),4); % rhop default is 2,500 kg/m3

out = deposition_velocity(ustar, T, p, L, zr, Dp_m, 'vdMethod','MCzhang01', ...
    'wsMethod','stokes');
wd(:,1) = out.wd;

out = deposition_velocity(ustar, T, p, L, zr, Dp_m, 'vdMethod','MCzhang01', ...
    'wsMethod','cliftgauvin');
wd(:,2) = out.wd;

out = deposition_velocity(ustar, T, p, L, zr, Dp_m, 'vdMethod','MCzhang01', ...
    'wsMethod','stokes_asph');
wd(:,3) = out.wd;

out = deposition_velocity(ustar, T, p, L, zr, Dp_m, 'vdMethod','MCzhang01', ...
    'wsMethod','bb');
wd(:,4) = out.wd;

figure(2); clf

subplot(2,2,1)
plot(Dp_um,wd*100)
set(gca,'ticklabelinterpreter','latex','fontsize',12, ...
    'xscale','linear','yscale','linear','ylim',[0 100],'xtick',0:20:100)
% hold on
% plot(Dp_um,ws(:,3)*100,'--k')
% hold off
xlabel('$D_p$ ($\mu$m)','interpreter','latex')
ylabel('$w_d$ (cm s$^{-1}$)','interpreter','latex')
legend({'Stokes','Clift \& Gauvin','Huang 2020','Triaxial Ellipsoid'}, ...
    ... %'$w_{s,\mathrm{Ellipsoid}}$'}, ...
    'interpreter','latex','location','northwest','edgecolor','none')
text(-20,100,'a.','interpreter','latex','fontsize',12)

subplot(2,2,2)
%plot(Dp_um,wd./wd(:,1))
plot(Dp_um,wd(:,1)./wd)
set(gca,'ticklabelinterpreter','latex','fontsize',12,'xscale','linear', ...
    'YLim',[1 1.8],'xtick',0:20:100)
xlabel('$D_p$ ($\mu$m)','interpreter','latex')
ylabel('$w_d$ ratio','interpreter','latex')
ylabel('$w_{d,\mathrm{Stokes}}/w_d$','interpreter','latex')

% legend({'','Clift \& Gauvin','Huang 2020','Triaxial Ellipsoid'}, ...
%     'interpreter','latex','location','northwest','edgecolor','none')
text(-20,1.8,'b.','interpreter','latex','fontsize',12)


save_figure_set(gcf,outputDir,'wdws')
close


%% --- plot wde and wdo as a function of Dp - wdewd0.fig
figure(3); clf

% Lab fit: alpha = 3.3, gamma = 0.35
alphaWde = 1.4;
gammaWde = 0.34;
methe = 'zhang01'; %   'zhang01r1'; %  'MCzhang01'; % %
meths = 'bb';

ustar = 0.1:0.1:0.4;
wd = NaN(length(Dp_m),length(ustar));
wde = wd;
ra = wd;
rb = wd;
ws_ = wd;
for i = 1:numel(ustar)
    out = deposition_velocity(ustar(i), T, p, L, zr, Dp_m, ...
        'vdMethod','MCzhang01','wsMethod','bb');
    wd(:,i) = out.wd;

    out = deposition_velocity(ustar(i), T, p, L, zr, Dp_m,'vdMethod', ...
        methe,'wsMethod',meths,'alpha',alphaWde,'gamma',gammaWde);
    wde(:,i) = out.wd;

    ra(:,i) = out.ra;
    rb(:,i) = out.rb;
    ws_(:,i) = out.ws;
end

subplot(2,2,1)
plot(Dp_um,wd(:,3)*100,Dp_um,wde(:,3)*100)
hold on
plot(Dp_um,ws(:,4)*100,':k')
hold off
set(gca,'ticklabelinterpreter','latex','fontsize',12,'xscale','log', ...
    'yscale','log','ylim',[1e-1,1e2],'ytick',[1e-1 1e0 1e1 1e2])
xlabel('$D_p$ ($\mu$m)','interpreter','latex')
ylabel('Velocity (cm s$^{-1}$)','interpreter','latex')
legend({'$w_d$','$w_d^e$','$w_s$'},'fontsize',12, ...
    'interpreter','latex','location','northwest','edgecolor','none')
text(2.5e-2,1e2,'a.','interpreter','latex','fontsize',12)
title('$u_*=0.3$','Interpreter','Latex','fontsize',12);

subplot(2,2,2)
for i = 1:4
    plot(Dp_um,wde(:,i)./wd(:,i))
    if i==1; hold on; end
end
hold off
set(gca,'ticklabelinterpreter','latex','fontsize',12,'xscale','log', ...
    'YLim',[0.6 1.4],'ytick',0.6:.2:1.4)
hold on; plot([.1 100],[1 1],'-k','linewidth',0.1); hold off
xlabel('$D_p$ ($\mu$m)','interpreter','latex')
ylabel('$w_d^e/w_d$','interpreter','latex')
text(2.5e-2,1.3,'b.','interpreter','latex','fontsize',12)

legend({'$u_*=0.1$','$u_*=0.2$','$u_*=0.3$','$u_*=0.4$'}, ...
    'Interpreter','Latex','edgecolor','none','location','southwest', ...
    'fontsize',12)


% Panel c
methe = 'MCzhang01'; % 'zhang01'; %   'zhang01r1'; %  'MCzhang01'; % %
meths = 'bb';
alphaWde = 1.4;
gammaWde = 0.35;

ustar = 0.1:0.1:0.4;
wd  = NaN(length(Dp_m),length(ustar));
wde = wd;
Eim  = wd;
Eime = wd;
for i = 1:numel(ustar)
    out = deposition_velocity(ustar(i), T, p, L, zr, Dp_m, ...
        'vdMethod','MCzhang01','wsMethod','bb');
    wd(:,i) = out.wd;
    Eim(:,i) = out.EIM;
    out = deposition_velocity(ustar(i), T, p, L, zr, Dp_m,'vdMethod', ...
        methe,'wsMethod',meths,'alpha',alphaWde,'gamma',gammaWde);
    wde(:,i) = out.wd;
    Eime(:,i) = out.EIM;
end

subplot(2,2,3)
for i = 1:4
    plot(Dp_um,wde(:,i)./wd(:,i))
    if i==1; hold on; end
end
hold off
set(gca,'ticklabelinterpreter','latex','fontsize',12,'xscale','log', ...
    'YLim',[1 1.6],'ytick',0.8:.2:1.6)
hold on; plot([.1 100],[1 1],':k'); hold off
xlabel('$D_p$ ($\mu$m)','interpreter','latex')
ylabel('$w_{d,\,\alpha=1.4,\,\gamma=0.34}/w_d$','interpreter','latex')
text(2.5e-2,1.6,'c.','interpreter','latex','fontsize',12)


% Panel d.
alphaWde = 3.3; %
gammaWde = 0.35; %
methe = 'zhang01r1'; % 'zhang01'; %   'zhang01r1'; %  'MCzhang01'; % %
meths = 'bb';

ustar = 0.1:0.1:0.4;
wd = NaN(length(Dp_m),length(ustar));
wde = wd;
for i = 1:numel(ustar)
    out = deposition_velocity(ustar(i), T, p, L, zr, Dp_m, ...
        'vdMethod','MCzhang01','wsMethod','bb');
    wd(:,i) = out.wd;
    out = deposition_velocity(ustar(i), T, p, L, zr, Dp_m,'vdMethod', ...
        methe,'wsMethod',meths,'alpha',alphaWde,'gamma',gammaWde);
    wde(:,i) = out.wd;
end

subplot(2,2,4)
for i = 1:4
    plot(Dp_um,wde(:,i)./wd(:,i))
    if i==1; hold on; end
end
hold off
set(gca,'ticklabelinterpreter','latex','fontsize',12,'xscale','log', ...
    'YLim',[1 1.3],'ytick',0.8:.1:1.4)
% hold on; plot([.1 100],[1 1],':k'); hold off
xlabel('$D_p$ ($\mu$m)','interpreter','latex')
ylabel('$w_d^e/w_d$','interpreter','latex')
text(2.5e-2,1.3,'d.','interpreter','latex','fontsize',12)

text(2e-1,1.2,{'$w_d^e\,$:','  $\alpha = 3.3$','  $\gamma = 0.35$', ...
    '  $R_1\,\mathrm{Saltation}$'},'Interpreter','Latex','fontsize',11)


save_figure_set(gcf,outputDir,'wdewd0')
close



%% --- plot wde / wdo as a function of height @ Dp - wdewdz.fig
% Unstable case looks very much like the neutral case
% Dp_um = [1 3 10 63 100];
% Dp_um = [1 10 20 63 100];
Dp_um = [2.5 10 20 62.5];

% Dp_um = [0.3 1.6 5 10 25 80];

alphaWde = 3.3; % 1.4;
gammaWde = 0.35; % 0.34;

ustar = 0.2;
Dp_m = Dp_um * 1e-6;
L = [-10 0 10];

wde = NaN(numel(z),numel(Dp_m),numel(L));
wd0 = NaN(numel(z),numel(Dp_m),numel(L));
ra = NaN(numel(z),numel(Dp_m),numel(L));
ws = NaN(numel(Dp_m),1);
rb = NaN(numel(Dp_m),numel(L));
rbe = NaN(numel(Dp_m),numel(L));

for i = 1:numel(z)
    for j = 1:numel(L)

        out = deposition_velocity(ustar, T, p, L(j), z(i), Dp_m, ...
            'vdMethod','MCzhang01','wsMethod','bb');
        wd0(i,:,j) = out.wd;
        ra(i,:,j) = out.ra;
        if i==1; rb(:,j) = out.rb; ws(:) = out.ws; end

        out = deposition_velocity(ustar, T, p, L(j), z(i), Dp_m, ...
            'vdMethod','zhang01r1','wsMethod','bb', ...
            'alpha',alphaWde,'gamma',gammaWde);
        wde(i,:,j) = out.wd;
        if i==1; rbe(:,j) = out.rb; end

    end
end


figure(4); clf

subplot(2,2,1)
plot(wde(:,:,2)./wd0(:,:,2),z)
set(gca,'ticklabelinterpreter','latex','fontsize',12,'xlim',[1 1.3])
ylabel('$z_r$ (m)','interpreter','latex')
xlabel('$w_d^e/w_d$','interpreter','latex')
title('$L=0$, $u_*=0.3$','interpreter','latex')
text(0.94,50,'a.','interpreter','latex','fontsize',12)


% Stretch vertically
pos = get(gca,'Position');
set(gca,'Position',[pos(1) pos(2)+pos(4)-1.5*pos(4) pos(3) 1.5*pos(4)])

text(1.06,20,{'$w_d^e\,$:','  $\alpha = 3.3$','  $\gamma = 0.35$', ...
    '  $R_1\,\mathrm{Saltation}$'},'Interpreter','Latex','fontsize',11)



subplot(2,2,2)
plot(wde(:,:,3)./wd0(:,:,3),z)
set(gca,'ticklabelinterpreter','latex','fontsize',12,'xlim',[1 1.3])
xlabel('$w_d^e/w_d$','interpreter','latex')
ylabel('$z_r$ (m)','interpreter','latex')
title('$L=10$, $u_*=0.3$','interpreter','latex')
text(0.94,50,'b.','interpreter','latex','fontsize',12)


% Stretch vertically
pos = get(gca,'Position');
set(gca,'Position',[pos(1) pos(2)+pos(4)-1.5*pos(4) pos(3) 1.5*pos(4)])

% '$D_p=0.1$ $\mu$m','$D_p=1.0$ $\mu$m',
lg = legend({'$D_p=2.5$ $\mu$m', ...
    '$D_p=10$ $\mu$m','$D_p=20$ $\mu$m','$D_p=62.5$ $\mu$m'}, ...
    'interpreter','latex','location','northwest','edgecolor','none');
lg.Position = [0.1718    0.7294    0.1911    0.1686];


save_figure_set(gcf,outputDir,'wdewdz')
close


%% ---- "Combined effects" - TotalDiff.fig
figure(5); clf

Dp_um  = [0.1:0.01:1 1.1:.1:10 11:100]';   % [um]
Dp_m = Dp_um * 1e-6;                       % [m]
L = 0;

% alpha and gamma for the electric analogue simulations w/Stokes
alphaWde = 1.4; % 3.3; %
gammaWde = 0.34; % 0.35; %
methe = 'zhang01'; % 'MCzhang01'; %
meths = 'stokes'; % 'bb'; % 'stokes';

zr = 10;
ustar = 0.1:0.1:0.4;
wd = NaN(length(Dp_m),length(ustar));
wde = wd;
for i = 1:numel(ustar)
    out = deposition_velocity(ustar(i), T, p, L, zr, Dp_m, ...
        'vdMethod','MCzhang01','wsMethod','bb');
    wd(:,i) = out.wd;
    out = deposition_velocity(ustar(i), T, p, L, zr, Dp_m,'vdMethod', ...
        methe,'wsMethod',meths,'alpha',alphaWde,'gamma',gammaWde);
    wde(:,i) = out.wd;
end

subplot(2,2,1)
for i = 1:4
    plot(Dp_um,wde(:,i)./wd(:,i))
    if i==1; hold on; end
end
plot([0.1 100],[1 1],':k')
hold off

set(gca,'fontsize',12,'ticklabelinterpreter','latex','xscale','linear', ...
    'xlim',[0.1 100],'ylim',[0.8 1.7],'xscale','log')
legend({'$u_*=0.1$','$u_*=0.2$','$u_*=0.3$','$u_*=0.4$'},'location','northwest', ...
    'interpreter','latex','EdgeColor','none')
ylabel('$w_{d,\mathrm{Stokes}}^e/w_{d,\mathrm{triaxial}}$','interpreter','latex')
xlabel('$D_p$ ($\mu$m)','interpreter','latex')



save_figure_set(gcf,outputDir,'TotalDiff')
close
%%

function save_figure_set(fig,outputDir,baseName)
savefig(fig,fullfile(outputDir,[baseName '.fig']))
exportgraphics(fig,fullfile(outputDir,[baseName '.jpg']),'Resolution',600)
exportgraphics(fig,fullfile(outputDir,[baseName '.pdf']),'ContentType','vector')
end
