% plot elevation transects relevant to the forced flow:
%
% Based on the inset map describing the area (f_closeup.pdf):
%   -116.6 to -115.6 E
%	32.9 to 33.5 N
%
% location of NKX: -117.12 E, 32.85 N
% site location: -115.855930 E, 33.169230 N
clearvars
close all

% load up the colormap
% Unused terrain colormap omitted from the satellite-only figure.
thisDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(thisDir));


% 1-colum fig is 1/2 window width, 2-column fig is 1 window width
% fontsize for a 1-column (19 pc) fig is 9.5
% fontsize for a 2-column (33 pc) figure is 11.5
fs = 12.;


% load pm and wspd data
% map limits
xlim = [-116.5 -115.5];
ylim = [32.75 33.75];


% % Make a map
% subplot(2,2,3)
% imagesc(x,y,dem.z')
% set(gca,'YLim',ylim,'XLim',xlim,'ydir','normal','TickDir','out', ...
%     'TicklabelInterpreter','Latex','FontSize',fs);
%
% % shade sub sea level bluish
% demcmap(dem.z(dem.z>0))
% caxis([0 1.5])
%
% % plot location of research site
% hold on
% scatter(-115.855930,33.169230,40,'w','fill','MarkerEdgeColor','k')
% hold off
% grid on
%
% % plot location of NKX
% hold on
% scatter(-117.12,32.85,40,'w','fill','MarkerEdgeColor','k','Marker','^')
% hold off
% grid on
%
%
% % labels
% ylabel('Latitude ($^\circ$N)','Interpreter','Latex','FontSize',fs)
% xlabel('Longitude ($^\circ$E)','Interpreter','Latex','FontSize',fs)
%
% % get aspect ratio correct
% pos = get(gca,'position');
% ratio = 1.2*diff(ylim)/diff(xlim);
% pos = [pos(1) pos(2)+.12 pos(3) pos(4)*ratio];
% set(gca,'position',pos)
%
% % colorbar
% cb = colorbar;
% cb.Location = 'southoutside';
% cb.Position = cb.Position - [0 .17 0 0];
% set(cb,'TickLabelInterpreter','Latex','FontSize',10,'XTick',0:.25:1.5, ...
%     'FontSize',fs)
% xlab = get(cb,'XTickLabel');
% xlab{1} = '$\le$0';
% xlab{end} = '$\ge$1.5';
% set(cb,'XTickLabel',xlab)
%
% xlabel(cb,'Terrain height (km)','Interpreter','Latex','FontSize',fs)
%
% % shoreline
% shore = load('shoreline.mat');
% g = shore.x < 1 & shore.y < 90;
% shore.x = shore.x(g); shore.y = shore.y(g);
% hold on
% plot(shore.x,shore.y,'Color',ones(3,1)*0.5,'linewidth',1.5)
% hold off
%
% text(-117.42,34.05,'b.','Interpreter','Latex','FontSize',fs, ...
%     'BackgroundColor','w','EdgeColor','k')
%
%
% cm1 = mycm.g2brn;
% cm2 = cm1(1:20,3);
% for i = 1:3
%     cm2(:,i) = linspace(cm1(end,i),1,20);
% end
% cm = [cm1' cm2']';
% colormap(cm)
%
% % plot location of research site -- again
% hold on
% scatter(-115.855930,33.169230,40,'w','fill','MarkerEdgeColor','k')
% hold off
% grid on


%% Add satellite image

sat = 'NPP';
sat = 'NOAA-21';


% Read in the MODIS data
% [A,R] = geotiffread('VIIRS.tif');
S = load(fullfile(projectRoot,'assets','location','NOAA-21.mat'));
A = S.A;
R = S.R;


lat = linspace(R.LatitudeLimits(1),R.LatitudeLimits(2),R.RasterSize(1));
lon = linspace(R.LongitudeLimits(1),R.LongitudeLimits(2),R.RasterSize(2));
lat = fliplr(lat);

% Brighten the image and make the dust pop a bit more
A = A(:,:,1:3)*1.05;
A(:,:,2) = A(:,:,2)*1.02;
A(:,:,3) = A(:,:,3)*1.05;

% Make the figure
% subplot(2,2,1)

image(R.LongitudeLimits,fliplr(R.LatitudeLimits),A)
grid on

set(gca,'TickDir','out','ydir','normal','xlim',xlim,'ylim',ylim, ...
    'TickLabelInterpreter','Latex','FontSize',fs, ...
    'XTick',-117:.25:-115,'YTick',32:.25:35);

ylabel('Latitude ($^\circ$N)','Interpreter','Latex','FontSize',fs)
xlabel('Longitude ($^\circ$E)','Interpreter','Latex','FontSize',fs)

% proper aspect ratio
% posM = get(gca,'Position');
% set(gca,'Position',[posM(1) posM(2) pos(3) pos(4)]);

% get aspect ratio correct
pos = get(gca,'position');
ratio = 1.2*diff(ylim)/diff(xlim);
pos = [pos(1) pos(2) pos(3)*0.75 pos(4)*ratio*0.75];
set(gca,'position',pos)

% plot site locations
hold on
scatter(-116.034444,33.033056,80,'wo','fill','MarkerEdgeColor','k')
scatter(-115.832778,33.190556,80,'ws','fill','MarkerEdgeColor','k','Marker','^')
scatter(-115.86102,33.12624,80,'ro','fill','MarkerEdgeColor','k')
hold off
grid on

% title
ttl = title([sat ' VIIRS: 26-May-2026'],'Interpreter','Latex', ...
    'FontSize',fs);

% text(-117.42,34.05,'a.','Interpreter','Latex','FontSize',fs, ...
%     'BackgroundColor','w','EdgeColor','k')


%% save
outDir = fullfile(projectRoot,'output');
if ~isfolder(outDir), mkdir(outDir), end
exportgraphics(gca,fullfile(outDir,'figLocation.jpg'),'Resolution',300)
exportgraphics(gca,fullfile(outDir,'figLocation.pdf'),'ContentType','vector')
savefig(gcf,fullfile(outDir,'figLocation.fig'))
