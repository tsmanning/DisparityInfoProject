% Plot kernel-smoothed densities used for main figs in paper and controls

clear all 
close all

% Load in KSDs
splPath = regexp(which('KSDplots'),filesep,'split');
ksdDir  = [filesep,fullfile(splPath{1:numel(splPath)-2}),filesep,'SceneStatsAnalysis/savedKSDmatFiles_BORISdataset/'];
saveDir = [filesep,fullfile(splPath{1:numel(splPath)-2}),filesep,'V1V2MTPopulations/plots/RFLocations/'];

load([ksdDir,'/V1densityMat_BORIS.mat']);
load([ksdDir,'/V2densityMat_BORIS.mat']);
load([ksdDir,'/MTdensityMat_BORIS.mat']);
load([ksdDir,'/V1V2Rect_BORIS.mat']);

% Plot
tix  = [-10 -5 0 5 10];
lims = [-10 10];

cmap = parula;
gam  = 255*linspace(0,1,256).^0.25 + 1;
cmap = [interp1(1:256,cmap(:,1),gam)' interp1(1:256,cmap(:,2),gam)' interp1(1:256,cmap(:,3),gam)'];

% V1
f1 = figure;
f1.Position = [100 100 760 600];

hold on

imagesc(tix,tix,V1densityMat); 
axis image; 
axis xy;

scatter(V1.x_pos,V1.y_pos,30,'w','filled','markerfacealpha',0.5);

set(gca,'PlotBoxAspectRatio',[1 1 1],'xlim',lims,'ylim',lims,'xtick',tix,'ytick',tix,'fontsize',20);
xlabel('horizontal eccentricity (deg)');
ylabel('vertical eccentricity (deg)');
title('V1');
cb1 = colorbar;
cb1.Label.String = 'probability';
clim([0 0.06]);
colormap(cmap);

% V2
f2 = figure;
f2.Position = [700 100 760 600];

hold on

imagesc(tix,tix,V2densityMat); 
axis image; 
axis xy;

scatter(V2.x_pos,V2.y_pos,30,'w','filled','markerfacealpha',0.5);

set(gca,'PlotBoxAspectRatio',[1 1 1],'xlim',lims,'ylim',lims,'xtick',tix,'ytick',tix,'fontsize',20);
xlabel('horizontal eccentricity (deg)');
ylabel('vertical eccentricity (deg)');
title('V2');
cb2 = colorbar;
cb2.Label.String = 'probability';
clim([0 0.06]);
colormap(cmap);

% MT
f3 = figure;
f3.Position = [1400 100 760 600];

hold on

imagesc(tix,tix,MTdensityMat); 
axis image; 
axis xy;

scatter(MT.x_pos,MT.y_pos,30,'w','filled','markerfacealpha',0.5);

set(gca,'PlotBoxAspectRatio',[1 1 1],'xlim',lims,'ylim',lims,'xtick',tix,'ytick',tix,'fontsize',20);
xlabel('horizontal eccentricity (deg)');
ylabel('vertical eccentricity (deg)');
title('MT');
cb3 = colorbar;
cb3.Label.String = 'probability';
clim([0 0.06]);
colormap(cmap);

% Restricted MT
f4 = figure;
f4.Position = [2100 100 760 600];

hold on

imagesc(tix,tix,V1V2RectMat); 
axis image; 
axis xy;

scatter(MT.x_posBnd,MT.y_posBnd,20,'w','filled','markerfacealpha',0.5);

set(gca,'PlotBoxAspectRatio',[1 1 1],'xlim',lims,'ylim',lims,'xtick',tix,'ytick',tix,'fontsize',20);
xlabel('horizontal eccentricity (deg)');
ylabel('vertical eccentricity (deg)');
title('MT (restricted)');
cb3 = colorbar;
cb3.Label.String = 'probability';
clim([0 0.06]);
colormap(cmap);

% Save
saveas(f1,[saveDir,'V1_KSD.svg']);
saveas(f2,[saveDir,'V2_KSD.svg']);
saveas(f3,[saveDir,'MT_KSD.svg']);
saveas(f4,[saveDir,'MTrestricted_KSD.svg']);




