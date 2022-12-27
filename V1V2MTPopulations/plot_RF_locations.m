function plot_RF_locations(V1,V2,MT,flag,topDir)

%% RECEPTIVE FIELD LOCATIONS

saveOn = 1;

splPath = regexp(which('plot_RF_locations'),filesep,'split');
ksdDir  = [filesep,fullfile(splPath{1:numel(splPath)-2}),filesep,'SceneStatsAnalysis/savedKSDmatFiles_BORISdataset/'];

f1 = figure; 
f1.Position = [100 100 650 600];
hold on; 
title('Receptive field locations');
scatter(V1.x_pos,V1.y_pos,[],ColorIt('b'),'filled');
scatter(V2.x_pos,V2.y_pos,[],ColorIt('g'),'filled');
scatter(MT.x_pos,MT.y_pos,[],ColorIt('r'),'filled');
plot([-12 12],[0 0],'k:');
plot([0 0],[-12 12],'k:');
axis([-12 12 -12 12]); axis equal tight;
legend('V1','V2','MT');
xlabel('horizontal eccentricity (deg)'); 
ylabel('vertical eccentricity (deg)'); 
box on;
set(gca,'fontsize',20);

% % plot as ksdensity
% figure; hold on;
% subplot(1,3,1); hold on;
% [v1f,v1xi] = ksdensity([V1.x_pos',V1.y_pos']);
% imagesc(v1xi,v1f)

% warning('add in manual ksdensity -- get code from Tyler for how he made these')

% f2 = figure; hold on;
% subplot(1,3,1); hold on; title('V1 RF density')
% load([ksdDir,'V1densityMat_BORIS.mat']);
% imagesc([-10 10],[-10 10],V1densityMat);
% axis image;
% 
% subplot(1,3,2); hold on; title('V2 RF density')
% load([ksdDir,'V2densityMat_BORIS.mat']);
% imagesc([-10 10],[-10 10],V2densityMat);
% axis image;
% 
% subplot(1,3,3); hold on; title('V2 RF density')
% load([ksdDir,'MTdensityMat_BORIS.mat']);
% imagesc([-10 10],[-10 10],MTdensityMat);
% axis image;

% Save images
if saveOn
%     saveas(f1,[topDir,'plots/RFLocations/RFlocations' flag '.png']);
    saveas(f1,[topDir,'plots/RFLocations/RFlocations' flag '.svg']);

%     saveas(f2,[topDir,'plots/RFLocations/RFlocationsDensity' flag '.png']);
%     saveas(f2,[topDir,'plots/RFLocations/RFlocationsDensity' flag '.eps'],'epsc');
end

end
