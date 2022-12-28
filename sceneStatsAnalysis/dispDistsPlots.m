% Plot disparity probability distributions used for main figs in paper and controls

clear all 
close all

%% Load in KSDs

% Define directory structure
splPath = regexp(which('dispDistsPlots'),filesep,'split');
rootDir = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
datDir  = [rootDir,'savedImageStats_BORISdataset',filesep];
figDir  = [rootDir,'plots',filesep];

addpath([rootDir,'sharedTools']);

% Load in disparity histograms based on each KSD
load([datDir,'dispHistStruct.mat']);

% Define stats subsets
imgSet      = {'sando','walking'};
imRegionF   = {'UpperVF','LowerVF','Central','Peripheral','all'};
region      = {'Circ','V1','V2','MT','V1V2Rect'};

% Unused VF breakdowns
% imRegion    = {'Upper VF_','Lower VF_','Central_','Peripheral_',''};
% imRegionLCI = {'UpperVFLCI','LowerVFLCI','CentralLCI','PeripheralLCI','allLCI'};
% imRegionUCI = {'UpperVFUCI','LowerVFUCI','CentralUCI','PeripheralUCI','allUCI'};
% suff        = {'_m1','_m2','_m2','_m1',''};

% Toggle saving figures on/off
saveOn = 0;


%% Plot disparity distributions and CIs
res = 52;
ub  = 2;
lb  = -2;

edgesDisp = linspace(lb,ub,res);
cntrDisp  = edgesDisp(1:end-1) + diff(edgesDisp(1:2))/2;

% Define some more pleasant colors
blueCI  = ColorIt('b');
greenCI = ColorIt('g');
redCI   = ColorIt('r');

fCount = 1;
f = cell(4,1);

% Loop over behavioral image sets
for ii = 1:2

    % Plot stats using full MT dataset
    f{fCount} = figure;
    f{fCount}.Position = [100 700 650 600];
    hold on;

    plot([0 0],[0 5],'--k','linewidth',2);

    fill([cntrDisp fliplr(cntrDisp) cntrDisp(1)],...
        [dispDat.(imgSet{ii}).(region{1}).('allUCI') fliplr(dispDat.(imgSet{ii}).(region{1}).('allLCI'))...
        dispDat.(imgSet{ii}).(region{1}).('allUCI')(1)],'k','edgecolor','none','facealpha',0.5);
    fill([cntrDisp fliplr(cntrDisp) cntrDisp(1)],...
        [dispDat.(imgSet{ii}).(region{2}).('allUCI') fliplr(dispDat.(imgSet{ii}).(region{2}).('allLCI'))...
        dispDat.(imgSet{ii}).(region{2}).('allUCI')(1)],blueCI,'edgecolor','none','facealpha',0.5);
    fill([cntrDisp fliplr(cntrDisp) cntrDisp(1)],...
        [dispDat.(imgSet{ii}).(region{3}).('allUCI') fliplr(dispDat.(imgSet{ii}).(region{3}).('allLCI'))...
        dispDat.(imgSet{ii}).(region{3}).('allUCI')(1)],greenCI,'edgecolor','none','facealpha',0.5);
    fill([cntrDisp fliplr(cntrDisp) cntrDisp(1)],...
        [dispDat.(imgSet{ii}).(region{4}).('allUCI') fliplr(dispDat.(imgSet{ii}).(region{4}).('allLCI'))...
        dispDat.(imgSet{ii}).(region{4}).('allUCI')(1)],redCI,'edgecolor','none','facealpha',0.5);

    p(1) = plot(cntrDisp,realDat.(imgSet{ii}).(region{1}).('all'),'color',[0 0 0],'linewidth',2);
    p(2) = plot(cntrDisp,realDat.(imgSet{ii}).(region{2}).('all'),'color',blueCI,'linewidth',2);
    p(3) = plot(cntrDisp,realDat.(imgSet{ii}).(region{3}).('all'),'color',greenCI,'linewidth',2);
    p(4) = plot(cntrDisp,realDat.(imgSet{ii}).(region{4}).('all'),'color',redCI,'linewidth',2);

    set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1]);
    legend(p,{'Cent 10\circ','V1','V2','MT'});
    xlabel('Horizontal disparity (\circ)');
    ylabel('Probability density');
    title([imgSet{ii}]);

    fCount = fCount + 1;

    % Plot stats using subsampled MT dataset restricted to V1/V2 bounds
    f{fCount} = figure;
    f{fCount}.Position = [700 700 650 600];
    hold on;

    plot([0 0],[0 5],'--k','linewidth',2);

    fill([cntrDisp fliplr(cntrDisp) cntrDisp(1)],...
        [dispDat.(imgSet{ii}).(region{1}).('allUCI') fliplr(dispDat.(imgSet{ii}).(region{1}).('allLCI'))...
        dispDat.(imgSet{ii}).(region{1}).('allUCI')(1)],'k','edgecolor','none','facealpha',0.5);
    fill([cntrDisp fliplr(cntrDisp) cntrDisp(1)],...
        [dispDat.(imgSet{ii}).(region{2}).('allUCI') fliplr(dispDat.(imgSet{ii}).(region{2}).('allLCI'))...
        dispDat.(imgSet{ii}).(region{2}).('allUCI')(1)],'b','edgecolor','none','facealpha',0.5);
    fill([cntrDisp fliplr(cntrDisp) cntrDisp(1)],...
        [dispDat.(imgSet{ii}).(region{3}).('allUCI') fliplr(dispDat.(imgSet{ii}).(region{3}).('allLCI'))...
        dispDat.(imgSet{ii}).(region{3}).('allUCI')(1)],'g','edgecolor','none','facealpha',0.5);
    fill([cntrDisp fliplr(cntrDisp) cntrDisp(1)],...
        [dispDat.(imgSet{ii}).(region{5}).('allUCI') fliplr(dispDat.(imgSet{ii}).(region{5}).('allLCI'))...
        dispDat.(imgSet{ii}).(region{5}).('allUCI')(1)],'r','edgecolor','none','facealpha',1);

    p(1) = plot(cntrDisp,realDat.(imgSet{ii}).(region{1}).('all'),'color',[0 0 0],'linewidth',2);
    p(2) = plot(cntrDisp,realDat.(imgSet{ii}).(region{2}).('all'),'color',blueCI,'linewidth',2);
    p(3) = plot(cntrDisp,realDat.(imgSet{ii}).(region{3}).('all'),'color',greenCI,'linewidth',2);
    p(4) = plot(cntrDisp,realDat.(imgSet{ii}).(region{5}).('all'),'color',redCI,'linewidth',2);

    set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1]);
    legend(p,{'Cent 10\circ','V1','V2','MT (restricted)'});
    xlabel('Horizontal disparity (\circ)');
    ylabel('Probability density');
    title([imgSet{ii}]);

    fCount = fCount + 1;

end


%% Save figs

if saveOn

    saveas(f{1},[figDir,'sando.svg']);
    saveas(f{2},[figDir,'sandoCntl.svg']);
    saveas(f{3},[figDir,'walking.svg']);
    saveas(f{4},[figDir,'walkingCntl.svg']);

end

