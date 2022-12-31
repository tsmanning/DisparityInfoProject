function [] = dispDistsPlots(saveOn)
% Plot disparity probability distributions used for main figs in paper and controls

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

    p(1) = shadeplot(realDat.(imgSet{ii}).(region{2}).('all'),[dispDat.(imgSet{ii}).(region{2}).('allLCI'); ...
                     dispDat.(imgSet{ii}).(region{2}).('allUCI')],cntrDisp,blueCI,0.5,4);
    p(2) = shadeplot(realDat.(imgSet{ii}).(region{3}).('all'),[dispDat.(imgSet{ii}).(region{3}).('allLCI'); ...
                     dispDat.(imgSet{ii}).(region{3}).('allUCI')],cntrDisp,greenCI,0.5,4);
    p(3) = shadeplot(realDat.(imgSet{ii}).(region{4}).('all'),[dispDat.(imgSet{ii}).(region{4}).('allLCI'); ...
                     dispDat.(imgSet{ii}).(region{4}).('allUCI')],cntrDisp,redCI,0.5,4);

    set(gca,'fontsize',25,'plotboxaspectratio',[1 1 1],'ylim',[0 5],'xtick',[-2:1:2]);
    legend(p,{'V1','V2','MT'});
    xlabel('Horizontal disparity (\circ)');
    ylabel('Probability density');
    title([imgSet{ii}]);

    fCount = fCount + 1;

    % Plot stats using subsampled MT dataset restricted to V1/V2 bounds
    f{fCount} = figure;
    f{fCount}.Position = [700 700 650 600];
    hold on;

    plot([0 0],[0 5],'--k','linewidth',2);

    p(1) = shadeplot(realDat.(imgSet{ii}).(region{2}).('all'),[dispDat.(imgSet{ii}).(region{2}).('allLCI'); ...
                     dispDat.(imgSet{ii}).(region{2}).('allUCI')],cntrDisp,blueCI,0.5,4);
    p(2) = shadeplot(realDat.(imgSet{ii}).(region{3}).('all'),[dispDat.(imgSet{ii}).(region{3}).('allLCI'); ...
                     dispDat.(imgSet{ii}).(region{3}).('allUCI')],cntrDisp,greenCI,0.5,4);
    p(3) = shadeplot(realDat.(imgSet{ii}).(region{5}).('all'),[dispDat.(imgSet{ii}).(region{5}).('allLCI'); ...
                     dispDat.(imgSet{ii}).(region{5}).('allUCI')],cntrDisp,redCI,0.5,4);

    set(gca,'fontsize',25,'plotboxaspectratio',[1 1 1],'ylim',[0 5],'xtick',[-2:1:2]);
    legend(p,{'V1','V2','MT (restricted)'});
    xlabel('Horizontal disparity (\circ)');
    ylabel('Probability density');
    title([imgSet{ii}]);

    fCount = fCount + 1;

end


%% Save figs

if saveOn

    if ~exist(figDir)
        mkdir(figDir)
    end

    saveas(f{1},[figDir,'sando.svg']);
    saveas(f{2},[figDir,'sandoCntl.svg']);
    saveas(f{3},[figDir,'walking.svg']);
    saveas(f{4},[figDir,'walkingCntl.svg']);

end

end
