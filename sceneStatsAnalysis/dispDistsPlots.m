% Plot disparity probability distributions used for main figs in paper and controls

clear all 
close all

%% Load in KSDs'allUCI'
ksdDir = ['/home/tyler/Documents/MATLAB/cooperLab/2-Modeling_Simulations/disparityMTModel/',...
          '5-sceneStatsAnalysis/'];
realDatDir = [ksdDir,'imgStats/'];
datDir = [ksdDir,'borisStats/'];

imgSet      = {'sando','walking'};
% imRegion    = {'Upper VF_','Lower VF_','Central_','Peripheral_',''};
imRegionF   = {'UpperVF','LowerVF','Central','Peripheral','all'};
% imRegionLCI = {'UpperVFLCI','LowerVFLCI','CentralLCI','PeripheralLCI','allLCI'};
% imRegionUCI = {'UpperVFUCI','LowerVFUCI','CentralUCI','PeripheralUCI','allUCI'};
region      = {'Circ','V1','V2','MT','V1V2Rect'};
% suff        = {'_m1','_m2','_m2','_m1',''};


%% Plot
res = 52;
ub  = 2;
lb  = -2;

edgesDisp = linspace(lb,ub,res);
cntrDisp  = edgesDisp(1:end-1) + diff(edgesDisp(1:2))/2;

% Main figs
ii = 1;

% Sando
thisCirc   = load([realDatDir,'dispHist',region{1},'_','',imgSet{ii}]);
circDat    = thisCirc.(['dispHistCirc','']);
thisV1     = load([realDatDir,'dispHist',region{2},'_','',imgSet{ii}]);
v1Dat      = thisV1.(['dispHistV1','']);
thisV2     = load([realDatDir,'dispHist',region{3},'_','',imgSet{ii}]);
v2Dat      = thisV2.(['dispHistV2','']);
thisMT     = load([realDatDir,'dispHist',region{4},'_','',imgSet{ii}]);
MTDat      = thisMT.(['dispHistMT','']);
thisMTrest = load([realDatDir,'dispHist',region{5},'_','',imgSet{ii}]);
MTDatrest  = thisMTrest.(['dispHistV1V2Rect','']);

load([datDir,'dispHistStruct.mat']);

f1 = figure;
f1.Position = [100 100 650 600];
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
     [dispDat.(imgSet{ii}).(region{4}).('allUCI') fliplr(dispDat.(imgSet{ii}).(region{4}).('allLCI'))...
      dispDat.(imgSet{ii}).(region{4}).('allUCI')(1)],'r','edgecolor','none','facealpha',0.5);

p(1) = plot(cntrDisp,circDat,'color',[0 0 0],'linewidth',2);
p(2) = plot(cntrDisp,v1Dat,'color',[0 0 1],'linewidth',2);
p(3) = plot(cntrDisp,v2Dat,'color',[0 1 0],'linewidth',2);
p(4) = plot(cntrDisp,MTDat,'color',[1 0 0],'linewidth',2);

set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1]);
legend(p,{'Cent 10\circ','V1','V2','MT'});
xlabel('Horizontal disparity (\circ)');
ylabel('Probability density');
title([imgSet{ii}]);

load([datDir,'dispHistStructCntl.mat']);

f2 = figure;
f2.Position = [700 100 650 600];
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
      dispDat.(imgSet{ii}).(region{5}).('allUCI')(1)],'r','edgecolor','none','facealpha',0.5);

p(1) = plot(cntrDisp,circDat,'color',[0 0 0],'linewidth',2);
p(2) = plot(cntrDisp,v1Dat,'color',[0 0 1],'linewidth',2);
p(3) = plot(cntrDisp,v2Dat,'color',[0 1 0],'linewidth',2);
p(4) = plot(cntrDisp,MTDatrest,'color',[1 0 0],'linewidth',2);

set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1]);
legend(p,{'Cent 10\circ','V1','V2','MT (restricted)'});
xlabel('Horizontal disparity (\circ)');
ylabel('Probability density');
title([imgSet{ii}]);

saveas(f1,[ksdDir,'2022_sceneStatsCode_Tyler/sando.svg']);
saveas(f2,[ksdDir,'2022_sceneStatsCode_Tyler/sandoCntl.svg']);

if 1
% Walking
ii = 2;

thisCirc   = load([realDatDir,'dispHist',region{1},'_','',imgSet{ii}]);
circDat    = thisCirc.(['dispHistCirc','']);
thisV1     = load([realDatDir,'dispHist',region{2},'_','',imgSet{ii}]);
v1Dat      = thisV1.(['dispHistV1','']);
thisV2     = load([realDatDir,'dispHist',region{3},'_','',imgSet{ii}]);
v2Dat      = thisV2.(['dispHistV2','']);
thisMT     = load([realDatDir,'dispHist',region{4},'_','',imgSet{ii}]);
MTDat      = thisMT.(['dispHistMT','']);
thisMTrest = load([realDatDir,'dispHist',region{5},'_','',imgSet{ii}]);
MTDatrest  = thisMTrest.(['dispHistV1V2Rect','']);

load([datDir,'dispHistStruct.mat']);

f3 = figure;
f3.Position = [100 700 650 600];
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
     [dispDat.(imgSet{ii}).(region{4}).('allUCI') fliplr(dispDat.(imgSet{ii}).(region{4}).('allLCI'))...
      dispDat.(imgSet{ii}).(region{4}).('allUCI')(1)],'r','edgecolor','none','facealpha',0.5);

p(1) = plot(cntrDisp,circDat,'color',[0 0 0],'linewidth',2);
p(2) = plot(cntrDisp,v1Dat,'color',[0 0 1],'linewidth',2);
p(3) = plot(cntrDisp,v2Dat,'color',[0 1 0],'linewidth',2);
p(4) = plot(cntrDisp,MTDat,'color',[1 0 0],'linewidth',2);

set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1]);
legend(p,{'Cent 10\circ','V1','V2','MT'});
xlabel('Horizontal disparity (\circ)');
ylabel('Probability density');
title([imgSet{ii}]);

load([datDir,'dispHistStructCntl.mat']);

f4 = figure;
f4.Position = [700 700 650 600];
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

p(1) = plot(cntrDisp,circDat,'color',[0 0 0],'linewidth',2);
p(2) = plot(cntrDisp,v1Dat,'color',[0 0 1],'linewidth',2);
p(3) = plot(cntrDisp,v2Dat,'color',[0 1 0],'linewidth',2);
p(4) = plot(cntrDisp,MTDatrest,'color',[1 0 0],'linewidth',2);

set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1]);
legend(p,{'Cent 10\circ','V1','V2','MT (restricted)'});
xlabel('Horizontal disparity (\circ)');
ylabel('Probability density');
title([imgSet{ii}]);


saveas(f3,[ksdDir,'2022_sceneStatsCode_Tyler/walking.svg']);
saveas(f4,[ksdDir,'2022_sceneStatsCode_Tyler/walkingCntl.svg']);
end
