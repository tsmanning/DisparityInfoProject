%% Plot example plots from Figure 1

clear all
close all

% Setup dirs
splPath  = regexp(which('Figure1plots.m'),filesep,'split');
topDir   = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
plotsDir = [topDir,'plots',filesep,'Figure1',filesep];

% Toggle save on/off
saveOn = 1;


%% 1A (bottom)
fig1A = figure;
fig1A.Position = [100 100 650 600];
hold on;

x     = linspace(-2,2,100);
tcMus = linspace(-1.65,1.65,5);
sig   = 0.4;

% Plot 5 example Gaussian tuning functions
for ii = 1:5

    tcTint = [1 1 1]*0.75*(ii-1)/5;
    tc = normpdf(x,tcMus(ii),sig);

    plot(x,tc,'color',tcTint,'linewidth',8);

end

set(gca,'xlim',[-2 2],'xtick',-2:2,'ylim',[0 1.1],'ytick',[],'fontsize',30,'plotboxaspectratio',[1 1 1]);
xlabel('Binocular disparity (\circ)'); 
ylabel('Firing rate');


%% 1B 

% These plots are generated in main_disparityStats using the 10deg circular
% sampling distribution


%% 1C

res = 52;
lb  = -2;
ub  = 2;edges_disp = linspace(lb,ub,res);
cntr_disp  = edges_disp(1:end-1) + diff(edges_disp(1:2))/2;

s = 0.3;
p = 1;

p1 = exp( -(abs(cntr_disp)/s).^p );
p1 = p1 / sum(p1);

p2 = p1.^2;
p2 = p2 / sum(p2);

ppt5 = p1.^0.5;
ppt5 = ppt5 / sum(ppt5);

% Plot probabilities to 0.5 and 2 powers
fig1C = figure;
fig1C.Position = [800 100 650 600];
hold on;

plot(cntr_disp, p1,'-','color',[0 0 0],'linewidth',2);
plot(cntr_disp, p2,'--','color',[0 0 0],'linewidth',2);
plot(cntr_disp, ppt5,':','color',[0 0 0],'linewidth',2);

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'xtick',-2:2,'fontsize',30)
axis square; box on;
ylim([0 .27])
legend('p=1','infomax','discrimax');
xlabel('Binocular disparity (\circ)'); 
ylabel('Fisher information (FI)');


%% 1D

dispTix = [-2 0 2];

s = RandStream('mt19937ar','Seed',94720);
RandStream.setGlobalStream(s);

% Define simulated data points
numDat = 80;

xVals = linspace(-1.8,1.8,numDat);

% Set total reconstruction error
totalError = 2;

% Generate some random data for each error function
% - L0 norm: partition into correct and incorrect estimates
nailedIt              = datasample(1:numDat,numDat/2,'replace',false);
didntNailIt           = ones(1,numDat,'logical');
didntNailIt(nailedIt) = didntNailIt(nailedIt) == 0;

l0dat              = 4*rand([1 numDat]) - 2;
l0dat(nailedIt)    = xVals(nailedIt);

% - L2 norm: define error in reconstruction as inverse of squared disparity
randSign              = double(rand([1,numDat])>0.5);
randSign(randSign==0) = -1;
% - [ground truth + (random error magnitude).*(random sign)]
l2dat                 = xVals + (rand([1,numDat]).*xVals.^2) .* (randSign);

% - L2 norm: scale to define total error
% l2dat = totalError*l2dat/sum(l2dat-xVals);

% - L0 norm: scale to define total error
l0dat(didntNailIt) = sum(l2dat-xVals)*l0dat(didntNailIt)/sum(l0dat(didntNailIt)-xVals(didntNailIt));
% l0dat(didntNailIt) = totalError*l0dat(didntNailIt)/sum(l0dat(didntNailIt)-xVals(didntNailIt));

% L0 norm
fig1Di = figure;
fig1Di.Position = [300 300 650 600];
hold on;

plot(x,ones(1,100),'k','linewidth',4);
plot([0 0],[0 1],'k','linewidth',4);

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'xtick',dispTix,'ylim',[0 1],'ytick',[0 1],'fontsize',30);
ylabel('Error Penalty'); 
xlabel('Reconstruction error');

% L0 reconstruction
fig1Dii = figure;
fig1Dii.Position = [980 300 650 600];
hold on;

scatter(xVals,l0dat,200,'k','filled');
plot([-2 2],[-2 2],'--k','linewidth',4);

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'xtick',dispTix,'ylim',[-2 2],'ytick',dispTix,'fontsize',30);
ylabel('Reconstruction'); 
xlabel('Binocular disparity (\circ)');

% L2 norm
fig1Diii = figure;
fig1Diii.Position = [300 800 650 600];
hold on;

plot(x,x.^2,'k','linewidth',4);

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'xtick',dispTix,'ylim',[0 1],'ytick',[0 1],'fontsize',30);
ylabel('Error Penalty'); 
xlabel('Reconstruction error');

% L2 reconstruction
fig1Div = figure;
fig1Div.Position = [980 800 650 600];
hold on;

scatter(xVals,l2dat,200,'k','filled');
plot([-2 2],[-2 2],'--k','linewidth',4);

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'xtick',dispTix,'ylim',[-2 2],'ytick',dispTix,'fontsize',30);
ylabel('Reconstruction'); 
xlabel('Binocular disparity (\circ)');


%% Save figs

if saveOn

    if ~exist(plotsDir)
        mkdir(plotsDir);
    end

    saveas(fig1A,[plotsDir,'fig1A.svg']);

    saveas(fig1C,[plotsDir,'fig1C.svg']);

    saveas(fig1Di,[plotsDir,'fig1Di.svg']);
    saveas(fig1Dii,[plotsDir,'fig1Dii.svg']);
    saveas(fig1Diii,[plotsDir,'fig1Diii.svg']);
    saveas(fig1Div,[plotsDir,'fig1Div.svg']);

end
