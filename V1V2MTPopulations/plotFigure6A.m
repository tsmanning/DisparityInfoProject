% Generate a pair of example Gabor tuning curves (original and resampled)

clear all
close all

splPath = regexp(which('main_reparameterizeV1.m'),filesep,'split');
topDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
saveDir = [topDir,'plots/Fig6A/'];

saveOn = 1;

% Define Gabor tuning function
TC = @(disp,r0,A,prPosDisp,sig,prFreq,prPhase) ...
            r0 + A*exp(-0.5*((disp - prPosDisp)/sig).^2 ).*...
                   cos(2*pi*prFreq*(disp - prPosDisp) + prPhase);

% Original V1 Cell
d   = linspace(-2,2,1000);
r0  = 30;
A   = 40;
pd  = 0;
sig = 1;
pf1 = 0.4;
pPh = 0.8;

% Resampled V1 Cell (only Disp. F. from MT distribution)
pf2 = 0.25;

% Define tuning functions
gabFxn1 = TC(d,r0,A,pd,sig,pf1,pPh);
gabFxn2 = TC(d,r0,A,pd,sig,pf2,pPh);


%% Plot
ymax = 80;

fig1 = figure;
fig1.Position = [100 500 650 600];

hold on;

plot(d,gabFxn1,'k','linewidth',4);
plot([-2 2],r0*[1 1],'--k','linewidth',2);

set(gca,'xlim',[-2 2],'ylim',[0 ymax],'xtick',[-2 -1 0 1 2],'ytick',[],'fontsize',20);
xlabel('Disparity (\circ)');
ylabel('Response (sp/s)');

fig2 = figure;
fig2.Position = [800 500 650 600];

hold on;

plot(d,gabFxn2,'k','linewidth',4);
plot([-2 2],r0*[1 1],'--k','linewidth',2);

set(gca,'xlim',[-2 2],'ylim',[0 ymax],'xtick',[-2 -1 0 1 2],'ytick',[],'fontsize',20);
xlabel('Disparity (\circ)');
ylabel('Response (sp/s)');


%% Save plots
if saveOn

    if ~exist(saveDir)
        mkdir(saveDir)
    end

    saveas(fig1,[saveDir,'originalTC.svg']);
    saveas(fig2,[saveDir,'resampledTC.svg']);
    
end