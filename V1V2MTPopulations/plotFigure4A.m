% Generate an example Gabor tuning curve + component functions

clear all
close all

splPath = regexp(which('main_reparameterizeV1.m'),filesep,'split');
topDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
saveDir = [topDir,'plots/Fig4A/'];

saveOn = 1;

% Define Gabor tuning function and component functions
TC = @(disp,r0,A,prPosDisp,sig,prFreq,prPhase) ...
            r0 + A*exp(-0.5*((disp - prPosDisp)/sig).^2 ).*...
                   cos(2*pi*prFreq*(disp - prPosDisp) + prPhase);

Gaus = @(disp,r0,A,prPosDisp,sig) ...
            r0 + A*exp(-0.5*((disp - prPosDisp)/sig).^2 );

cosF = @(disp,r0,A,prPosDisp,prFreq,prPhase) ...
            r0 + A*cos(2*pi*prFreq*(disp - prPosDisp) + prPhase);

% Define example tuning curve parameters
d   = linspace(-2,2,1000);
r0  = 30;
A   = 40;
pd  = 0;
sig = 1;
pf  = 0.2;
pPh = 0.8;

% Generate tuning curves
gabFxn = TC(d,r0,A,pd,sig,pf,pPh);
gauFxn = Gaus(d,0,1,pd,sig);
cosFxn = cosF(d,0,1,pd,pf,pPh);


%% Plot
ymax = 80;

fig1 = figure;
fig1.Position = [100 100 650 600];

hold on;

plot(d,gabFxn,'k','linewidth',4);
plot([-2 2],r0*[1 1],'--k','linewidth',2);

set(gca,'xlim',[-2 2],'ylim',[0 ymax],'xtick',[-2 -1 0 1 2],'ytick',[],'fontsize',20);
xlabel('Disparity (\circ)');
ylabel('Response (sp/s)');

fig2 = figure;
fig2.Position = [100 100 650 600];

hold on;

plot(d,cosFxn,'--k','linewidth',4);
plot(-pPh/(2*pi*pf)*[1 1],[-1 1],'--k','linewidth',2);
plot(0*[1 1],[-1 1],'--k','linewidth',2);
plot([-2 2],0*[1 1],'--k','linewidth',2);

set(gca,'xlim',[-2 2],'ylim',[-1 1],'xtick',[-2 -1 0 1 2],'ytick',[],'fontsize',20);
xlabel('Disparity (\circ)');
ylabel('Response (sp/s)');

fig3 = figure;
fig3.Position = [100 100 650 600];

hold on;

plot(d,gauFxn,'--k','linewidth',4);
plot(pd*[1 1],[0 1],'--k','linewidth',2);
plot(sig*[1 1],[0 1],'--k','linewidth',2);

set(gca,'xlim',[-2 2],'ylim',[0 1],'xtick',[-2 -1 0 1 2],'ytick',[],'fontsize',20);
xlabel('Disparity (\circ)');
ylabel('Response (sp/s)');


%% Save plots
if saveOn

    if ~exist(saveDir)
        mkdir(saveDir)
    end

    saveas(fig1,[saveDir,'fullTC.svg']);
    saveas(fig2,[saveDir,'cosTC.svg']);
    saveas(fig3,[saveDir,'gauTC.svg']);

end