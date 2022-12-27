% Generate an example Gabor tuning curve 

clear all
close all

TC = @(disp,r0,A,prPosDisp,sig,prFreq,prPhase) ...
            r0 + A*exp(-0.5*((disp - prPosDisp)/sig).^2 ).*...
                   cos(2*pi*prFreq*(disp - prPosDisp) + prPhase);

Gaus = @(disp,r0,A,prPosDisp,sig) ...
            r0 + A*exp(-0.5*((disp - prPosDisp)/sig).^2 );

cosF = @(disp,r0,A,prPosDisp,prFreq,prPhase) ...
            r0 + A*cos(2*pi*prFreq*(disp - prPosDisp) + prPhase);

d   = linspace(-2,2,1000);
r0  = 30;
A   = 40;
pd  = 0;
sig = 1;
pf  = 0.2;
pPh = 0.8;

gabFxn = TC(d,r0,A,pd,sig,pf,pPh);
gauFxn = Gaus(d,0,1,pd,sig);
cosFxn = cosF(d,0,1,pd,pf,pPh);

ymax = 80;

fig = figure;
fig.Position = [100 100 650 600];

hold on;

plot(d,gabFxn,'k','linewidth',4);
plot([-2 2],r0*[1 1],'--k','linewidth',2);

set(gca,'xlim',[-2 2],'ylim',[0 ymax],'xtick',[-2 -1 0 1 2],'ytick',[],'fontsize',20);
xlabel('Disparity (\circ)');
ylabel('Response (sp/s)');

fig = figure;
fig.Position = [100 100 650 600];

hold on;

plot(d,cosFxn,'--k','linewidth',4);
plot(-pPh/(2*pi*pf)*[1 1],[-1 1],'--k','linewidth',2);
plot(0*[1 1],[-1 1],'--k','linewidth',2);
plot([-2 2],0*[1 1],'--k','linewidth',2);

set(gca,'xlim',[-2 2],'ylim',[-1 1],'xtick',[-2 -1 0 1 2],'ytick',[],'fontsize',20);
xlabel('Disparity (\circ)');
ylabel('Response (sp/s)');

fig = figure;
fig.Position = [100 100 650 600];

hold on;

plot(d,gauFxn,'--k','linewidth',4);
plot(pd*[1 1],[0 1],'--k','linewidth',2);
plot(sig*[1 1],[0 1],'--k','linewidth',2);

set(gca,'xlim',[-2 2],'ylim',[0 1],'xtick',[-2 -1 0 1 2],'ytick',[],'fontsize',20);
xlabel('Disparity (\circ)');
ylabel('Response (sp/s)');