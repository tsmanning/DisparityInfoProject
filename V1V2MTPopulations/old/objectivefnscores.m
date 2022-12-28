%%%%%%%%%%%%%%%%%%%%%%%
% FI figures of merit %
%%%%%%%%%%%%%%%%%%%%%%%
% calculate scores for objective functions across brain regions
% EA Sept 2022
% known issues:
% - double check FI inputs
% - double check scene stats
% - try several normalizations


clear all;
close all;
clc;

%% load in FI data
dat = load('results_GaussFit.mat');
x = dat.x;
v1 = dat.v1;
v2 = dat.v2;
mt = dat.mt;

figure; hold on;
plot(x,v1,'color','b');
plot(x,v2,'color','g');
plot(x,mt,'color','r');
legend('V1','V2','MT');
xlabel('disparity (D)');
ylabel('FI (normalized)');

%% load in probability distributions

 probdir = '/Users/emily/Dropbox/FisherDisparityProject/FisherDisparityProjectShared/SceneStatsAnalysis/';
%probdir = 'C:/Users/alexa/Dropbox/FisherDisparityProjectShared/SceneStatsAnalysis/2022_sceneStatsCode_Tyler/savedImageStats_BORISdataset/';

% resolution used inTyler's scene stats
res = 52;
lb  = -2;
ub  = 2;edges_disp = linspace(lb,ub,res);
cntr_disp  = edges_disp(1:end-1) + diff(edges_disp(1:2))/2;

% load in Tyler's scene stats analysis
sampling = {'Circ','V1','V2','MT'};

for s = 1:length(sampling)
    
    % load data
    sndo = load([probdir,'savedImageStats_BORISdataset/dispHist' sampling{s} '_sando.mat']);
    wlk = load([probdir,'savedImageStats_BORISdataset/dispHist' sampling{s} '_walking.mat']);
       
    if s == 1
        Sando(s,:)  = sndo.dispHistCirc./sum(sndo.dispHistCirc);
        Walk(s,:)   = wlk.dispHistCirc./sum(wlk.dispHistCirc);
    elseif s == 2
        Sando(s,:)  = sndo.dispHistV1./sum(sndo.dispHistV1);    % sando scene stats
        Walk(s,:)   = wlk.dispHistV1./sum(wlk.dispHistV1);      % walking scene stats
        FI(s,:)     = dat.v1;                                   % V1 FI
        lfi(s,:,:)  = dat.lv1';                                  % bootstrapped V1 FI samples
        clr = 'b';
    elseif s == 3
        Sando(s,:)  = sndo.dispHistV2./sum(sndo.dispHistV2);
        Walk(s,:)   = wlk.dispHistV2./sum(wlk.dispHistV2);
        FI(s,:)     = dat.v2;
        lfi(s,:,:)  = dat.lv2';
        clr = 'g';
    elseif s == 4
        Sando(s,:)  = sndo.dispHistMT./sum(sndo.dispHistMT);
        Walk(s,:)   = wlk.dispHistMT./sum(wlk.dispHistMT);
        FI(s,:)     = dat.mt;
        lfi(s,:,:)  = dat.lmt';
        clr = 'r';
    end
    
    % if using the resampled V1, and MT data, all disparity stats should match V2
    if strcmp(flag,'_resampledV1MT') && s > 1
        sndo = load([probdir,'dispHistV2_sando.mat']);
        wlk = load([probdir,'dispHistV2_walking.mat']);
        Sando(s,:)  = sndo.dispHistV2./sum(sndo.dispHistV2);
        Walk(s,:)   = wlk.dispHistV2./sum(wlk.dispHistV2);
    end
   
    
end

%% normalize and score

close all;

% FI sums to 1,1000
for c = [1,1000]
    figure; 
    plotscores(c*v1,c*v2,c*mt,Sando,Walk);
    title(['FI sums to ',string(c)]);
    csums = [sum(c*v1),sum(c*v2),sum(c*mt)]
end

tmp = load('/Users/emily/Dropbox/FisherDisparityProject/FisherDisparityProjectShared/OLDSceneStatsAnalysis/DataFromTyler/dispHistV2_sando.mat');

%tmp.dispHistV2 = (tmp.dispHistV2/sum(tmp.dispHistV2);

figure; hold on;
subplot(1,3,1); hold on;
plot(x,v1,'k')
plot(x,(Walk(2,:).^2)/sum(Walk(2,:).^2),'r')
ylim([0 .25])
subplot(1,3,2); hold on;
plot(x,v2,'k')
plot(x,(Walk(3,:).^2)/sum(Walk(3,:).^2),'r')
%plot(x,tmp.dispHistV2.,'b')
ylim([0 .25])
subplot(1,3,3); hold on;
plot(x,mt,'k')
plot(x,(Walk(4,:).^2)/sum(Walk(4,:).^2),'r')
ylim([0 .25])



keyboard

%% sqrt FI sums to constant
sqrtsumsorig = [sum(sqrt(v1)),sum(sqrt(v2)),sum(sqrt(mt))]
%for c = [1,10,100,1000,10000]
c = 100;
    cv1 = v1*c^2/sqrtsumsorig(1)^2;
    cv2 = v2*c^2/sqrtsumsorig(2)^2;
    cmt = mt*c^2/sqrtsumsorig(3)^2;
    figure; 
    plotscores(cv1,cv2,cmt,Sando,Walk);
    title(['sqrt FI sums to ',string(c)]);
    sqrtsums = [sum(sqrt(cv1)),sum(sqrt(cv2)),sum(sqrt(cmt))]
%end

%% scoring functions

function score = infoscore(FI,p) % this is fake, change!
    score = sum(p.*log(FI));
end

function score = infoscore_forreal(FI,p)
    score = sum(p.*log(FI));
end

function score = discscore(FI,p)
    score = sum(-p./FI);
end

function plotscores(FI1,FI2,FI3,Sando,Walk)
    subplot(1,2,1);
    hold on;
    h1 = plot([1,2,3],[infoscore(FI1,Sando(2,:)),infoscore(FI2,Sando(3,:)),infoscore(FI3,Sando(4,:))],'-k');
    h2 = plot([1,2,3],[infoscore(FI1,Walk(2,:)),infoscore(FI2,Walk(3,:)),infoscore(FI3,Walk(4,:))],'--k');
    plot(1,infoscore(FI1,Sando(2,:)),'og')
    plot(2,infoscore(FI2,Sando(3,:)),'ob')
    plot(3,infoscore(FI3,Sando(4,:)),'or')
    plot(1,infoscore(FI1,Walk(2,:)),'og')
    plot(2,infoscore(FI2,Walk(3,:)),'ob')
    plot(3,infoscore(FI3,Walk(4,:)),'or')
    xticks([1,2,3])
    xticklabels(['V1';'V2';'MT'])
    ylabel('infomax score (higher is better)')
    
    subplot(1,2,2);
    hold on
    h1 = plot([1,2,3],[discscore(FI1,Sando(2,:)),discscore(FI2,Sando(3,:)),discscore(FI3,Sando(4,:))],'-k');
    h2 = plot([1,2,3],[discscore(FI1,Walk(2,:)),discscore(FI2,Walk(3,:)),discscore(FI3,Walk(4,:))],'--k');
    plot(1,discscore(FI1,Sando(2,:)),'og')
    plot(2,discscore(FI2,Sando(3,:)),'ob')
    plot(3,discscore(FI3,Sando(4,:)),'or')
    plot(1,discscore(FI1,Walk(2,:)),'og')
    plot(2,discscore(FI2,Walk(3,:)),'ob')
    plot(3,discscore(FI3,Walk(4,:)),'or')
    legend([h1,h2],{'sando','walk'})
    xticks([1,2,3])
    xticklabels(['V1';'V2';'MT'])
    ylabel('discrimax score (higher is better)')
end

