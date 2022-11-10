% Assess whether we can account for FI differences between V1 and MT based
% on substituting a single parameter distribution from one to the other

clear all
close all

%% Load in saved Gabor fits
splPath = regexp(which('main_reparameterizeV1.m'),filesep,'split');
topDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
saveDir = [topDir,'plots/reparameterization/'];

addpath(genpath([topDir,'/helper_functions']));

load([topDir,'resultsALL_final.mat']);

% Clip negative FIs (this is kinda problematic?)
V1.FI(V1.FI < 0) = 0;
V2.FI(V2.FI < 0) = 0;
MT.FI(MT.FI < 0) = 0;


%% Resample V1 parameter fits from the MT distribution
numResamps = 100;

for rs = 1:numResamps

    if mod(rs,10) == 0
        disp(['Running resampling run: ',num2str(rs),'/',num2str(numResamps)]);
    end

    V1_off   = permuteFits(V1,V2,MT,'V1','MT','offset');
    V1_amp   = permuteFits(V1,V2,MT,'V1','MT','amplitude');
    V1_envM  = permuteFits(V1,V2,MT,'V1','MT','envMean');
    V1_envW  = permuteFits(V1,V2,MT,'V1','MT','envWidth');
    V1_freq  = permuteFits(V1,V2,MT,'V1','MT','frequency');
    V1_phase = permuteFits(V1,V2,MT,'V1','MT','phase');

    V1_reparams = {V1_off.P, V1_amp.P, V1_envM.P, V1_envW.P, V1_freq.P, V1_phase.P};

    % Recalulate individual cell FIs under reparameterization
    dTC = @(disp,r0,A,prPosDisp,sig,prFreq,prPhase) ...
        -A*exp(-0.5*((disp - prPosDisp)/sig).^2 ).*( ...
        (disp - prPosDisp)/(sig^2).*cos(2*pi*prFreq*(disp - prPosDisp) + prPhase)...
        + (2*pi*prFreq)*sin(2*pi*prFreq*(disp - prPosDisp) + prPhase));
    rng = 2;
    xg1 = [-rng : 0.01 : rng]; %fixed spatial range

    % For each parameter
    for ii = 1:6

        numCells = size(V1_reparams{ii},1);

        % For each cell
        for jj = 1:numCells

            p = V1_reparams{ii}(jj,:);

            % Calculate this tuning curve
            g1 = p(1) + p(2)*exp( -(xg1-p(3)).^2 / (2*p(4)^2) ) .* cos(2*pi*p(5)*(xg1-p(3))+p(6));

            % Calculate tuning curve derivative
            dg1 = dTC(xg1,p(1),p(2),p(3),p(4),p(5),p(6));

            % Calculate Fisher information
            fi = dg1.^2 ./ g1;

            % Remove FI values less than 0
            fi(fi<0) = 0;

            V1_reparamsFI{ii,rs}(jj,:) = fi;

        end

    end
end


%% Calculate FI for each reparameterized V1 distribution & MT (+ Laplacian fits)

% disparities corresponding to each column in FI matrices
x = [-2 : 0.01 : 2];

% resolution used in Tyler's scene stats
res = 52;
lb  = -2;
ub  = 2;
edges_disp = linspace(lb,ub,res);
cntr_disp  = edges_disp(1:end-1) + diff(edges_disp(1:2))/2;

% find indices of x that intersect with x
[~,ia,~] = intersect(round(x,2),round(cntr_disp,2));

% grab just those indices
x = x(ia);
V1.FI = V1.FI(:,ia);
V2.FI = V2.FI(:,ia);
MT.FI = MT.FI(:,ia);

% sum it up to get population FI
v1 = sum(V1.FI); v1 = v1 / sum(v1);
v2 = sum(V2.FI); v2 = v2 / sum(v2);
mt = sum(MT.FI); mt = mt / sum(mt);

% Fit a generalized Laplacian to the population FI
glV1 = fitGeneralizedLaplacian(x,v1);
glV2 = fitGeneralizedLaplacian(x,v2);
glMT = fitGeneralizedLaplacian(x,mt);

% Get summary stats for mt
varMT     = (glMT(1).^2 .*gamma(3./glMT(2)) ) ./gamma(1./glMT(2));
exkurtMT  = ( (gamma(5./glMT(2)).*gamma(1./glMT(2))) ./ (gamma(3./glMT(2)).^2) ) - 3;
entropyMT = 1./glMT(2) - log( glMT(2)./ (2*glMT(1).*gamma(1./glMT(2))) );
    
% Downsample, sum each reparameterized V1 dist, & fit with Laplacian
for rs = 1:numResamps
    for ii = 1:6

        V1_reparamsFIds{ii,rs} = V1_reparamsFI{ii,rs}(:,ia);

        v1_RP(rs,ii,:) = sum(V1_reparamsFIds{ii,rs});
        v1_RP(rs,ii,:) = v1_RP(rs,ii,:) / sum(v1_RP(rs,ii,:));

        glV1_RP(rs,ii,:) = fitGeneralizedLaplacian(x,squeeze(v1_RP(rs,ii,:))');

        varV1(rs,ii)     = (glV1_RP(rs,ii,1).^2 .*gamma(3./glV1_RP(rs,ii,2)) ) ./gamma(1./glV1_RP(rs,ii,2));
        exkurtV1(rs,ii)  = ( (gamma(5./glV1_RP(rs,ii,2)).*gamma(1./glV1_RP(rs,ii,2))) ./ (gamma(3./glV1_RP(rs,ii,2)).^2) ) - 3;
        entropyV1(rs,ii) = 1./glV1_RP(rs,ii,2) - log( glV1_RP(rs,ii,2)./ (2*glV1_RP(rs,ii,1).*gamma(1./glV1_RP(rs,ii,2))) );

        JSD(rs,ii)       = getJSDiv(squeeze(v1_RP(rs,ii,:))',mt);

    end
end


%% Find the mean and std (95% CIs?) of the bootstrapped FIs, Laplace fits and JS divergence

meanV1RP = squeeze(mean(v1_RP,1));
stdV1RP  = squeeze(std(v1_RP,[],1));

meanglV1RP = squeeze(mean(glV1_RP,1));
stdglV1RP  = squeeze(std(glV1_RP,[],1));

medianVar     = median(varV1,1);
iqrVar        = prctile(varV1,[25 75],1);
medianexkurt  = median(exkurtV1,1);
iqrexkurt     = prctile(exkurtV1,[25 75],1);
medianEntropy = median(entropyV1,1);
iqrEntropy    = prctile(entropyV1,[25 75],1);
medianJSD     = median(JSD,1);
iqrJSD        = prctile(JSD,[25 75],1);

summaryStats.meanV1RP   = meanV1RP;
summaryStats.stdV1RP    = stdV1RP;
summaryStats.meanglV1RP = meanglV1RP;
summaryStats.stdglV1RP  = stdglV1RP;

summaryStats.meanVar     = medianVar;
summaryStats.stdVar      = iqrVar;
summaryStats.meanexkurt  = medianexkurt;
summaryStats.stdexkurt   = iqrexkurt;
summaryStats.meanEntropy = medianEntropy;
summaryStats.stdEntropy  = iqrEntropy;
summaryStats.meanJSD     = medianJSD;
summaryStats.stdJSD      = iqrJSD;

%% Save reparameterized data

save([topDir,'reparameterizedFits.mat'],'V1_reparams','v1_RP','glV1_RP','summaryStats');


%% Plot each reparameterized V1 distribution against MT (+ Laplacian fits)
titles = {'pedestal','amplitude','env. mean','env. width','frequency','phase'};

f1 = figure; 
f1.Position = [100 100 1500 1000];

for ii = 1:6
    subplot(2,3,ii);
    hold on;

    for rs = 1:numResamps
        plot(x, squeeze(v1_RP(rs,ii,:)),'color',ColorIt('b')+[0.2 0.2 0.2],'linewidth',1);
    end
    pl(1) = plot(x, meanV1RP(ii,:),'color',ColorIt('b'),'linewidth',2.5);
    pl(2) = plot(x, mt,'color',ColorIt('r'),'linewidth',2.5);
    
    glV1F = exp( -(abs(x)/meanglV1RP(ii,1)).^meanglV1RP(ii,2) );
    glV1F = glV1F / sum(glV1F);
    plot( x, glV1F, '--','color',ColorIt('b'),'linewidth',2 );

    glmtF = exp( -(abs(x)/glMT(1)).^glMT(2) );
    glmtF = glmtF / sum(glmtF);
    plot( x, glmtF, '--','color',ColorIt('r'),'linewidth',2 );

    legend(pl,{'V1','MT'});
    set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'xtick',-2:0.5:2);
    xlabel('Disparity');
    ylabel('Normalized Fisher information');
    title(titles{ii});

end


%% Compare Laplace parameters between reparameterizations
statTitles = {'Variance','Excess Kurtosis','Entropy'};
allStats   = {varV1,exkurtV1,entropyV1};
allMeds    = {medianVar,medianexkurt,medianEntropy};
allIQRs    = {iqrVar,iqrexkurt,iqrEntropy};
allMT      = {varMT,exkurtMT,entropyMT};

f2 = figure; 
f2.Position = [400 100 1500 460];

for jj = 1:3

    subplot(1,3,jj);
    hold on;

    thisStat = allStats{jj};
    thisMed = allMeds{jj};
    thisIQR  = allIQRs{jj};

    for ii = 1:6

        swarmchart(ii*ones(size(thisStat,1),1),thisStat(:,ii),50,[0.6 0.6 0.6]);
        scatter(ii,thisMed(ii),100,'k','filled');
        errorbar(ii,thisMed(ii),thisMed(ii)-thisIQR(1,ii),thisIQR(2,ii)-thisMed(ii),'k','linewidth',2);

    end

    plot([0 7],allMT{jj}*[1 1],'--k','linewidth',1.5);

    set(gca,'plotboxaspectratio',[1 1 1],'xlim',[0 7],'xtick',1:6,'xticklabels',titles,'yscale','log');
    xtickangle(45);
    title(statTitles{jj});

end


%% Compare JS divergence between V1 reparameterizations and MT FI distribution

f3 = figure; 
f3.Position = [800 100 900 600];
hold on;

for ii = 1:6
    
    swarmchart(ii*ones(size(JSD,1),1),JSD(:,ii),50,[0.6 0.6 0.6]);
    scatter(ii,medianJSD(ii),100,'k','filled');
    errorbar(ii,medianJSD(ii),medianJSD(ii)-iqrJSD(1,ii),iqrJSD(2,ii)-medianJSD(ii),'k','linewidth',2);

end

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[0 7],'xtick',1:6,'xticklabels',titles,'yscale','log');
xtickangle(45);
ylabel('JSD');
title('V1-MT Population FI Divergence');

%% Save plots

if ~exist(saveDir)
    mkdir(saveDir)
end

saveas


