% Assess whether we can account for FI differences between V1 and MT based
% on substituting a single parameter distribution from one to the other

clear all
close all

%% Load in saved Gabor fits

% Define directory structure
splPath = regexp(which('main_reparameterizeV1.m'),filesep,'split');
topDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
dataDir = [topDir,'analysisFiles',filesep];
saveDir = [topDir,'plots/reparameterization/'];

addpath(genpath([topDir,'helper_functions']));

% Load in neuronal data and indices of neurons in final subsample
load([dataDir,'resultsALL_final.mat']);
load([dataDir,'eccInds.mat']);

% Toggle save on/off
saveOn = 1;


%% Resample V1 parameter fits from the MT distribution

% Clip negative FIs
V1.FI(V1.FI < 0) = 0;
V2.FI(V2.FI < 0) = 0;
MT.FI(MT.FI < 0) = 0;

% Select only final subsample
V1.P  = V1.P(eccInds.V1,:);
V1.FI = V1.FI(eccInds.V1,:);
V2.P  = V2.P(eccInds.V2,:);
V2.FI = V2.FI(eccInds.V2,:);
MT.P  = MT.P(eccInds.MT,:);
MT.FI = MT.FI(eccInds.MT,:);

numResamps = 500;

if saveOn
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
v1   = sum(V1.FI); 
v1UN = sum(v1);
v1   = v1 / sum(v1);

v2   = sum(V2.FI); 
v2UN = sum(v2);
v2   = v2 / sum(v2);

mt   = sum(MT.FI); 
mtUN = sum(mt);
mt   = mt / sum(mt);

% Fit a generalized Laplacian to the population FI
glV1 = fitGeneralizedLaplacian(x,v1);
glV2 = fitGeneralizedLaplacian(x,v2);
glMT = fitGeneralizedLaplacian(x,mt);

% Get summary stats for mt
varMT     = (glMT(1).^2 .*gamma(3./glMT(2)) ) ./gamma(1./glMT(2));
exkurtMT  = ( (gamma(5./glMT(2)).*gamma(1./glMT(2))) ./ (gamma(3./glMT(2)).^2) ) - 3;
entropyMT = 1./glMT(2) - log( glMT(2)./ (2*glMT(1).*gamma(1./glMT(2))) );

if saveOn
% Downsample, sum each reparameterized V1 dist, & fit with Laplacian
for ii = 1:6
    for rs = 1:numResamps

        V1_reparamsFIds{ii,rs} = V1_reparamsFI{ii,rs}(:,ia);

        v1_RP(rs,ii,:)  = sum(V1_reparamsFIds{ii,rs});
        v1_TFI(rs,ii,:) = sum(v1_RP(rs,ii,:));
        v1_RP(rs,ii,:)  = v1_RP(rs,ii,:) / sum(v1_RP(rs,ii,:));

        glV1_RP(rs,ii,:) = fitGeneralizedLaplacian(x,squeeze(v1_RP(rs,ii,:))');

        varV1(rs,ii)     = (glV1_RP(rs,ii,1).^2 .*gamma(3./glV1_RP(rs,ii,2)) ) ./gamma(1./glV1_RP(rs,ii,2));
        exkurtV1(rs,ii)  = ( (gamma(5./glV1_RP(rs,ii,2)).*gamma(1./glV1_RP(rs,ii,2))) ./ (gamma(3./glV1_RP(rs,ii,2)).^2) ) - 3;
        entropyV1(rs,ii) = 1./glV1_RP(rs,ii,2) - log( glV1_RP(rs,ii,2)./ (2*glV1_RP(rs,ii,1).*gamma(1./glV1_RP(rs,ii,2))) );

        JSD(rs,ii)       = getJSDiv(squeeze(v1_RP(rs,ii,:))',mt);

    end

    % Find 95% CIs with 1D KSD
    alpha = 0.05;

    for jj = 1:size(v1_RP,3)
        [thisKSD,theseVals] = ksdensity(v1_RP(:,ii,jj));

        thisCumDens = cumsum(thisKSD/sum(thisKSD));

        [~,LCIind] = min(abs(thisCumDens - alpha/2));
        [~,UCIind] = min(abs(thisCumDens - (1 - alpha/2)));

        v1_med(ii,jj) = median(squeeze(v1_RP(:,ii,jj)),'omitnan');
        v1_LCI(ii,jj) = theseVals(LCIind);
        v1_UCI(ii,jj) = theseVals(UCIind);
    end

end


%% Find the mean and std (95% CIs?) of the bootstrapped FIs, Laplace fits and JS divergence

meanV1RP = squeeze(mean(v1_RP,1));
stdV1RP  = squeeze(std(v1_RP,[],1));
meanglV1RP = squeeze(mean(glV1_RP,1));
stdglV1RP  = squeeze(std(glV1_RP,[],1));

summaryStats.meanV1RP   = meanV1RP;
summaryStats.stdV1RP    = stdV1RP;
summaryStats.meanglV1RP = meanglV1RP;
summaryStats.stdglV1RP  = stdglV1RP;
summaryStats.medV1RP    = v1_med;
summaryStats.v1_LCI     = v1_LCI;
summaryStats.v1_UCI     = v1_UCI;

medianVar     = median(varV1,1);
iqrVar        = prctile(varV1,[25 75],1);
medianexkurt  = median(exkurtV1,1);
iqrexkurt     = prctile(exkurtV1,[25 75],1);
medianEntropy = median(entropyV1,1);
iqrEntropy    = prctile(entropyV1,[25 75],1);
medianJSD     = median(JSD,1);
iqrJSD        = prctile(JSD,[25 75],1);
medianTFI     = median(v1_TFI,1);
iqrTFI        = prctile(v1_TFI,[25 75],1);

summaryStats.medianVar     = medianVar;
summaryStats.iqrVar        = iqrVar;
summaryStats.medianexkurt  = medianexkurt;
summaryStats.iqrexkurt     = iqrexkurt;
summaryStats.medianEntropy = medianEntropy;
summaryStats.iqrEntropy    = iqrEntropy;
summaryStats.medianJSD     = medianJSD;
summaryStats.iqrJSD        = iqrJSD;
end


%% Save reparameterized data
if saveOn
    save([dataDir,'reparameterizedFits.mat'],'V1_reparams','v1_RP','glV1_RP','summaryStats');
else
    load([dataDir,'reparameterizedFits.mat']);
    
    v1_med = summaryStats.medV1RP;
    v1_LCI = summaryStats.v1_LCI;
    v1_UCI = summaryStats.v1_UCI;
    meanglV1RP = summaryStats.meanglV1RP;

    medianVar     = summaryStats.medianVar;
    iqrVar        = summaryStats.iqrVar;
    medianexkurt  = summaryStats.medianexkurt;
    iqrexkurt     = summaryStats.iqrexkurt;
    medianEntropy = summaryStats.medianEntropy;
    iqrEntropy    = summaryStats.iqrEntropy;
    medianJSD     = summaryStats.medianJSD;
    iqrJSD        = summaryStats.iqrJSD;
end


%% Plot each reparameterized V1 distribution against MT (+ Laplacian fits)
titles = {'pedestal','amplitude','env. mean','env. width','frequency','phase'};

% Define nicer looking color palette
blueCI  = ColorIt('b');
greenCI = ColorIt('g');
redCI   = ColorIt('r');

f1 = figure; 
f1.Position = [100 100 1500 1000];

for ii = 1:6
    subplot(2,3,ii);
    hold on;

    if 0
        % Plot bootstrapped FI dists
        for rs = 1:numResamps
            plot(x, squeeze(v1_RP(rs,ii,:)),'color',ColorIt('b')+[0.2 0.2 0.2],'linewidth',1);
        end

        % Plot mean/true FI dists
        pl(1) = plot(x, meanV1RP(ii,:),'color',ColorIt('b'),'linewidth',2.5);
        pl(2) = plot(x, mt,'color',ColorIt('r'),'linewidth',2.5);
    else
        % Plot shadeplot
        pl(1) = shadeplot(v1_med(ii,:),[v1_LCI(ii,:); v1_UCI(ii,:)],x,[0 0 0],0.4,2.5);
        pl(2) = plot(x, mt,'color',ColorIt('r'),'linewidth',2.5);
    end


    % Plot mean/true Laplace fits
    glV1F = exp( -(abs(x)/meanglV1RP(ii,1)).^meanglV1RP(ii,2) );
    glV1F = glV1F / sum(glV1F);
    pl(3) = plot( x, glV1F, '--','color','k','linewidth',2 );

    glmtF = exp( -(abs(x)/glMT(1)).^glMT(2) );
    glmtF = glmtF / sum(glmtF);
    pl(4) = plot( x, glmtF, '--','color',ColorIt('r'),'linewidth',2 );

    if ii == 1
    legend(pl,{'V1','MT','V1 Fit','MT Fit'},'location','northwest');
    end

    set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'xtick',-2:1:2,'ylim',[0 0.4],'fontsize',15);
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

        swarmchart(ii*ones(size(thisStat,1),1),thisStat(:,ii),50,[0.6 0.6 0.6],'XJitterWidth',0.7);
        p(1) = scatter(ii,thisMed(ii),100,'k','filled');
        errorbar(ii,thisMed(ii),thisMed(ii)-thisIQR(1,ii),thisIQR(2,ii)-thisMed(ii),'k','linewidth',2);

    end

    pl(2) = plot([0 7],allMT{jj}*[1 1],'--r','linewidth',2.5);

    set(gca,'plotboxaspectratio',[1 1 1],'xlim',[0 7],'xtick',1:6,'xticklabels',titles,'yscale','log','fontsize',15);
    xtickangle(45);
    title(statTitles{jj});
    
    if jj == 1
        legend(pl,{'V1','MT'});
    end

end


%% Compare JS divergence between V1 reparameterizations and MT FI distribution

f3 = figure; 
f3.Position = [800 100 625 600];
hold on;

for ii = 1:6
    
    swarmchart(ii*ones(size(JSD,1),1),JSD(:,ii),50,[0.6 0.6 0.6],'XJitterWidth',0.7);
    scatter(ii,medianJSD(ii),100,'k','filled');
    errorbar(ii,medianJSD(ii),medianJSD(ii)-iqrJSD(1,ii),iqrJSD(2,ii)-medianJSD(ii),'k','linewidth',2);

end

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[0 7],'xtick',1:6,'xticklabels',titles,'yscale','log','fontsize',15);
xtickangle(45);
ylabel('JSD');
title('V1-MT Population FI Divergence');


%% Plot AUC ("total fisher info") for the true and reparameterized distributions
f4 = figure; 
f4.Position = [800 800 650 600];
hold on;

for ii = 1:6
    
    swarmchart(ii*ones(size(v1_TFI,1),1),v1_TFI(:,ii),50,[0.6 0.6 0.6],'XJitterWidth',0.7);
    pl(1) = scatter(ii,medianTFI(ii),100,'k','filled');
    errorbar(ii,medianTFI(ii),medianTFI(ii)-iqrTFI(1,ii),iqrTFI(2,ii)-medianTFI(ii),'k','linewidth',2);

    pl(2) = plot([0 7],mtUN*[1 1],'--r','linewidth',2.5);

    if jj == 1
        legend(p,{'V1','MT'});
    end

end

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[0 7],'xtick',1:6,'xticklabels',titles,'yscale','log',...
    'ylim',[1e5 1e8],'fontsize',15);
xtickangle(45);
ylabel('AUC');
title('V1-MT Total Population FI');


%% Save plots

if saveOn
    if ~exist(saveDir)
        mkdir(saveDir)
    end

    saveas(f1,[saveDir,'reparameterizationBootstraps.svg']);
    saveas(f2,[saveDir,'reparStats.svg']);
    saveas(f3,[saveDir,'reparJSD.svg']);
    saveas(f4,[saveDir,'reparTFI.svg']);
end


