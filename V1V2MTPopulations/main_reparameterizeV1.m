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


%% Resample V1 parameter fits from the MT distribution + downsample

% Define parameter names
titles = {'pedestal','amplitude','env. mean','env. width','frequency','phase'};

% Define function for Gabor tuning curve
TC = @(x,p) p(1) + p(2)*exp( -(x-p(3)).^2 / (2*p(4)^2) ) .* cos(2*pi*p(5)*(x-p(3))+p(6));

% Define function for derivative of Gabor tuning curve
dTC = @(x,p) -p(2)*exp(-0.5*((x - p(3))/p(4)).^2 ).* ...
             ( (x - p(3))/(p(4)^2).*cos(2*pi*p(5)*(x - p(3)) + p(6))...
               + (2*pi*p(5))*sin(2*pi*p(5)*(x - p(3)) + p(6)) );

% Clip negative FIs
V1.FI(V1.FI < 0) = 0;
MT.FI(MT.FI < 0) = 0;

% Select only final subsample
V1.P  = V1.P(eccInds.V1,:);
V1.FI = V1.FI(eccInds.V1,:);
MT.P  = MT.P(eccInds.MT,:);
MT.FI = MT.FI(eccInds.MT,:);

% Define disparity support corresponding to each column in FI matrices
xFull = [-2 : 0.01 : 2];

% Define disparity support that matches support from scene statistics
res = 52;
lb  = -2;
ub  = 2;
edges_disp = linspace(lb,ub,res);
cntr_disp  = edges_disp(1:end-1) + diff(edges_disp(1:2))/2;

% Find indices of x that intersect with x
[~,ia,~] = intersect(round(xFull,2),round(cntr_disp,2));

% Downsample MT FI distributions by grabbing just those indices
x     = xFull(ia);
MT.FI = MT.FI(:,ia);

numResamps = 500;
% numResamps = 250;

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

    V1_offCntl   = permuteFits(V1,V2,MT,'V1','V1','offset');
    V1_ampCntl   = permuteFits(V1,V2,MT,'V1','V1','amplitude');
    V1_envMCntl  = permuteFits(V1,V2,MT,'V1','V1','envMean');
    V1_envWCntl  = permuteFits(V1,V2,MT,'V1','V1','envWidth');
    V1_freqCntl  = permuteFits(V1,V2,MT,'V1','V1','frequency');
    V1_phaseCntl = permuteFits(V1,V2,MT,'V1','V1','phase');

    V1_reparams     = {V1_off.P, V1_amp.P, V1_envM.P, V1_envW.P, V1_freq.P, V1_phase.P};
    V1_reparamsCntl = {V1_offCntl.P, V1_ampCntl.P, V1_envMCntl.P, V1_envWCntl.P, V1_freqCntl.P, V1_phaseCntl.P};

    % Recalulate individual cell FIs under reparameterization

    % For each parameter
    for ii = 1:6

        numCells = size(V1_reparams{ii},1);

        % For each cell
        for jj = 1:numCells

            p  = V1_reparams{ii}(jj,:);
            pC = V1_reparamsCntl{ii}(jj,:);

            % Calculate Fisher information
            fi  = dTC(xFull,p).^2 ./ TC(xFull,p);
            fiC = dTC(xFull,pC).^2 ./ TC(xFull,pC);

            % Remove FI values less than 0
            lt0Cnt_fi  = sum(fi<0);
            lt0Cnt_fiC = sum(fiC<0);

            fi(fi<0)   = 0;
            fiC(fiC<0) = 0;

            % Downsample
            fi  = fi(ia);
            fiC = fiC(ia);

            eq0Cnt_fi  = sum(fi==0);
            eq0Cnt_fiC = sum(fiC==0);

            % Collect into cell array
            V1_reparamsFI{ii,rs}(jj,:)  = fi;
            V1_reparamsFIC{ii,rs}(jj,:) = fiC;

            V1_rplt0{ii,rs}(jj)  = lt0Cnt_fi;
            V1_rplt0C{ii,rs}(jj) = lt0Cnt_fiC;

            V1_rpDSeq0{ii,rs}(jj)  = eq0Cnt_fi;
            V1_rpDSeq0C{ii,rs}(jj) = eq0Cnt_fiC;

        end

    end
end
end

%% Figure out how bad we messed up the TCs by making non-zero firing rates
perclt0Elems = cellfun(@(x) sum(x),V1_rplt0)/(numCells*(res-1)); % number of <0 elements for each bootstrap/param
perclt0Cells = cellfun(@(x) sum(x~=0),V1_rplt0)/numCells;        % number of cells with <0 elements for each 

meanlt0Elems = mean(perclt0Elems,2);
stdlt0Elems  = std(perclt0Elems,[],2);

meanlt0Cells = mean(perclt0Cells,2);
stdlt0Cells  = std(perclt0Cells,[],2);

numeq0ElemsDS = cellfun(@(x) sum(x),V1_rpDSeq0)/(numCells*(res-1));
numeq0CellsDS = cellfun(@(x) sum(x~=0),V1_rpDSeq0)/numCells;

meaneq0ElemsDS = mean(numeq0ElemsDS,2);
stdeq0ElemsDS  = std(numeq0ElemsDS,[],2);

meaneq0CellsDS = mean(numeq0CellsDS,2);
stdeq0CellsDS  = std(numeq0CellsDS,[],2);

f0 = figure;
f0.Position = [100 100 900 900];

subplot(2,2,1);
hold on;
scatter(1:6,meanlt0Elems,100,'k','filled');
errorbar(meanlt0Elems,stdlt0Elems,'linestyle','none','color','k','linewidth',3);
set(gca,'plotboxaspectratio',[1 1 1],'xtick',1:6,'xticklabels',titles,'xlim',[0.5 6.5],'ylim',[0 1]);
title('% of TC elements < 0');

subplot(2,2,2);
hold on;
scatter(1:6,meaneq0ElemsDS,100,'k','filled');
errorbar(meaneq0ElemsDS,stdeq0ElemsDS,'linestyle','none','color','k','linewidth',3);
set(gca,'plotboxaspectratio',[1 1 1],'xtick',1:6,'xticklabels',titles,'xlim',[0.5 6.5],'ylim',[0 1]);
title('% of Downsampled TC elements = 0');

subplot(2,2,3);
hold on;
scatter(1:6,meanlt0Cells,100,'k','filled');
errorbar(meanlt0Cells,stdlt0Cells,'linestyle','none','color','k','linewidth',3);
set(gca,'plotboxaspectratio',[1 1 1],'xtick',1:6,'xticklabels',titles,'xlim',[0.5 6.5],'ylim',[0 1]);
title('% of TC cells >=1 element FI < 0');

subplot(2,2,4);
hold on;
scatter(1:6,meaneq0CellsDS,100,'k','filled');
errorbar(meaneq0CellsDS,stdeq0CellsDS,'linestyle','none','color','k','linewidth',3);
set(gca,'plotboxaspectratio',[1 1 1],'xtick',1:6,'xticklabels',titles,'xlim',[0.5 6.5],'ylim',[0 1]);
title('% of Downsampled TC cells >=1 element FI = 0');


%% Calculate FI for each reparameterized V1 distribution & MT (+ Laplacian fits)

% Get unnormalized population FI
mt         = sum(MT.FI); 

% Get AUC for pop FI + normalize by cell count
mtUN       = sum(mt)/size(MT.FI,1);

% Normalize pop FI only by cell count
mtcellNorm = mt / size(MT.FI,1);
% Normalize pop FI by AUC
mt         = mt / sum(mt);

% % Fit a generalized Laplacian to the population FI
% glV1 = fitGeneralizedLaplacian(x,v1);
% glV2 = fitGeneralizedLaplacian(x,v2);
% glMT = fitGeneralizedLaplacian(x,mt);

% Get summary stats for mt
% varMT     = (glMT(1).^2 .*gamma(3./glMT(2)) ) ./gamma(1./glMT(2));
% exkurtMT  = ( (gamma(5./glMT(2)).*gamma(1./glMT(2))) ./ (gamma(3./glMT(2)).^2) ) - 3;
% entropyMT = 1./glMT(2) - log( glMT(2)./ (2*glMT(1).*gamma(1./glMT(2))) );

if saveOn
% Downsample, sum each reparameterized V1 dist, & fit with Laplacian
JSD         = nan(numResamps,6);
JSDcellNorm = nan(numResamps,6);
JSDC        = nan(numResamps,6);

for ii = 1:6
    
    theseRSPopFIdist         = nan(numResamps,res-1);
    theseRSPopFIdistCellNorm = nan(numResamps,res-1);
    theseCntlRSPopFIdist     = nan(numResamps,res-1);

    for rs = 1:numResamps

        % Collect downsampled FI distributions for each cell in this
        % bootstrapped sample from MT distribution
        theseRSFIdists     = V1_reparamsFI{ii,rs};
        theseCntlRSFIdists = V1_reparamsFIC{ii,rs};

        sampleCellCnt      = size(theseRSFIdists,1);

        % Sum over all the cells
        thisRSPopFIdist         = sum(theseRSFIdists);

        % Sum the area under this un-normalized population FI and divide by
        % total cell count
        theseRSTFI(rs,ii)       = sum(thisRSPopFIdist) / sampleCellCnt;

        % Normalize the population FI only by number of cells in sample
        thisRSPopFIdistCellNorm = thisRSPopFIdist / sampleCellCnt;

        % Normalize the population FI so AUC sums to 1
        thisRSPopFIdist         = thisRSPopFIdist / sum(thisRSPopFIdist);

        % Repeat above steps for the bootstraps sampled from V1
        % distribution
        thisCntlRSPopFIdist         = sum(theseCntlRSFIdists);
        theseCntlRSTFI(rs,ii)       = sum(thisCntlRSPopFIdist) / sampleCellCnt;
        thisCntlRSPopFIdistCellNorm = thisCntlRSPopFIdist / sampleCellCnt;
        thisCntlRSPopFIdist         = thisCntlRSPopFIdist / sum(thisCntlRSPopFIdist);

%         glV1_RP(rs,ii,:)  = fitGeneralizedLaplacian(x,squeeze(v1_RP(rs,ii,:))');
%         glV1_RPC(rs,ii,:) = fitGeneralizedLaplacian(x,squeeze(v1_RPC(rs,ii,:))');

%         varV1(rs,ii)     = (glV1_RP(rs,ii,1).^2 .*gamma(3./glV1_RP(rs,ii,2)) ) ./gamma(1./glV1_RP(rs,ii,2));
%         exkurtV1(rs,ii)  = ( (gamma(5./glV1_RP(rs,ii,2)).*gamma(1./glV1_RP(rs,ii,2))) ./ (gamma(3./glV1_RP(rs,ii,2)).^2) ) - 3;
%         entropyV1(rs,ii) = 1./glV1_RP(rs,ii,2) - log( glV1_RP(rs,ii,2)./ (2*glV1_RP(rs,ii,1).*gamma(1./glV1_RP(rs,ii,2))) );

        % Calculate JSD for resampling from MT
        JSD(rs,ii)         = getJSDiv(squeeze(thisRSPopFIdist),mt);
        JSDcellNorm(rs,ii) = getJSDiv(squeeze(thisRSPopFIdistCellNorm),mtcellNorm);

%         varV1C(rs,ii)     = (glV1_RPC(rs,ii,1).^2 .*gamma(3./glV1_RPC(rs,ii,2)) ) ./gamma(1./glV1_RPC(rs,ii,2));
%         exkurtV1C(rs,ii)  = ( (gamma(5./glV1_RPC(rs,ii,2)).*gamma(1./glV1_RPC(rs,ii,2))) ./ (gamma(3./glV1_RPC(rs,ii,2)).^2) ) - 3;
%         entropyV1C(rs,ii) = 1./glV1_RPC(rs,ii,2) - log( glV1_RPC(rs,ii,2)./ (2*glV1_RPC(rs,ii,1).*gamma(1./glV1_RPC(rs,ii,2))) );

        % Calculate JSD for resampling from V1
        JSDC(rs,ii)       = getJSDiv(squeeze(thisCntlRSPopFIdist),mt);

        % Collect into matrices for stats
        theseRSPopFIdist(rs,:)         = thisRSPopFIdist;
        theseRSPopFIdistCellNorm(rs,:) = thisRSPopFIdistCellNorm;
        theseCntlRSPopFIdist(rs,:)     = thisCntlRSPopFIdist;

    end

    % Define alpha (set here to find 95% CIs)
    alpha = 0.05;

    % Loop over elements in disparity support
    for jj = 1:(res-1)

        if 0

        % Find kernel-smoothed density of FI values for each element of disparity support
        [thisKSD,theseVals]                 = ksdensity(theseRSPopFIdist(:,jj),...
                                                'Support','positive','BoundaryCorrection','reflection');
        [thisKSDcellNorm,theseValscellNorm] = ksdensity(theseRSPopFIdistCellNorm(:,jj),...
                                                'Support','positive','BoundaryCorrection','reflection');
        [thisKSDC,theseValsC]               = ksdensity(theseCntlRSPopFIdist(:,jj),...
                                                'Support','positive','BoundaryCorrection','reflection');

        % Zero out values 

        % Convert density to normalized cumulative distribution
        thisCumDens          = cumsum(thisKSD/sum(thisKSD));
        thisCumDenscellNorm  = cumsum(thisKSDcellNorm/sum(thisKSDcellNorm));
        thisCumDensC         = cumsum(thisKSDC/sum(thisKSDC));

        % Find indices of CDF closest to 2.5% and 97.5% (or whatever alpha is set to)
        [~,LEBind]         = min(abs(thisCumDens - alpha/2));
        [~,UEBind]         = min(abs(thisCumDens - (1 - alpha/2)));

        [~,LEBindcellNorm] = min(abs(thisCumDenscellNorm - alpha/2));
        [~,UEBindcellNorm] = min(abs(thisCumDenscellNorm - (1 - alpha/2)));

        [~,LEBindC]        = min(abs(thisCumDensC - alpha/2));
        [~,UEBindC]        = min(abs(thisCumDensC - (1 - alpha/2)));

        % Find median FI for this support element and index into KSD support 
        % (Sampling from MT)
        v1_med(ii,jj) = median(squeeze(theseRSPopFIdist(:,jj)),'omitnan');
        v1_LEB(ii,jj) = theseVals(LEBind);
        v1_UEB(ii,jj) = theseVals(UEBind);

        v1_medcellNorm(ii,jj) = median(squeeze(theseRSPopFIdistCellNorm(:,jj)),'omitnan');
        v1_LEBcellNorm(ii,jj) = theseValscellNorm(LEBindcellNorm);
        v1_UEBcellNorm(ii,jj) = theseValscellNorm(UEBindcellNorm);

        % (Sampling from V1)
        v1_medC(ii,jj) = median(squeeze(theseCntlRSPopFIdist(:,jj)),'omitnan');
        v1_LEBC(ii,jj) = theseVals(LEBindC);
        v1_UEBC(ii,jj) = theseVals(UEBindC);

        else

        % Find median FI for this support element and index into KSD support 
        % (Sampling from MT)
        v1_med(ii,jj) = median(squeeze(theseRSPopFIdist(:,jj)),'omitnan');
        v1_LEB(ii,jj) = prctile(squeeze(theseRSPopFIdist(:,jj)),[25]);
        v1_UEB(ii,jj) = prctile(squeeze(theseRSPopFIdist(:,jj)),[75]);

        v1_medcellNorm(ii,jj) = median(squeeze(theseRSPopFIdistCellNorm(:,jj)),'omitnan');
        v1_LEBcellNorm(ii,jj) = prctile(squeeze(theseRSPopFIdistCellNorm(:,jj)),[25]);
        v1_UEBcellNorm(ii,jj) = prctile(squeeze(theseRSPopFIdistCellNorm(:,jj)),[75]);

        % (Sampling from V1)
        v1_medC(ii,jj) = median(squeeze(theseCntlRSPopFIdist(:,jj)),'omitnan');
        v1_LEBC(ii,jj) = prctile(squeeze(theseCntlRSPopFIdist(:,jj)),[25]);
        v1_UEBC(ii,jj) = prctile(squeeze(theseCntlRSPopFIdist(:,jj)),[75]);

        end

    end

end


%% Calculate metastatistics

% Sampled from MT
medianJSD     = median(JSD,1);
iqrJSD        = prctile(JSD,[25 75],1);
medianTFI     = median(theseRSTFI,1);
iqrTFI        = prctile(theseRSTFI,[25 75],1);
% medianVar     = median(varV1,1);
% iqrVar        = prctile(varV1,[25 75],1);
% medianexkurt  = median(exkurtV1,1);
% iqrexkurt     = prctile(exkurtV1,[25 75],1);
% medianEntropy = median(entropyV1,1);
% iqrEntropy    = prctile(entropyV1,[25 75],1);

medianJSDcellNorm     = median(JSDcellNorm,1);
iqrJSDcellNorm        = prctile(JSDcellNorm,[25 75],1);

% Sampled from V1
medianJSDC     = median(JSDC,1);
iqrJSDC        = prctile(JSDC,[25 75],1);
medianTFIC     = median(theseCntlRSTFI,1);
iqrTFIC        = prctile(theseCntlRSTFI,[25 75],1);
% medianVarC     = median(varV1C,1);
% iqrVarC        = prctile(varV1C,[25 75],1);
% medianexkurtC  = median(exkurtV1C,1);
% iqrexkurtC     = prctile(exkurtV1C,[25 75],1);
% medianEntropyC = median(entropyV1C,1);
% iqrEntropyC    = prctile(entropyV1C,[25 75],1);

% Omnibus test for differences in median JSD between parameters
values   = JSD(:);
groupVec = repelem(1:6,numResamps);

[kwTestp,thisTable] = kruskalwallis(values,groupVec,'off');
kwTestStats         = [thisTable{2,3} thisTable{2,5}];

% Pairwise tests for significant differences
numPars = 6;

[a,b]   = meshgrid(1:(numPars-1),1:(numPars-1));
c       = a + b ;
d       = c(c <= numPars);

pairs   = repelem(1:numPars-1,fliplr(1:numPars-1))';
pairMat = [pairs d];

for ii = 1:numel(pairs)

    setA = JSD(:,pairMat(ii,1));
    setB = JSD(:,pairMat(ii,2));

    theseVals                    = [setA setB];
    [rsTestp(ii,1),~,thisStruct] = ranksum(theseVals(:,1),theseVals(:,2));
    rsTestStats(ii,1)            = thisStruct.zval;

end

% Collect stats
summaryStats.kwTestp     = kwTestp;
summaryStats.kwTestStats = kwTestStats;
summaryStats.rsTestp     = [pairMat rsTestp];
summaryStats.rsTestStats = [pairMat rsTestStats];


%% Collect stats into structure
summaryStats.medV1RP = v1_med;
summaryStats.v1_LEB  = v1_LEB;
summaryStats.v1_UEB  = v1_UEB;

summaryStats.medV1RPcellNorm = v1_medcellNorm;
summaryStats.v1_LEBcellNorm  = v1_LEBcellNorm;
summaryStats.v1_UEBcellNorm  = v1_UEBcellNorm;

summaryStats.medV1RPC = v1_medC;
summaryStats.v1_LEBC  = v1_LEBC;
summaryStats.v1_UEBC  = v1_UEBC;

% summaryStats.medianVar     = medianVar;
% summaryStats.iqrVar        = iqrVar;
% summaryStats.medianexkurt  = medianexkurt;
% summaryStats.iqrexkurt     = iqrexkurt;
% summaryStats.medianEntropy = medianEntropy;
% summaryStats.iqrEntropy    = iqrEntropy;
summaryStats.medianJSD     = medianJSD;
summaryStats.iqrJSD        = iqrJSD;

summaryStats.medianJSDcellNorm     = medianJSDcellNorm;
summaryStats.iqrJSDcellNorm        = iqrJSDcellNorm;

% summaryStats.medianVarC     = medianVarC;
% summaryStats.iqrVarC        = iqrVarC;
% summaryStats.medianexkurtC  = medianexkurtC;
% summaryStats.iqrexkurtC     = iqrexkurtC;
% summaryStats.medianEntropyC = medianEntropyC;
% summaryStats.iqrEntropyC    = iqrEntropyC;
summaryStats.medianJSDC     = medianJSDC;
summaryStats.iqrJSDC        = iqrJSDC;

end


%% Save reparameterized data
if saveOn

    save([dataDir,'reparameterizedFits.mat'],'summaryStats');

else

    load([dataDir,'reparameterizedFits.mat']);
    
    v1_med = summaryStats.medV1RP;
    v1_LEB = summaryStats.v1_LEB;
    v1_UEB = summaryStats.v1_UEB;

    v1_med = summaryStats.medV1RPcellNorm;
    v1_LEB = summaryStats.v1_LEBcellNorm;
    v1_UEB = summaryStats.v1_UEBcellNorm;

    v1_medC = summaryStats.medV1RPC;
    v1_LEBC = summaryStats.v1_LEBC;
    v1_UEBC = summaryStats.v1_UEBC;

    medianJSD = summaryStats.medianJSD;
    iqrJSD    = summaryStats.iqrJSD;

    medianJSDcellNorm = summaryStats.medianJSDcellNorm;
    iqrJSDcellNorm    = summaryStats.iqrJSDcellNorm;
    
    medianJSDC     = summaryStats.medianJSDC;
    iqrJSDC        = summaryStats.iqrJSDC;

end


%% Plot each reparameterized V1 distribution against MT (+ Laplacian fits)
f1 = figure; 
f1.Position = [100 100 1500 1000];

for ii = 1:6
    subplot(2,3,ii);
    hold on;

    % Plot shadeplot
    pl(1) = shadeplot(v1_med(ii,:),[v1_LEB(ii,:); v1_UEB(ii,:)],x,[0 0 0],0.4,2.5);
    pl(2) = plot(x, mt,'color',ColorIt('r'),'linewidth',2.5);

    if ii == 1
    legend(pl,{'V1','MT'},'location','northwest');
    end

    set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'xtick',-2:1:2,'ylim',[0 0.2],'fontsize',15);
    xlabel('Disparity');
    ylabel('Normalized Fisher information');
    title(titles{ii},'AUC normalized');

end

f2 = figure; 
f2.Position = [1100 100 1500 1000];

for ii = 1:6
    subplot(2,3,ii);
    hold on;

    % Since we want to plot on log scale here due to huge differences in
    % linear scale between pars, replace lower CI zero values with eps
    v1_LEBcellNorm(v1_LEBcellNorm==0) = eps;

    % Plot shadeplot
    pl(1) = shadeplot(v1_medcellNorm(ii,:),[v1_LEBcellNorm(ii,:); v1_UEBcellNorm(ii,:)],x,[0 0 0],0.4,2.5);
    pl(2) = plot(x, mtcellNorm,'color',ColorIt('r'),'linewidth',2.5);
%     plot(x, v1_medcellNorm(ii,:),'color',[0 0 0],'linewidth',2.5);
%     plot(x, v1_LEBcellNorm(ii,:),'color',0.5*[1 1 1],'linewidth',2.5);
%     plot(x, v1_UEBcellNorm(ii,:),'color',0.5*[1 1 1],'linewidth',2.5);

    if ii == 1
    legend(pl,{'V1','MT'},'location','northwest');
    end

    set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'xtick',-2:1:2,'ylim',[1E-1 1.9E4],'yscale','log','fontsize',15);
    xlabel('Disparity');
    ylabel('Normalized Fisher information');
    title(titles{ii},' Cell count normalized');

end


%% Compare Laplace parameters between reparameterizations
% statTitles = {'Variance','Excess Kurtosis','Entropy'};
% allStats   = {varV1,exkurtV1,entropyV1};
% allMeds    = {medianVar,medianexkurt,medianEntropy};
% allIQRs    = {iqrVar,iqrexkurt,iqrEntropy};
% allMT      = {varMT,exkurtMT,entropyMT};
% 
% f2 = figure; 
% f2.Position = [400 100 1500 460];
% 
% for jj = 1:3
% 
%     subplot(1,3,jj);
%     hold on;
% 
%     thisStat = allStats{jj};
%     thisMed = allMeds{jj};
%     thisIQR  = allIQRs{jj};
% 
%     for ii = 1:6
% 
%         swarmchart(ii*ones(size(thisStat,1),1),thisStat(:,ii),50,[0.6 0.6 0.6],'XJitterWidth',0.7);
%         p(1) = scatter(ii,thisMed(ii),100,'k','filled');
%         errorbar(ii,thisMed(ii),thisMed(ii)-thisIQR(1,ii),thisIQR(2,ii)-thisMed(ii),'k','linewidth',2);
% 
%     end
% 
%     pl(2) = plot([0 7],allMT{jj}*[1 1],'--r','linewidth',2.5);
% 
%     set(gca,'plotboxaspectratio',[1 1 1],'xlim',[0 7],'xtick',1:6,'xticklabels',titles,'yscale','log','fontsize',15);
%     xtickangle(45);
%     title(statTitles{jj});
%     
%     if jj == 1
%         legend(pl,{'V1','MT'});
%     end
% 
% end


%% Compare JS divergence between V1 reparameterizations and MT FI distribution

f3 = figure; 
f3.Position = [800 100 625 600];
hold on;

for ii = 1:6

    [p(ii),~,thisStruct] = ranksum(JSD(:,ii),JSDC(:,ii));
    statsJSD(ii,1)       = thisStruct.zval;
    
    if p(ii) < 0.05
        scatter(ii,0.9,100,'k','*');
    end

    swarmchart(ii*ones(size(JSD,1),1)+0.125,JSD(:,ii),50,[.3 .3 .3] + 0.7*ColorIt('o'),'XJitterWidth',0.4);
    swarmchart(ii*ones(size(JSDC,1),1)-0.125,JSDC(:,ii),50,[0.6 0.6 0.6],'XJitterWidth',0.4);

    pl(1) = scatter(ii+0.125,medianJSD(ii),100,ColorIt('o'),'filled');
    e = errorbar(ii+0.125,medianJSD(ii),medianJSD(ii)-iqrJSD(1,ii),iqrJSD(2,ii)-medianJSD(ii),'k','linewidth',2);
    e.Color = ColorIt('o');

    pl(2) = scatter(ii-0.125,medianJSDC(ii),100,'k','filled');
    errorbar(ii-0.125,medianJSDC(ii),medianJSDC(ii)-iqrJSDC(1,ii),iqrJSDC(2,ii)-medianJSDC(ii),'k','linewidth',2);
    
end

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[0 7],'xtick',1:6,'xticklabels',titles,'yscale','log','fontsize',15);
xtickangle(45);
ylabel('JSD');
title('V1-MT Population FI Divergence');
legend(pl,{'MT','V1 (shuffled)'},'Location','southeast');


%% Plot AUC ("total fisher info") for the true and reparameterized distributions
f4 = figure; 
f4.Position = [800 800 650 600];
hold on;

for ii = 1:6
    
    [p2(ii),~,thisStruct] = ranksum(theseRSTFI(:,ii),theseCntlRSTFI(:,ii));
    statsTFI(ii,1)        = thisStruct.zval;

    if p2(ii) < 0.05
        scatter(ii,4*1e7,100,'k','*');
    end

    swarmchart(ii*ones(size(theseRSTFI,1),1)+0.125,theseRSTFI(:,ii),50,[.3 .3 .3] + 0.7*ColorIt('o'),'XJitterWidth',0.4);
    swarmchart(ii*ones(size(theseCntlRSTFI,1),1)-0.125,theseCntlRSTFI(:,ii),50,[0.6 0.6 0.6],'XJitterWidth',0.4);

    pl(1) = scatter(ii+0.125,medianTFI(ii),100,ColorIt('o'),'filled');
    e = errorbar(ii+0.125,medianTFI(ii),medianTFI(ii)-iqrTFI(1,ii),iqrTFI(2,ii)-medianTFI(ii),'o','linewidth',2);
    e.Color = ColorIt('o');

    pl(2) = scatter(ii-0.125,medianTFIC(ii),100,'k','filled');
    errorbar(ii-0.125,medianTFIC(ii),medianTFIC(ii)-iqrTFIC(1,ii),iqrTFIC(2,ii)-medianTFIC(ii),'k','linewidth',2);

    plot([0 7],mtUN*[1 1],'--r','linewidth',2.5);

end

set(gca,'plotboxaspectratio',[1 1 1],'xlim',[0 7],'xtick',1:6,'xticklabels',titles,'yscale','log',...
    'ylim',5*[1e2 1e6],'fontsize',15);
xtickangle(45);
ylabel('AUC');
title('V1-MT Total Population FI');
legend(pl,{'MT','V1 (shuffled)'},'Location','northeast');


%% Save plots

summaryStats.pJSD     = p;
summaryStats.pTFI     = p2;
summaryStats.statsJSD = statsJSD;
summaryStats.statsTFI = statsTFI;

if saveOn
    
    save([dataDir,'reparameterizedFits.mat'],'summaryStats');

    if ~exist(saveDir)
        mkdir(saveDir)
    end

    saveas(f1,[saveDir,'reparameterizationBootstraps.svg']);
    saveas(f2,[saveDir,'reparameterizationBootstrapsCellNumNorm.svg']);
    saveas(f3,[saveDir,'reparJSD.svg']);
    saveas(f4,[saveDir,'reparTFI.svg']);
    
end


