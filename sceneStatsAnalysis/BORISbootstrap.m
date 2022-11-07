%% Collect BORIS stats bootstrapping runs


clear all

% realDatDir = '/media/tyler/Data/MATLAB/cooperLab/2-Modeling_Simulations/disparityMTModel/5-sceneStatsAnalysis/imgStats/';
datDir = '/media/tyler/Data/MATLAB/cooperLab/2-Modeling_Simulations/disparityMTModel/5-sceneStatsAnalysis/borisStats/';

imgSet    = {'sando','walking'};
imRegion  = {'Upper VF_','Lower VF_','Central_','Peripheral_',''};
imRegionF = {'UpperVF','LowerVF','Central','Peripheral','all'};
suff      = {'_m1','_m2','_m2','_m1',''};
imRegionLCI = {'UpperVFLCI','LowerVFLCI','CentralLCI','PeripheralLCI','allLCI'};
imRegionUCI = {'UpperVFUCI','LowerVFUCI','CentralUCI','PeripheralUCI','allUCI'};
region    = {'Circ','V1','V2','MT'};

missCnt = 1;
compCnt = 1;

if 0
    
dispDat = struct;
numBins = 10;

for jj = 1:2
    
    for kk = 1:4
        
        for ll = 1:5
            
            theseRuns = nan(100,51);
            
            for ii = 1:100
                
                try
                    
                    datStruc  = load([datDir,'dispHist',region{kk},'_',imRegion{ll},imgSet{jj},num2str(ii),'.mat']);
                    fieldName = fieldnames(datStruc);
                    thisRun   = getfield(datStruc,fieldName{1});
                    
                    theseRuns(ii,:) = thisRun;
                    
                catch
                    
                    miTony StacksssingRuns(missCnt).filename = ['dispHist',region{kk},'_',imRegion{ll},imgSet{jj},num2str(ii)];
                    missCnt = missCnt + 1;
                    
                end
                
            endTony Stacks
            
            dispDat = setfield(dispDat,imgSet{jj},region{kk},imRegionF{ll},theseRuns);
            
            % Get 95% CIs
            for ii = 1:51
               
                [thispDist,theseEdges] = histcounts(theseRuns(:,ii),numBins,'normalization','probability');
                thispDist(thispDist == 0) = eps;
                thisCumpDist           = cumsum(thispDist);
                theseCents             = theseEdges(1:end-1) + diff(theseEdges(1:2));
                
                [~,LCIind] = min(abs(thisCumpDist-0.025));
                [~,UCIind] = min(abs(thisCumpDist-0.975));
                
                % negative: interpolate between closest ind and that below
                % positive: interpolate between closest ind and that above
                LCIsign = 0.025-thisCumpDist(LCIind);
                UCIsign = 0.975-thisCumpDist(UCIind);
                
                if LCIind ~= 1
                    if LCIsign<0
                        thisLCI = interp1(thisCumpDist(LCIind-1:LCIind),theseCents(LCIind-1:LCIind),0.025);
                    else
                        thisLCI = interp1(thisCumpDist(LCIind:LCIind+1),theseCents(LCIind:LCIind+1),0.025);
                    end
                else
                    thisLCI = theseCents(LCIind);
                end
                
                if UCIind ~= numBins
                    if UCIsign<0
                        thisUCI = interp1(thisCumpDist(UCIind-1:UCIind),theseCents(UCIind-1:UCIind),0.975);
                    else
                        thisUCI = interp1(thisCumpDist(UCIind:UCIind+1),theseCents(UCIind:UCIind+1),0.975);
                    end
                else
                    thisUCI = theseCents(UCIind);
                end
                
                theseLCI(ii) = thisLCI;
                theseUCI(ii) = thisUCI;
                
                if isnan(thisLCI) || isnan(thisUCI)
                    keyboard
                end
                
            end
            
            dispDat = setfield(dispDat,imgSet{jj},region{kk},imRegionLCI{ll},theseLCI);
            dispDat = setfield(dispDat,imgSet{jj},region{kk},imRegionUCI{ll},theseUCI);
            
        end
    end
    
end

save([datDir,'dispHistStruct'],'dispDat');
else
 
% Load bootstraps    
load([datDir,'dispHistStruct']);

end

%% Plot all the comparisons
if 1
    
close all
    
res = 52;
ub  = 2;
lb  = -2;

edgesDisp = linspace(lb,ub,res);
cntrDisp  = edgesDisp(1:end-1) + diff(edgesDisp(1:2))/2;

if ~exist([datDir,'figures/'])
    mkdir([datDir,'figures/']);
end

for ii = 1:2
    
    for jj = 1:5
        
        % Load real stats
        thisCirc = load([realDatDir,'dispHist',region{1},'_',imRegion{jj},imgSet{ii}]);
        circDat = thisCirc.(['dispHistCirc',suff{jj}]);
        thisV1 = load([realDatDir,'dispHist',region{2},'_',imRegion{jj},imgSet{ii}]);
        v1Dat = thisV1.(['dispHistV1',suff{jj}]);
        thisV2 = load([realDatDir,'dispHist',region{3},'_',imRegion{jj},imgSet{ii}]);
        v2Dat = thisV2.(['dispHistV2',suff{jj}]);
        thisMT = load([realDatDir,'dispHist',region{4},'_',imRegion{jj},imgSet{ii}]);
        MTDat = thisMT.(['dispHistMT',suff{jj}]);
        
        % Plot
        fig(ii,jj) = figure;
        fig(ii,jj).Position = [100+100*(ii-1) 100+100*(jj-1) 650 650];
        hold on;
        
        plot([0 0],[0 5],'--k','linewidth',2);
        
        fill([cntrDisp fliplr(cntrDisp) cntrDisp(1)],...
             [dispDat.(imgSet{ii}).(region{1}).(imRegionUCI{jj}) fliplr(dispDat.(imgSet{ii}).(region{1}).(imRegionLCI{jj}))...
              dispDat.(imgSet{ii}).(region{1}).(imRegionUCI{jj})(1)],'k','edgecolor','none','facealpha',0.5);
        fill([cntrDisp fliplr(cntrDisp) cntrDisp(1)],...
             [dispDat.(imgSet{ii}).(region{2}).(imRegionUCI{jj}) fliplr(dispDat.(imgSet{ii}).(region{2}).(imRegionLCI{jj}))...
              dispDat.(imgSet{ii}).(region{2}).(imRegionUCI{jj})(1)],'b','edgecolor','none','facealpha',0.5);
        fill([cntrDisp fliplr(cntrDisp) cntrDisp(1)],...
             [dispDat.(imgSet{ii}).(region{3}).(imRegionUCI{jj}) fliplr(dispDat.(imgSet{ii}).(region{3}).(imRegionLCI{jj}))...
              dispDat.(imgSet{ii}).(region{3}).(imRegionUCI{jj})(1)],'g','edgecolor','none','facealpha',0.5);
        fill([cntrDisp fliplr(cntrDisp) cntrDisp(1)],...
             [dispDat.(imgSet{ii}).(region{4}).(imRegionUCI{jj}) fliplr(dispDat.(imgSet{ii}).(region{4}).(imRegionLCI{jj}))...
              dispDat.(imgSet{ii}).(region{4}).(imRegionUCI{jj})(1)],'r','edgecolor','none','facealpha',0.5);

%         plot(cntrDisp,dispDat.(imgSet{ii}).(region{1}).(imRegionF{jj}),'color',[0.75 0.75 0.75],'linewidth',1);
%         plot(cntrDisp,dispDat.(imgSet{ii}).(region{2}).(imRegionF{jj}),'color',[0.75 0.75 1],'linewidth',1);
%         plot(cntrDisp,dispDat.(imgSet{ii}).(region{3}).(imRegionF{jj}),'color',[0.75 1 0.75],'linewidth',1);
%         plot(cntrDisp,dispDat.(imgSet{ii}).(region{4}).(imRegionF{jj}),'color',[1 0.75 0.75],'linewidth',1);
%         
        p(1) = plot(cntrDisp,circDat,'color',[0 0 0],'linewidth',2);
        p(2) = plot(cntrDisp,v1Dat,'color',[0 0 1],'linewidth',2);
        p(3) = plot(cntrDisp,v2Dat,'color',[0 1 0],'linewidth',2);
        p(4) = plot(cntrDisp,MTDat,'color',[1 0 0],'linewidth',2);
        
        set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1]);
        legend(p,{'Cent 10\circ','V1','V2','MT'});
        xlabel('Horizontal disparity (\circ)');
        ylabel('Probability density');
        title([imRegionF{jj},' - ',imgSet{ii}]);
        
        saveas(fig(ii,jj),[datDir,'figures/',imRegion{jj},'_',imgSet{ii},'.svg']);
        
    end
    
end

end