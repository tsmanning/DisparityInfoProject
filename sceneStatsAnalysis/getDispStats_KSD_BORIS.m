function [edges_disp,f1,varargout] = getDispStats_KSD_BORIS(imSet,subset,plotOn,resampIter)
% For image set what are the disparity statistics?

% Where am I?
splPath  = regexp(which('getDispStats'),filesep,'split');
rootDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
imDir    = [rootDir,'NaturalImageDB/BORISimageSet/'];
statsDir = [rootDir,'5-sceneStatsAnalysis/'];

%% Collect disparities within masks and histogram

res        = 52;
% res        = 25;

ub = 2;
lb = -2;

edges_disp = linspace(lb,ub,res);
cntr_disp  = edges_disp(1:end-1) + diff(edges_disp(1:2))/2;

% Labels for differential masking
switch subset
    case 'ecc'
        lab{1} = 'Peripheral';
        lab{2} = 'Central';
        
    case 'vertPos'
        lab{1} = 'Upper VF';
        lab{2} = 'Lower VF';
end

pmax = 5;

% Grab kernel-smoothed RF probability densities from datasets
load([statsDir,'V1densityMat_BORIS.mat'])
load([statsDir,'V2densityMat_BORIS.mat'])
load([statsDir,'MTdensityMat_BORIS.mat'])
load([statsDir,'circDensityMat_BORIS.mat'])

% Grab disparity images from BORIS dataset
switch imSet
    case 'sando'
        load([imDir,'Making_sandwich_disparity.mat']);
        sando = horizontal_disparity;
        numIms = size(sando,3);
        
    case 'walking'
        load([imDir,'Walking_outside_disparity.mat']);
        walk  = horizontal_disparity;
        numIms  = size(walk,3);
end

switch subset
    case 'all'
        dispHistV1 = nan(numIms,res-1);
        dispHistV2 = nan(numIms,res-1);
        dispHistMT = nan(numIms,res-1);
        dispHistCirc = nan(numIms,res-1);
        
    otherwise
        dispHistV1_m1 = nan(numIms,res-1);
        dispHistV2_m1 = nan(numIms,res-1);
        dispHistMT_m1 = nan(numIms,res-1);
        dispHistCirc_m1 = nan(numIms,res-1);
        
        dispHistV1_m2 = nan(numIms,res-1);
        dispHistV2_m2 = nan(numIms,res-1);
        dispHistMT_m2 = nan(numIms,res-1);
        dispHistCirc_m2 = nan(numIms,res-1);
        
end


for ii = 1:numIms
        
    if mod(ii,50) == 0
    disp(['Running image ',num2str(ii),'/',num2str(numIms)]);
    end
    
    % For bootstrapping, select a random image from set
    imInd = randi(numIms);
    
    %% Get image size
    % in BORIS data, the matrix is 10 deg wide and 10 tall; sampled in 207x207
    supp1D = -103:103;
    
    [xgrid,ygrid] = meshgrid(supp1D,fliplr(supp1D));
    ecc           = sqrt(xgrid.^2 + ygrid.^2);
    ecc           = ecc*(10/103);
    
    
    %% Make a mask of the visual field subregion of interest
    switch subset
        case 'all'
            
            mask = ecc<10;
            
        case 'ecc'
            
            medEcc      = 5;
            lowerCutoff = 1;
            upperCutoff = 10;
            
            % Peripheral
            mask1 = (ecc>medEcc) & (ecc<upperCutoff);
            
            % Central
            mask2 = (ecc<medEcc) & (ecc>lowerCutoff);
            
        case 'vertPos'
            
            % Upper
            mask1 = ygrid<0;
            
            % Lower
            mask2 = ygrid>0;
            
    end
    
    % BORIS image set has row 1 = -10.3deg, row 2 = -10.3deg, etc so let's
    % flip it to match RF probability densities
    dispCrop = flipud(squeeze(horizontal_disparity(:,:,imInd)));

    
    %% Resample disparities based on RF location probability densities    
    imSize   = size(MTdensityMat,1);
    
    numSamps = 1000;
    
    switch subset
        case 'all'
            xSampV1    = nan(numSamps,1);
            ySampV1    = nan(numSamps,1);
            xSampV2    = nan(numSamps,1);
            ySampV2    = nan(numSamps,1);
            xSampMT    = nan(numSamps,1);
            ySampMT    = nan(numSamps,1);
            xSampCirc  = nan(numSamps,1);
            ySampCirc  = nan(numSamps,1);
        otherwise
            xSampV1_m1    = nan(numSamps,1);
            ySampV1_m1    = nan(numSamps,1);
            xSampV2_m1    = nan(numSamps,1);
            ySampV2_m1    = nan(numSamps,1);
            xSampMT_m1    = nan(numSamps,1);
            ySampMT_m1    = nan(numSamps,1);
            xSampCirc_m1  = nan(numSamps,1);
            ySampCirc_m1  = nan(numSamps,1);
            
            xSampV1_m2    = nan(numSamps,1);
            ySampV1_m2    = nan(numSamps,1);
            xSampV2_m2    = nan(numSamps,1);
            ySampV2_m2    = nan(numSamps,1);
            xSampMT_m2    = nan(numSamps,1);
            ySampMT_m2    = nan(numSamps,1);
            xSampCirc_m2  = nan(numSamps,1);
            ySampCirc_m2  = nan(numSamps,1);
    end
    
    % First modify KSD mats so we mask out pixels with undefined
    % disparities and/or unwanted ROIs in the VF
    undefImMask = ~isnan(dispCrop);
    
    switch subset
        case 'all'
            V1densityMatMasked   = V1densityMat.*mask.*undefImMask;
            V2densityMatMasked   = V2densityMat.*mask.*undefImMask;
            MTdensityMatMasked   = MTdensityMat.*mask.*undefImMask;
            circDensityMatMasked = circDensityMat.*mask.*undefImMask;
            
        otherwise
            V1densityMatMasked1   = V1densityMat.*mask1.*undefImMask;
            V2densityMatMasked1   = V2densityMat.*mask1.*undefImMask;
            MTdensityMatMasked1   = MTdensityMat.*mask1.*undefImMask;
            circDensityMatMasked1 = circDensityMat.*mask1.*undefImMask;
            
            V1densityMatMasked2   = V1densityMat.*mask2.*undefImMask;
            V2densityMatMasked2   = V2densityMat.*mask2.*undefImMask;
            MTdensityMatMasked2   = MTdensityMat.*mask2.*undefImMask;
            circDensityMatMasked2 = circDensityMat.*mask2.*undefImMask;
            
    end
    
    % Sometimes the KSD plots are all zero after masking if some of the
    % images are full of undefined regions. If one of these images is
    % encountered, just skip it
    switch subset
        case 'all'
            check(1) = sum(V1densityMatMasked(:));
            check(2) = sum(V2densityMatMasked(:));
            check(3) = sum(MTdensityMatMasked(:));
            check(4) = sum(circDensityMatMasked(:));
            
            allCheck = sum([check(1) == 0; check(2) == 0; check(3) == 0; check(4) == 0]);
            
        otherwise
            check(1) = sum(V1densityMatMasked1(:));
            check(2) = sum(V2densityMatMasked1(:));
            check(3) = sum(MTdensityMatMasked1(:));
            check(4) = sum(circDensityMatMasked1(:));
            
            check(5) = sum(V1densityMatMasked2(:));
            check(6) = sum(V2densityMatMasked2(:));
            check(7) = sum(MTdensityMatMasked2(:));
            check(8) = sum(circDensityMatMasked2(:));
            
            allCheck = sum([check(1) == 0; check(2) == 0; check(3) == 0; check(4) == 0;...
                            check(5) == 0; check(6) == 0; check(7) == 0; check(8) == 0]);
    end
    
    if allCheck
        continue
    end
    
    % Sample indices of disparity plots based on KSD
    for jj = 1:numSamps
        switch subset
            case 'all'
                [ySampV1(jj),xSampV1(jj)]     = pinky(1:imSize,1:imSize,V1densityMatMasked);
                [ySampV2(jj),xSampV2(jj)]     = pinky(1:imSize,1:imSize,V2densityMatMasked);
                [ySampMT(jj),xSampMT(jj)]     = pinky(1:imSize,1:imSize,MTdensityMatMasked);
                [ySampCirc(jj),xSampCirc(jj)] = pinky(1:imSize,1:imSize,circDensityMatMasked);
            otherwise
                [ySampV1_m1(jj),xSampV1_m1(jj)]     = pinky(1:imSize,1:imSize,V1densityMatMasked1);
                [ySampV2_m1(jj),xSampV2_m1(jj)]     = pinky(1:imSize,1:imSize,V2densityMatMasked1);
                [ySampMT_m1(jj),xSampMT_m1(jj)]     = pinky(1:imSize,1:imSize,MTdensityMatMasked1);
                [ySampCirc_m1(jj),xSampCirc_m1(jj)] = pinky(1:imSize,1:imSize,circDensityMatMasked1);
                
                [ySampV1_m2(jj),xSampV1_m2(jj)]     = pinky(1:imSize,1:imSize,V1densityMatMasked2);
                [ySampV2_m2(jj),xSampV2_m2(jj)]     = pinky(1:imSize,1:imSize,V2densityMatMasked2);
                [ySampMT_m2(jj),xSampMT_m2(jj)]     = pinky(1:imSize,1:imSize,MTdensityMatMasked2);
                [ySampCirc_m2(jj),xSampCirc_m2(jj)] = pinky(1:imSize,1:imSize,circDensityMatMasked2);
        end
        
    end
    
    % convert to linear indices
    switch subset
        case 'all'
            sampIndV1   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV1,ySampV1);
            sampIndV2   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV2,ySampV2);
            sampIndMT   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampMT,ySampMT);
            sampIndCirc = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampCirc,ySampCirc);
        otherwise
            sampIndV1_m1   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV1_m1,ySampV1_m1);
            sampIndV2_m1   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV2_m1,ySampV2_m1);
            sampIndMT_m1   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampMT_m1,ySampMT_m1);
            sampIndCirc_m1 = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampCirc_m1,ySampCirc_m1);
            
            sampIndV1_m2   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV1_m2,ySampV1_m2);
            sampIndV2_m2   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV2_m2,ySampV2_m2);
            sampIndMT_m2   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampMT_m2,ySampMT_m2);
            sampIndCirc_m2 = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampCirc_m2,ySampCirc_m2);
            
    end
    
    % Select disparities using these indices
    switch subset
        case 'all'
            dispResampV1   = dispCrop(sampIndV1);
            dispResampV2   = dispCrop(sampIndV2);
            dispResampMT   = dispCrop(sampIndMT);
            dispResampCirc = dispCrop(sampIndCirc);
        otherwise
            dispResampV1_m1   = dispCrop(sampIndV1_m1);
            dispResampV2_m1   = dispCrop(sampIndV2_m1);
            dispResampMT_m1   = dispCrop(sampIndMT_m1);
            dispResampCirc_m1 = dispCrop(sampIndCirc_m1);
            
            dispResampV1_m2   = dispCrop(sampIndV1_m2);
            dispResampV2_m2   = dispCrop(sampIndV2_m2);
            dispResampMT_m2   = dispCrop(sampIndMT_m2);
            dispResampCirc_m2 = dispCrop(sampIndCirc_m2);
    end
    
    % Calculate final histograms
    switch subset
        case 'all'
            dispHistV1(ii,:)   = histcounts(dispResampV1,edges_disp);
            dispHistV2(ii,:)   = histcounts(dispResampV2,edges_disp);
            dispHistMT(ii,:)   = histcounts(dispResampMT,edges_disp);
            dispHistCirc(ii,:) = histcounts(dispResampCirc,edges_disp);
        otherwise
            dispHistV1_m1(ii,:)   = histcounts(dispResampV1_m1,edges_disp);
            dispHistV2_m1(ii,:)   = histcounts(dispResampV2_m1,edges_disp);
            dispHistMT_m1(ii,:)   = histcounts(dispResampMT_m1,edges_disp);
            dispHistCirc_m1(ii,:) = histcounts(dispResampCirc_m1,edges_disp);
            
            dispHistV1_m2(ii,:)   = histcounts(dispResampV1_m2,edges_disp);
            dispHistV2_m2(ii,:)   = histcounts(dispResampV2_m2,edges_disp);
            dispHistMT_m2(ii,:)   = histcounts(dispResampMT_m2,edges_disp);
            dispHistCirc_m2(ii,:) = histcounts(dispResampCirc_m2,edges_disp);
    end
    
end


%% Plotting


% Collapse across images
switch subset
    case 'all'
        dispHistV1   = sum(dispHistV1,1,'omitnan');
        dispHistV2   = sum(dispHistV2,1,'omitnan');
        dispHistMT   = sum(dispHistMT,1,'omitnan');
        dispHistCirc = sum(dispHistCirc,1,'omitnan');
    otherwise
        dispHistV1_m1   = sum(dispHistV1_m1,1,'omitnan');
        dispHistV2_m1   = sum(dispHistV2_m1,1,'omitnan');
        dispHistMT_m1   = sum(dispHistMT_m1,1,'omitnan');
        dispHistCirc_m1 = sum(dispHistCirc_m1,1,'omitnan');
        
        dispHistV1_m2   = sum(dispHistV1_m2,1,'omitnan');
        dispHistV2_m2   = sum(dispHistV2_m2,1,'omitnan');
        dispHistMT_m2   = sum(dispHistMT_m2,1,'omitnan');
        dispHistCirc_m2 = sum(dispHistCirc_m2,1,'omitnan');
end

% Normalize to get PDF
switch subset
    case 'all'
        dispHistV1   = dispHistV1/(sum(dispHistV1)*diff(edges_disp(1:2)));
        dispHistV2   = dispHistV2/(sum(dispHistV2)*diff(edges_disp(1:2)));
        dispHistMT   = dispHistMT/(sum(dispHistMT)*diff(edges_disp(1:2)));
        dispHistCirc = dispHistCirc/(sum(dispHistCirc)*diff(edges_disp(1:2)));
    otherwise
        dispHistV1_m1   = dispHistV1_m1/(sum(dispHistV1_m1)*diff(edges_disp(1:2)));
        dispHistV2_m1   = dispHistV2_m1/(sum(dispHistV2_m1)*diff(edges_disp(1:2)));
        dispHistMT_m1   = dispHistMT_m1/(sum(dispHistMT_m1)*diff(edges_disp(1:2)));
        dispHistCirc_m1 = dispHistCirc_m1/(sum(dispHistCirc_m1)*diff(edges_disp(1:2)));
        
        dispHistV1_m2   = dispHistV1_m2/(sum(dispHistV1_m2)*diff(edges_disp(1:2)));
        dispHistV2_m2   = dispHistV2_m2/(sum(dispHistV2_m2)*diff(edges_disp(1:2)));
        dispHistMT_m2   = dispHistMT_m2/(sum(dispHistMT_m2)*diff(edges_disp(1:2)));
        dispHistCirc_m2 = dispHistCirc_m2/(sum(dispHistCirc_m2)*diff(edges_disp(1:2)));
end

if plotOn
% Colors to use in plotting
colorMat = colororder;  

switch subset
    case 'ecc'
        colors = [3 4];
        labels = {['Central VF (<',num2str(round(medEcc)),'\circ)'],...
                  ['Peripheral VF (>',num2str(round(medEcc)),'\circ)']};
        
    case 'vertPos'
        colors = [5 6];
        labels = {'Upper VF','Lower VF'};
        
    case 'figure'
        colors = [1 2];
        labels = {'Figure','Ground'};    
        
end
    
f1 = figure;
f1.Position = [300 300 720 770];
hold on;

switch subset
    case 'all' 
        % V1-based KSD sampling
        f2 = figure;
        f2.Position = [100 100 720 770];
        plot([0 0],[0 pmax],'--k','linewidth',2);
        p2 = plot(cntr_disp,dispHistV1,'color',[0 0 0],'linewidth',4);
        title('Disparity probability');
        set(gca,'fontsize',20,'xlim',[lb ub],'ylim',[0 pmax],'plotboxaspectratio',[1 1 1]);
        xlabel('Horizontal disparity (\circ)');
        ylabel('Probability density');
        legend([p2],'V1','location','northeast');
        
        % V2-based KSD sampling
        f3 = figure;
        f3.Position = [700 100 720 770];
        plot([0 0],[0 pmax],'--k','linewidth',2);
        p3 = plot(cntr_disp,dispHistV2,'color',[0 0 0],'linewidth',4);
        title('Disparity probability');
        set(gca,'fontsize',20,'xlim',[lb ub],'ylim',[0 pmax],'plotboxaspectratio',[1 1 1]);
        xlabel('Horizontal disparity (\circ)');
        ylabel('Probability density');
        legend([p3],'V2','location','northeast');
        
        % MT-based KSD sampling
        f4 = figure;
        f4.Position = [700 700 720 770];
        plot([0 0],[0 pmax],'--k','linewidth',2);
        p4 = plot(cntr_disp,dispHistMT,'color',[0 0 0],'linewidth',4);
        title('Disparity probability');
        set(gca,'fontsize',20,'xlim',[lb ub],'ylim',[0 pmax],'plotboxaspectratio',[1 1 1]);
        xlabel('Horizontal disparity (\circ)');
        ylabel('Probability density');
        legend([p4],'MT','location','northeast');
        
        % Circ-based KSD sampling
        f4 = figure;
        f4.Position = [700 700 720 770];
        plot([0 0],[0 pmax],'--k','linewidth',2);
        p4 = plot(cntr_disp,dispHistCirc,'color',[0 0 0],'linewidth',4);
        title('Disparity probability');
        set(gca,'fontsize',20,'xlim',[lb ub],'ylim',[0 pmax],'plotboxaspectratio',[1 1 1]);
        xlabel('Horizontal disparity (\circ)');
        ylabel('Probability density');
        legend([p4],'Cent 10\circ','location','northeast');
        
    otherwise
        % V1-based KSD sampling
        f2 = figure;
        f2.Position = [100 100 720 770];
        hold on
        plot([0 0],[0 pmax],'--k','linewidth',2);
        p1 = plot(cntr_disp,dispHistV1_m1,'color',[0 0 0],'linewidth',4);
        p1b = plot(cntr_disp,dispHistV1_m2,'color',[0 0 0],'linewidth',4,'linestyle','--');
        title('Disparity probability');
        set(gca,'fontsize',20,'xlim',[lb ub],'ylim',[0 pmax],'plotboxaspectratio',[1 1 1]);
        xlabel('Horizontal disparity (\circ)');
        ylabel('Probability density');
        legend([p1,p1b],{['V1 - ',lab{1}],['V1 - ',lab{2}]},'location','northeast');
        
        % V2-based KSD sampling
        f3 = figure;
        f3.Position = [700 100 720 770];
        hold on
        plot([0 0],[0 pmax],'--k','linewidth',2);
        p2 = plot(cntr_disp,dispHistV2_m1,'color',[0 0 0],'linewidth',4);
        p2b = plot(cntr_disp,dispHistV2_m2,'color',[0 0 0],'linewidth',4,'linestyle','--');
        title('Disparity probability');
        set(gca,'fontsize',20,'xlim',[lb ub],'ylim',[0 pmax],'plotboxaspectratio',[1 1 1]);
        xlabel('Horizontal disparity (\circ)');
        ylabel('Probability density');
        legend([p2,p2b],{['V2 - ',lab{1}],['V2 - ',lab{2}]},'location','northeast');
        
        % MT-based KSD sampling
        f4 = figure;
        f4.Position = [700 700 720 770];
        hold on
        plot([0 0],[0 pmax],'--k','linewidth',2);
        p3 = plot(cntr_disp,dispHistMT_m1,'color',[0 0 0],'linewidth',4);
        p3b = plot(cntr_disp,dispHistMT_m2,'color',[0 0 0],'linewidth',4,'linestyle','--');
        title('Disparity probability');
        set(gca,'fontsize',20,'xlim',[lb ub],'ylim',[0 pmax],'plotboxaspectratio',[1 1 1]);
        xlabel('Horizontal disparity (\circ)');
        ylabel('Probability density');
        legend([p3,p3b],{['MT - ',lab{1}],['MT - ',lab{2}]},'location','northeast');
        
        % Circ-based KSD sampling
        f4 = figure;
        f4.Position = [700 700 720 770];
        hold on
        plot([0 0],[0 pmax],'--k','linewidth',2);
        p4 = plot(cntr_disp,dispHistCirc_m1,'color',[0 0 0],'linewidth',4);
        p4b = plot(cntr_disp,dispHistCirc_m2,'color',[0 0 0],'linewidth',4,'linestyle','--');
        title('Disparity probability');
        set(gca,'fontsize',20,'xlim',[lb ub],'ylim',[0 pmax],'plotboxaspectratio',[1 1 1]);
        xlabel('Horizontal disparity (\circ)');
        ylabel('Probability density');
        legend([p4,p4b],{['Cent 10\circ - ',lab{1}],['Cent 10\circ - ',lab{2}]},'location','northeast');
end
end

% Define output arguments
switch subset
    case 'all' 
        varargout{1} = dispHistV1;
        varargout{2} = dispHistV2;
        varargout{3} = dispHistMT;
        varargout{4} = dispHistCirc;
        
    otherwise
        varargout{1} = dispHistV1_m1;
        varargout{2} = dispHistV2_m1;
        varargout{3} = dispHistMT_m1;
        varargout{4} = dispHistCirc_m1;
        
        varargout{5} = dispHistV1_m2;
        varargout{6} = dispHistV2_m2;
        varargout{7} = dispHistMT_m2;
        varargout{8} = dispHistCirc_m2;
        
end

%% Save
if nargin > 3
    suffix = num2str(resampIter);
else
    suffix = '';
end

if ~exist([statsDir,'borisStats/'])
    mkdir([statsDir,'borisStats/']);
end
    
switch subset
    
    case 'all'
        save([statsDir,'borisStats/dispHistV1_',imSet,suffix],'dispHistV1');
        save([statsDir,'borisStats/dispHistV2_',imSet,suffix],'dispHistV2');
        save([statsDir,'borisStats/dispHistMT_',imSet,suffix],'dispHistMT');
        save([statsDir,'borisStats/dispHistCirc_',imSet,suffix],'dispHistCirc');
    otherwise
        save([statsDir,'borisStats/dispHistV1_',lab{1},'_',imSet,suffix],'dispHistV1_m1');
        save([statsDir,'borisStats/dispHistV2_',lab{1},'_',imSet,suffix],'dispHistV2_m1');
        save([statsDir,'borisStats/dispHistMT_',lab{1},'_',imSet,suffix],'dispHistMT_m1');
        save([statsDir,'borisStats/dispHistCirc_',lab{1},'_',imSet,suffix],'dispHistCirc_m1');
        
        save([statsDir,'borisStats/dispHistV1_',lab{2},'_',imSet,suffix],'dispHistV1_m2');
        save([statsDir,'borisStats/dispHistV2_',lab{2},'_',imSet,suffix],'dispHistV2_m2');
        save([statsDir,'borisStats/dispHistMT_',lab{2},'_',imSet,suffix],'dispHistMT_m2');
        save([statsDir,'borisStats/dispHistCirc_',lab{2},'_',imSet,suffix],'dispHistCirc_m2');
end

end


