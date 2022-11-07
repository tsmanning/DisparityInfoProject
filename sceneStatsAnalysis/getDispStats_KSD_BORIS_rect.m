function [edges_disp,f1,varargout] = getDispStats_KSD_BORIS_rect(imSet,subset,plotOn,resampIter)
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
load([statsDir,'V1V2Rect_BORIS.mat'])

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
        dispHistV1V2Rect = nan(numIms,res-1);
        
    otherwise
        dispHistV1V2Rect_m1 = nan(numIms,res-1);
        
        dispHistV1V2Rect_m2 = nan(numIms,res-1);
        
end


for ii = 1:numIms
        
    if mod(ii,50) == 0
    disp(['Running image ',num2str(ii),'/',num2str(numIms)]);
    end
    
    % For bootstrapping, select a random image from set
%     imInd = randi(numIms);
    imInd = ii;


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
    imSize   = size(V1V2RectMat,1);
    
    numSamps = 1000;
    
    switch subset
        case 'all'
            xSampV1V2Rect    = nan(numSamps,1);
            ySampV1V2Rect    = nan(numSamps,1);

        otherwise
            xSampV1V2Rect_m1    = nan(numSamps,1);
            ySampV1V2Rect_m1    = nan(numSamps,1);

            
            xSampV1V2Rect_m2    = nan(numSamps,1);
            ySampV1V2Rect_m2    = nan(numSamps,1);

    end
    
    % First modify KSD mats so we mask out pixels with undefined
    % disparities and/or unwanted ROIs in the VF
    undefImMask = ~isnan(dispCrop);
    
    switch subset
        case 'all'
            V1V2RectMatMasked   = V1V2RectMat.*mask.*undefImMask;
            
        otherwise
            V1V2RectMatMasked1   = V1V2RectMat.*mask1.*undefImMask;
            V1V2RectMatMasked2   = V1V2RectMat.*mask2.*undefImMask;

    end
    
    % Sometimes the KSD plots are all zero after masking if some of the
    % images are full of undefined regions. If one of these images is
    % encountered, just skip it
    switch subset
        case 'all'
            check(1) = sum(V1V2RectMatMasked(:));
            
            allCheck = sum([check(1) == 0]);
            
        otherwise
            check(1) = sum(V1V2RectMatMasked1(:));
            
            check(5) = sum(V1V2RectMatMasked2(:));
            
            allCheck = sum([check(1) == 0; check(5) == 0;]);
    end
    
    if allCheck
        continue
    end
    
    % Sample indices of disparity plots based on KSD
    for jj = 1:numSamps
        switch subset
            case 'all'
                [ySampV1V2Rect(jj),xSampV1V2Rect(jj)]     = pinky(1:imSize,1:imSize,V1V2RectMatMasked);

            otherwise
                [ySampV1V2Rect_m1(jj),xSampV1V2Rect_m1(jj)]     = pinky(1:imSize,1:imSize,V1V2RectMatMasked1);
                
                [ySampV1V2Rect_m2(jj),xSampV1V2Rect_m2(jj)]     = pinky(1:imSize,1:imSize,V1V2RectMatMasked2);

        end
        
    end
    
    % convert to linear indices
    switch subset
        case 'all'
            sampIndV1V2Rect   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV1V2Rect,ySampV1V2Rect);

        otherwise
            sampIndV1V2Rect_m1   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV1V2Rect_m1,ySampV1V2Rect_m1);
            
            sampIndV1V2Rect_m2   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV1V2Rect_m2,ySampV1V2Rect_m2);
            
    end
    
    % Select disparities using these indices
    switch subset
        case 'all'
            dispResampV1V2Rect   = dispCrop(sampIndV1V2Rect);

        otherwise
            dispResampV1V2Rect_m1   = dispCrop(sampIndV1V2Rect_m1);
            
            dispResampV1V2Rect_m2   = dispCrop(sampIndV1V2Rect_m2);

    end
    
    % Calculate final histograms
    switch subset
        case 'all'
            dispHistV1V2Rect(ii,:)   = histcounts(dispResampV1V2Rect,edges_disp);

        otherwise
            dispHistV1V2Rect_m1(ii,:)   = histcounts(dispResampV1V2Rect_m1,edges_disp);
            
            dispHistV1V2Rect_m2(ii,:)   = histcounts(dispResampV1V2Rect_m2,edges_disp);

    end
    
end


%% Plotting


% Collapse across images
switch subset
    case 'all'
        dispHistV1V2Rect   = sum(dispHistV1V2Rect,1,'omitnan');

    otherwise
        dispHistV1V2Rect_m1   = sum(dispHistV1V2Rect_m1,1,'omitnan');
        
        dispHistV1V2Rect_m2   = sum(dispHistV1V2Rect_m2,1,'omitnan');

end

% Normalize to get PDF
switch subset
    case 'all'
        dispHistV1V2Rect   = dispHistV1V2Rect/(sum(dispHistV1V2Rect)*diff(edges_disp(1:2)));

    otherwise
        dispHistV1V2Rect_m1   = dispHistV1V2Rect_m1/(sum(dispHistV1V2Rect_m1)*diff(edges_disp(1:2)));
        
        dispHistV1V2Rect_m2   = dispHistV1V2Rect_m2/(sum(dispHistV1V2Rect_m2)*diff(edges_disp(1:2)));

end

if plotOn
% % Colors to use in plotting
% colorMat = colororder;  
% 
% switch subset
%     case 'ecc'
%         colors = [3 4];
%         labels = {['Central VF (<',num2str(round(medEcc)),'\circ)'],...
%                   ['Peripheral VF (>',num2str(round(medEcc)),'\circ)']};
%         
%     case 'vertPos'
%         colors = [5 6];
%         labels = {'Upper VF','Lower VF'};
%         
%     case 'figure'
%         colors = [1 2];
%         labels = {'Figure','Ground'};    
%         
% end
%     
% f1 = figure;
% f1.Position = [300 300 720 770];
% hold on;
% 
% switch subset
%     case 'all' 
%         % V1-based KSD sampling
%         f2 = figure;
%         f2.Position = [100 100 720 770];
%         plot([0 0],[0 pmax],'--k','linewidth',2);
%         p2 = plot(cntr_disp,dispHistV1,'color',[0 0 0],'linewidth',4);
%         title('Disparity probability');
%         set(gca,'fontsize',20,'xlim',[lb ub],'ylim',[0 pmax],'plotboxaspectratio',[1 1 1]);
%         xlabel('Horizontal disparity (\circ)');
%         ylabel('Probability density');
%         legend([p2],'V1','location','northeast');
%         
%         % V2-based KSD sampling
%         f3 = figure;
%         f3.Position = [700 100 720 770];
%         plot([0 0],[0 pmax],'--k','linewidth',2);
%         p3 = plot(cntr_disp,dispHistV2,'color',[0 0 0],'linewidth',4);
%         title('Disparity probability');
%         set(gca,'fontsize',20,'xlim',[lb ub],'ylim',[0 pmax],'plotboxaspectratio',[1 1 1]);
%         xlabel('Horizontal disparity (\circ)');
%         ylabel('Probability density');
%         legend([p3],'V2','location','northeast');
%         
%         % MT-based KSD sampling
%         f4 = figure;
%         f4.Position = [700 700 720 770];
%         plot([0 0],[0 pmax],'--k','linewidth',2);
%         p4 = plot(cntr_disp,dispHistMT,'color',[0 0 0],'linewidth',4);
%         title('Disparity probability');
%         set(gca,'fontsize',20,'xlim',[lb ub],'ylim',[0 pmax],'plotboxaspectratio',[1 1 1]);
%         xlabel('Horizontal disparity (\circ)');
%         ylabel('Probability density');
%         legend([p4],'MT','location','northeast');
%         
%         % Circ-based KSD sampling
%         f4 = figure;
%         f4.Position = [700 700 720 770];
%         plot([0 0],[0 pmax],'--k','linewidth',2);
%         p4 = plot(cntr_disp,dispHistCirc,'color',[0 0 0],'linewidth',4);
%         title('Disparity probability');
%         set(gca,'fontsize',20,'xlim',[lb ub],'ylim',[0 pmax],'plotboxaspectratio',[1 1 1]);
%         xlabel('Horizontal disparity (\circ)');
%         ylabel('Probability density');
%         legend([p4],'Cent 10\circ','location','northeast');
%         
%     otherwise
%         % V1-based KSD sampling
%         f2 = figure;
%         f2.Position = [100 100 720 770];
%         hold on
%         plot([0 0],[0 pmax],'--k','linewidth',2);
%         p1 = plot(cntr_disp,dispHistV1_m1,'color',[0 0 0],'linewidth',4);
%         p1b = plot(cntr_disp,dispHistV1_m2,'color',[0 0 0],'linewidth',4,'linestyle','--');
%         title('Disparity probability');
%         set(gca,'fontsize',20,'xlim',[lb ub],'ylim',[0 pmax],'plotboxaspectratio',[1 1 1]);
%         xlabel('Horizontal disparity (\circ)');
%         ylabel('Probability density');
%         legend([p1,p1b],{['V1 - ',lab{1}],['V1 - ',lab{2}]},'location','northeast');
%         
%         % V2-based KSD sampling
%         f3 = figure;
%         f3.Position = [700 100 720 770];
%         hold on
%         plot([0 0],[0 pmax],'--k','linewidth',2);
%         p2 = plot(cntr_disp,dispHistV2_m1,'color',[0 0 0],'linewidth',4);
%         p2b = plot(cntr_disp,dispHistV2_m2,'color',[0 0 0],'linewidth',4,'linestyle','--');
%         title('Disparity probability');
%         set(gca,'fontsize',20,'xlim',[lb ub],'ylim',[0 pmax],'plotboxaspectratio',[1 1 1]);
%         xlabel('Horizontal disparity (\circ)');
%         ylabel('Probability density');
%         legend([p2,p2b],{['V2 - ',lab{1}],['V2 - ',lab{2}]},'location','northeast');
%         
%         % MT-based KSD sampling
%         f4 = figure;
%         f4.Position = [700 700 720 770];
%         hold on
%         plot([0 0],[0 pmax],'--k','linewidth',2);
%         p3 = plot(cntr_disp,dispHistMT_m1,'color',[0 0 0],'linewidth',4);
%         p3b = plot(cntr_disp,dispHistMT_m2,'color',[0 0 0],'linewidth',4,'linestyle','--');
%         title('Disparity probability');
%         set(gca,'fontsize',20,'xlim',[lb ub],'ylim',[0 pmax],'plotboxaspectratio',[1 1 1]);
%         xlabel('Horizontal disparity (\circ)');
%         ylabel('Probability density');
%         legend([p3,p3b],{['MT - ',lab{1}],['MT - ',lab{2}]},'location','northeast');
%         
%         % Circ-based KSD sampling
%         f4 = figure;
%         f4.Position = [700 700 720 770];
%         hold on
%         plot([0 0],[0 pmax],'--k','linewidth',2);
%         p4 = plot(cntr_disp,dispHistCirc_m1,'color',[0 0 0],'linewidth',4);
%         p4b = plot(cntr_disp,dispHistCirc_m2,'color',[0 0 0],'linewidth',4,'linestyle','--');
%         title('Disparity probability');
%         set(gca,'fontsize',20,'xlim',[lb ub],'ylim',[0 pmax],'plotboxaspectratio',[1 1 1]);
%         xlabel('Horizontal disparity (\circ)');
%         ylabel('Probability density');
%         legend([p4,p4b],{['Cent 10\circ - ',lab{1}],['Cent 10\circ - ',lab{2}]},'location','northeast');
% end

end

% Define output arguments
switch subset
    case 'all' 
        varargout{1} = dispHistV1V2Rect;
        
    otherwise
        varargout{1} = dispHistV1V2Rect_m1;
        
        varargout{2} = dispHistV1V2Rect_m2;
        
end

%% Save
if resampIter > 0
    suffix = num2str(resampIter);
else
    suffix = '';
end

if ~exist([statsDir,'borisStats/'])
    mkdir([statsDir,'borisStats/']);
end
    
switch subset
    
    case 'all'
        save([statsDir,'borisStats/dispHistV1V2Rect_',imSet,suffix],'dispHistV1V2Rect');

    otherwise
        save([statsDir,'borisStats/dispHistV1V2Rect_',lab{1},'_',imSet,suffix],'dispHistV1V2Rect_m1');
        
        save([statsDir,'borisStats/dispHistV1V2Rect_',lab{2},'_',imSet,suffix],'dispHistV1V2Rect_m2');
end

end


