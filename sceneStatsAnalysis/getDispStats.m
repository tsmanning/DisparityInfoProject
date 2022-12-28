function [edges_disp,f1,varargout] = getDispStats(imSet,subset,res,resampIter)
% For image set what are the disparity statistics?

% Where am I?
splPath  = regexp(which('getDispStats'),filesep,'split');
rootDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
imDir    = [rootDir,'BORISimageSet',filesep];
KSDdir   = [rootDir,'savedKSDmatFiles_BORISdataset',filesep];
statsDir = [rootDir,'savedImageStats_BORISdataset',filesep];


%% Collect disparities within masks and histogram

% Define disparity boundaries of histogram (deg)
ub  = 2;
lb  = -2;

edges_disp = linspace(lb,ub,res);

% Labels for differential masking
switch subset
    case 'ecc'
        lab{1} = 'Peripheral';
        lab{2} = 'Central';
        
    case 'vertPos'
        lab{1} = 'Upper VF';
        lab{2} = 'Lower VF';
end

% Grab kernel-smoothed RF probability densities from datasets
load([KSDdir,'V1densityMat_BORIS.mat'])
load([KSDdir,'V2densityMat_BORIS.mat'])
load([KSDdir,'MTdensityMat_BORIS.mat'])
load([KSDdir,'circDensityMat_BORIS.mat'])
load([KSDdir,'V1V2Rect_BORIS.mat'])

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
        dispHistV1V2Rect = nan(numIms,res-1);

    otherwise
        dispHistV1_m1 = nan(numIms,res-1);
        dispHistV2_m1 = nan(numIms,res-1);
        dispHistMT_m1 = nan(numIms,res-1);
        dispHistCirc_m1 = nan(numIms,res-1);
        dispHistV1V2Rect_m1 = nan(numIms,res-1);
        
        dispHistV1_m2 = nan(numIms,res-1);
        dispHistV2_m2 = nan(numIms,res-1);
        dispHistMT_m2 = nan(numIms,res-1);
        dispHistCirc_m2 = nan(numIms,res-1);
        dispHistV1V2Rect_m2 = nan(numIms,res-1);
end

% Loop over all images in dataset
for ii = 1:numIms
        
    if mod(ii,50) == 0
    disp(['Running image ',num2str(ii),'/',num2str(numIms)]);
    end
    
    % For bootstrapping, select a random image from set
    if resampIter ~= 0
        imInd = randi(numIms);
    else
        imInd = ii;
    end


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
            xSampV1V2Rect = nan(numSamps,1);
            ySampV1V2Rect = nan(numSamps,1);

        otherwise
            xSampV1_m1    = nan(numSamps,1);
            ySampV1_m1    = nan(numSamps,1);
            xSampV2_m1    = nan(numSamps,1);
            ySampV2_m1    = nan(numSamps,1);
            xSampMT_m1    = nan(numSamps,1);
            ySampMT_m1    = nan(numSamps,1);
            xSampCirc_m1  = nan(numSamps,1);
            ySampCirc_m1  = nan(numSamps,1);
            xSampV1V2Rect_m1 = nan(numSamps,1);
            ySampV1V2Rect_m1 = nan(numSamps,1);

            xSampV1_m2    = nan(numSamps,1);
            ySampV1_m2    = nan(numSamps,1);
            xSampV2_m2    = nan(numSamps,1);
            ySampV2_m2    = nan(numSamps,1);
            xSampMT_m2    = nan(numSamps,1);
            ySampMT_m2    = nan(numSamps,1);
            xSampCirc_m2  = nan(numSamps,1);
            ySampCirc_m2  = nan(numSamps,1);
            xSampV1V2Rect_m2 = nan(numSamps,1);
            ySampV1V2Rect_m2 = nan(numSamps,1);
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
            V1V2RectMatMasked    = V1V2RectMat.*mask.*undefImMask;
            
        otherwise
            V1densityMatMasked1   = V1densityMat.*mask1.*undefImMask;
            V2densityMatMasked1   = V2densityMat.*mask1.*undefImMask;
            MTdensityMatMasked1   = MTdensityMat.*mask1.*undefImMask;
            circDensityMatMasked1 = circDensityMat.*mask1.*undefImMask;
            V1V2RectMatMasked1    = V1V2RectMat.*mask1.*undefImMask;
            
            V1densityMatMasked2   = V1densityMat.*mask2.*undefImMask;
            V2densityMatMasked2   = V2densityMat.*mask2.*undefImMask;
            MTdensityMatMasked2   = MTdensityMat.*mask2.*undefImMask;
            circDensityMatMasked2 = circDensityMat.*mask2.*undefImMask;
            V1V2RectMatMasked2    = V1V2RectMat.*mask2.*undefImMask;
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
            check(5) = sum(V1V2RectMatMasked(:));
            
            allCheck = sum([check(1) == 0; check(2) == 0; check(3) == 0; check(4) == 0;  check(5) == 0]);
            
        otherwise
            check(1) = sum(V1densityMatMasked1(:));
            check(2) = sum(V2densityMatMasked1(:));
            check(3) = sum(MTdensityMatMasked1(:));
            check(4) = sum(circDensityMatMasked1(:));
            check(5) = sum(V1V2RectMatMasked1(:));
            
            check(6) = sum(V1densityMatMasked2(:));
            check(7) = sum(V2densityMatMasked2(:));
            check(8) = sum(MTdensityMatMasked2(:));
            check(9) = sum(circDensityMatMasked2(:));
            check(10) = sum(V1V2RectMatMasked2(:));
            
            allCheck = sum([check(1) == 0; check(2) == 0; check(3) == 0; check(4) == 0;...
                            check(5) == 0; check(6) == 0; check(7) == 0; check(8) == 0;...
                            check(9) == 0; check(10) == 0]);
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
                [ySampV1V2Rect(jj),xSampV1V2Rect(jj)]     = pinky(1:imSize,1:imSize,V1V2RectMatMasked);

            otherwise
                [ySampV1_m1(jj),xSampV1_m1(jj)]     = pinky(1:imSize,1:imSize,V1densityMatMasked1);
                [ySampV2_m1(jj),xSampV2_m1(jj)]     = pinky(1:imSize,1:imSize,V2densityMatMasked1);
                [ySampMT_m1(jj),xSampMT_m1(jj)]     = pinky(1:imSize,1:imSize,MTdensityMatMasked1);
                [ySampCirc_m1(jj),xSampCirc_m1(jj)] = pinky(1:imSize,1:imSize,circDensityMatMasked1);
                [ySampV1V2Rect_m1(jj),xSampV1V2Rect_m1(jj)]     = pinky(1:imSize,1:imSize,V1V2RectMatMasked1);
                
                [ySampV1_m2(jj),xSampV1_m2(jj)]     = pinky(1:imSize,1:imSize,V1densityMatMasked2);
                [ySampV2_m2(jj),xSampV2_m2(jj)]     = pinky(1:imSize,1:imSize,V2densityMatMasked2);
                [ySampMT_m2(jj),xSampMT_m2(jj)]     = pinky(1:imSize,1:imSize,MTdensityMatMasked2);
                [ySampCirc_m2(jj),xSampCirc_m2(jj)] = pinky(1:imSize,1:imSize,circDensityMatMasked2);
                [ySampV1V2Rect_m2(jj),xSampV1V2Rect_m2(jj)]     = pinky(1:imSize,1:imSize,V1V2RectMatMasked2);
        end
        
    end
    
    % convert to linear indices
    switch subset
        case 'all'
            sampIndV1   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV1,ySampV1);
            sampIndV2   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV2,ySampV2);
            sampIndMT   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampMT,ySampMT);
            sampIndCirc = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampCirc,ySampCirc);
            sampIndV1V2Rect   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV1V2Rect,ySampV1V2Rect);
        otherwise
            sampIndV1_m1   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV1_m1,ySampV1_m1);
            sampIndV2_m1   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV2_m1,ySampV2_m1);
            sampIndMT_m1   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampMT_m1,ySampMT_m1);
            sampIndCirc_m1 = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampCirc_m1,ySampCirc_m1);
            sampIndV1V2Rect_m1   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV1V2Rect_m1,ySampV1V2Rect_m1);

            sampIndV1_m2   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV1_m2,ySampV1_m2);
            sampIndV2_m2   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV2_m2,ySampV2_m2);
            sampIndMT_m2   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampMT_m2,ySampMT_m2);
            sampIndCirc_m2 = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampCirc_m2,ySampCirc_m2);
            sampIndV1V2Rect_m2   = sub2ind([size(dispCrop,1) size(dispCrop,2)],xSampV1V2Rect_m2,ySampV1V2Rect_m2);
    end
    
    % Select disparities using these indices
    switch subset
        case 'all'
            dispResampV1   = dispCrop(sampIndV1);
            dispResampV2   = dispCrop(sampIndV2);
            dispResampMT   = dispCrop(sampIndMT);
            dispResampCirc = dispCrop(sampIndCirc);
            dispResampV1V2Rect   = dispCrop(sampIndV1V2Rect);
        otherwise
            dispResampV1_m1   = dispCrop(sampIndV1_m1);
            dispResampV2_m1   = dispCrop(sampIndV2_m1);
            dispResampMT_m1   = dispCrop(sampIndMT_m1);
            dispResampCirc_m1 = dispCrop(sampIndCirc_m1);
            dispResampV1V2Rect_m1   = dispCrop(sampIndV1V2Rect_m1);
            
            dispResampV1_m2   = dispCrop(sampIndV1_m2);
            dispResampV2_m2   = dispCrop(sampIndV2_m2);
            dispResampMT_m2   = dispCrop(sampIndMT_m2);
            dispResampCirc_m2 = dispCrop(sampIndCirc_m2);
            dispResampV1V2Rect_m2   = dispCrop(sampIndV1V2Rect_m2);
    end
    
    % Calculate final histograms
    switch subset
        case 'all'
            dispHistV1(ii,:)   = histcounts(dispResampV1,edges_disp);
            dispHistV2(ii,:)   = histcounts(dispResampV2,edges_disp);
            dispHistMT(ii,:)   = histcounts(dispResampMT,edges_disp);
            dispHistCirc(ii,:) = histcounts(dispResampCirc,edges_disp);
            dispHistV1V2Rect(ii,:) = histcounts(dispResampV1V2Rect,edges_disp);
        otherwise
            dispHistV1_m1(ii,:)   = histcounts(dispResampV1_m1,edges_disp);
            dispHistV2_m1(ii,:)   = histcounts(dispResampV2_m1,edges_disp);
            dispHistMT_m1(ii,:)   = histcounts(dispResampMT_m1,edges_disp);
            dispHistCirc_m1(ii,:) = histcounts(dispResampCirc_m1,edges_disp);
            dispHistV1V2Rect_m1(ii,:) = histcounts(dispResampV1V2Rect_m1,edges_disp);
            
            dispHistV1_m2(ii,:)   = histcounts(dispResampV1_m2,edges_disp);
            dispHistV2_m2(ii,:)   = histcounts(dispResampV2_m2,edges_disp);
            dispHistMT_m2(ii,:)   = histcounts(dispResampMT_m2,edges_disp);
            dispHistCirc_m2(ii,:) = histcounts(dispResampCirc_m2,edges_disp);
            dispHistV1V2Rect_m2(ii,:) = histcounts(dispResampV1V2Rect_m2,edges_disp);
    end
    
end


% Collapse across images
switch subset
    case 'all'
        dispHistV1   = sum(dispHistV1,1,'omitnan');
        dispHistV2   = sum(dispHistV2,1,'omitnan');
        dispHistMT   = sum(dispHistMT,1,'omitnan');
        dispHistCirc = sum(dispHistCirc,1,'omitnan');
        dispHistV1V2Rect   = sum(dispHistV1V2Rect,1,'omitnan');
    otherwise
        dispHistV1_m1   = sum(dispHistV1_m1,1,'omitnan');
        dispHistV2_m1   = sum(dispHistV2_m1,1,'omitnan');
        dispHistMT_m1   = sum(dispHistMT_m1,1,'omitnan');
        dispHistCirc_m1 = sum(dispHistCirc_m1,1,'omitnan');
        dispHistV1V2Rect_m1   = sum(dispHistV1V2Rect_m1,1,'omitnan');
        
        dispHistV1_m2   = sum(dispHistV1_m2,1,'omitnan');
        dispHistV2_m2   = sum(dispHistV2_m2,1,'omitnan');
        dispHistMT_m2   = sum(dispHistMT_m2,1,'omitnan');
        dispHistCirc_m2 = sum(dispHistCirc_m2,1,'omitnan');
        dispHistV1V2Rect_m2   = sum(dispHistV1V2Rect_m2,1,'omitnan');
end

% Normalize to get PDF
switch subset
    case 'all'
        dispHistV1   = dispHistV1/(sum(dispHistV1)*diff(edges_disp(1:2)));
        dispHistV2   = dispHistV2/(sum(dispHistV2)*diff(edges_disp(1:2)));
        dispHistMT   = dispHistMT/(sum(dispHistMT)*diff(edges_disp(1:2)));
        dispHistCirc = dispHistCirc/(sum(dispHistCirc)*diff(edges_disp(1:2)));
        dispHistV1V2Rect   = dispHistV1V2Rect/(sum(dispHistV1V2Rect)*diff(edges_disp(1:2)));
    otherwise
        dispHistV1_m1   = dispHistV1_m1/(sum(dispHistV1_m1)*diff(edges_disp(1:2)));
        dispHistV2_m1   = dispHistV2_m1/(sum(dispHistV2_m1)*diff(edges_disp(1:2)));
        dispHistMT_m1   = dispHistMT_m1/(sum(dispHistMT_m1)*diff(edges_disp(1:2)));
        dispHistCirc_m1 = dispHistCirc_m1/(sum(dispHistCirc_m1)*diff(edges_disp(1:2)));
        dispHistV1V2Rect_m1   = dispHistV1V2Rect_m1/(sum(dispHistV1V2Rect_m1)*diff(edges_disp(1:2)));
        
        dispHistV1_m2   = dispHistV1_m2/(sum(dispHistV1_m2)*diff(edges_disp(1:2)));
        dispHistV2_m2   = dispHistV2_m2/(sum(dispHistV2_m2)*diff(edges_disp(1:2)));
        dispHistMT_m2   = dispHistMT_m2/(sum(dispHistMT_m2)*diff(edges_disp(1:2)));
        dispHistCirc_m2 = dispHistCirc_m2/(sum(dispHistCirc_m2)*diff(edges_disp(1:2)));
        dispHistV1V2Rect_m2   = dispHistV1V2Rect_m2/(sum(dispHistV1V2Rect_m2)*diff(edges_disp(1:2)));
end

% Define output arguments
switch subset
    case 'all' 
        varargout{1} = dispHistV1;
        varargout{2} = dispHistV2;
        varargout{3} = dispHistMT;
        varargout{4} = dispHistCirc;
        varargout{5} = dispHistV1V2Rect;
        
    otherwise
        varargout{1} = dispHistV1_m1;
        varargout{2} = dispHistV2_m1;
        varargout{3} = dispHistMT_m1;
        varargout{4} = dispHistCirc_m1;
        varargout{5} = dispHistV1V2Rect_m1;
        
        varargout{6} = dispHistV1_m2;
        varargout{7} = dispHistV2_m2;
        varargout{8} = dispHistMT_m2;
        varargout{9} = dispHistCirc_m2;
        varargout{10} = dispHistV1V2Rect_m2;
end

%% Save
if resampIter ~= 0
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
        save([statsDir,'borisStats/dispHistV1V2Rect_',imSet,suffix],'dispHistV1V2Rect');
    otherwise
        save([statsDir,'borisStats/dispHistV1_',lab{1},'_',imSet,suffix],'dispHistV1_m1');
        save([statsDir,'borisStats/dispHistV2_',lab{1},'_',imSet,suffix],'dispHistV2_m1');
        save([statsDir,'borisStats/dispHistMT_',lab{1},'_',imSet,suffix],'dispHistMT_m1');
        save([statsDir,'borisStats/dispHistCirc_',lab{1},'_',imSet,suffix],'dispHistCirc_m1');
        save([statsDir,'borisStats/dispHistV1V2Rect_',lab{1},'_',imSet,suffix],'dispHistV1V2Rect_m1');
        
        save([statsDir,'borisStats/dispHistV1_',lab{2},'_',imSet,suffix],'dispHistV1_m2');
        save([statsDir,'borisStats/dispHistV2_',lab{2},'_',imSet,suffix],'dispHistV2_m2');
        save([statsDir,'borisStats/dispHistMT_',lab{2},'_',imSet,suffix],'dispHistMT_m2');
        save([statsDir,'borisStats/dispHistCirc_',lab{2},'_',imSet,suffix],'dispHistCirc_m2');
        save([statsDir,'borisStats/dispHistV1V2Rect_',lab{2},'_',imSet,suffix],'dispHistV1V2Rect_m2');
end

end


