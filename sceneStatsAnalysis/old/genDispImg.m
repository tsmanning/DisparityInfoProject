% Convert range maps from image dataset to cyclopean image disparity maps

clear all
close all

% Save images?
saveOn = 1;
% saveOn = 0;

% How many image pairs?
% numSamps = 200;
numSamps = 1;

% Set fixation type (random or on boundary)
fpType = 'rand';
% fpType = 'border';
% fpType = 'vergence';
% fpType = 'figure';

% Plot example pairs and disparity map?
% plotEx = 1;
plotEx = 0;

% Plot RF array?
plotArr = 0;


%% Setup dirs and plotting constants

% Where am I?
splPath = regexp(which('genDispImg'),filesep,'split');
rootDir = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
neuronRootDir = [filesep,fullfile(splPath{1:numel(splPath)-2}),filesep];

% Add helper functions
addpath([rootDir,'helper_functions',filesep]);

% Set dir for depth maps 
% scene_path = [filesep,fullfile(splPath{1:numel(splPath)-2}),filesep,...
%               'NaturalImageDB/McCannDataset/CPSRange_full/RangeDatabase1080p',filesep,];
scene_path = [filesep,fullfile(splPath{1:numel(splPath)-2}),filesep,...
              'NaturalImageDB/McCannDataset/CPSRange/RangeDatabase1080p',filesep,];
          
% Redblue colormap
colors  = gray(256);
redblue = [ [ones(length(colors),1) colors(:,1) colors(:,1)]; 
            flipud([colors(:,1) colors(:,1) ones(length(colors),1)]) ];
          
          
%% Setup viewer constants & scene constraints
dims                    = [1080, 1920];
degreePerPixel          = 0.02;
hFov                    = dims(2)*degreePerPixel;
vFov                    = dims(1)*degreePerPixel;

% define eye positions. IPD = 6cm. eye1 = left eye
global eye_1; 
global eye_2;

eye_1  = [-0.03, 0, 0];
eye_2  = [0.03, 0, 0];

% coordinates of each image pixel in pixel units
vax = 1:dims(1);
hax = 1:dims(2);
[col_vals,row_vals] = meshgrid(hax,vax);

% empty matrix the size of our depth maps
emp = nan(dims(1),dims(2));

% Get IDs of valid scenes & fix numbering
load([rootDir,'metadata/valid_scenes.mat']);
load([rootDir,'metadata/dmap_indices.mat']);

%%%%%%%%%%% FIX THIS WHEN YOU HAVE FULL IMAGE DATASET
valid_scenes = valid_scenes(1:8);
%%%%%%%%%%%%%%%5

% Set which subject's labeling data to use
subStr = {'All','CW', 'LR', 'AN', 'NY','OR','CA','FL','MN','ND'}; 
subID = 1;


%% Setup neural constraints

% Set size of cropped image to feed into model
vfPatchRad = 10; % Deg
vfBuffer   = 2*vfPatchRad/degreePerPixel;

% Just use one of saved image sets to grab RF locations for plotting demo
%%% might be smarter in future to either align RFs with contour of border
%%% or center many RFs with different sizes at border location and generate
%%% underlying V1 RFs from this (depending on how you want to add in
%%% non-linearities)
NeuronMat = load([neuronRootDir,'SimData/NeuronResps_2021-06-23_test1.mat']);


%% Image selection loop

if ~exist([neuronRootDir,'NaturalImageDB/imgPairs_',fpType,'/'])
    mkdir([neuronRootDir,'NaturalImageDB/imgPairs_',fpType,'/']);
end

% Setup expected vergence distribution if we're doing that
if strcmp(fpType,'vergence')
    
    % Define bin boundaries in diopters
    dioLow  = 0;
    dioHigh = 3;
    xB      = linspace(dioLow,dioHigh,11);
    
    % Define bin centers
    x = xB(1:end-1) + diff(xB(1:2));
    
    % Set expected vergence distribution
    f = exp(x-0.5);
    f = f/sum(f);
    
    % Get quota for number of distances (diopters)
    counts = round(numSamps*f);
    
    % Get bin boundaries in meters
    bounds = 1./xB;
    
    % Get boundaries for bins
    numReps    = 1;
    nCounts    = zeros(1,numel(x));
    theseDists = nan(numSamps,1);

end

for ns = 1:numSamps
    
    %% Select a depth map
    
    % Get depth map for a scene
    scene = randsample(valid_scenes,1);
    
    depthMatrixL = load([scene_path '/lRange' sprintf('%03d',dmap_indices(scene)) '.mat']);
    depthMatrixR = load([scene_path '/rRange' sprintf('%03d',dmap_indices(scene)) '.mat']);
    
    % Get views of scene in grayscale
    imL = imread([scene_path '/lImage' sprintf('%03d',dmap_indices(scene)) '.png']);
    imR = imread([scene_path '/rImage' sprintf('%03d',dmap_indices(scene)) '.png']);
    
    imLgs = mean(imL,3);
    imRgs = mean(imR,3);
    
    % Get figure/ground map and number of labels
    if subID > 1
        subInd = subID;
        subData = load([rootDir,'labeling_data/',subStr{subInd},'_data.mat']);
    else
        subInd = randi(9,1) + 1;
        subData = load([rootDir,'labeling_data/',subStr{subInd},'_data.mat']);
    end
    fgMap = subData.data{scene};
    labelM = categories(fgMap);
    
    % This is a matrix of coordinates of every pixel, as follows
    % dim1: - is forward,     + is backward
    % dim2: - is left,        + is right
    % dim3: - is down,        + is up
    currentRangeR = depthMatrixR.rangeMap;
    currentRangeL = depthMatrixL.rangeMap;
    
    % Set out of range distances (sky, weird reflections) to nan
    currentRangeR(currentRangeR(:,:,1) > -2) = NaN;
    currentRangeL(currentRangeL(:,:,1) > -2) = NaN;
    
    % Redefine coords to
    % dim1: - is left,        + is right
    % dim2: - is down,        + is up
    % dim3: - is into screen, + is out of screen
    currentRangeXR = flipud(currentRangeR(:,:,2));
    currentRangeYR = flipud(currentRangeR(:,:,3));
    currentRangeZR = flipud(-currentRangeR(:,:,1));
    
    currentRangeXL = currentRangeL(:,:,2);
    currentRangeYL = currentRangeL(:,:,3);
    currentRangeZL = -currentRangeL(:,:,1);
    
    % Flip the figure map to match this
    fgMap = flipud(fgMap);
    
    
    %% Convert range map to disparity map for a given fixation point
    
    % Select a random fixation, keep selecting until we don't have sky,
    % have a square cyclopean window, and we're matching natural vergence
    % distribution if desired
    skyCheck       = false;
    lookDirCheck   = false;
    vergenceCheck  = false;

    
    while ~skyCheck && ~lookDirCheck && ~vergenceCheck
        
        switch fpType
            case 'rand'
                
                % row and column indices of fixation (limited to fit RF array)
                r = randi(dims(1)-vfBuffer-1) + 0.5*vfBuffer;
                c = randi(dims(2)-vfBuffer-1) + 0.5*vfBuffer;
                vergenceCheck = true;
                
            case 'vergence'
                
                % row and column indices of fixation (limited to fit RF array)
                r = randi(dims(1)-vfBuffer-1) + 0.5*vfBuffer;
                c = randi(dims(2)-vfBuffer-1) + 0.5*vfBuffer;
                
            case 'figure'
                
                % Select object at random
                objInd = randi(size(labelM,1));
                
                % Find rows/columns of points within this object
                [allrows, allcolumns] = find(fgMap == labelM{objInd,1});
                
                % Pick random location within this object for fixation
                rand_ind = randi(size(allrows,1));
                r = allrows(rand_ind);
                c = allcolumns(rand_ind);
                vergenceCheck = true;
                
            case 'border'
                
                %%% might want to make a function wrapper for all this
                %%%% NEED TO MAKE SURE THIS ALSO ONLY CONSIDERS POINTS FAR
                %%%% ENOUGH AWAY FROM IMAGE BORDERS FOR NEURON MODEL
                
                % Select object at random
                objInd = randi(size(labelM,1));
                
                % Find rows/columns of points within this object
                [allrows, allcolumns] = find(fgMap == labelM{objInd,1});
                
                % Get coords of bounding rect around object
                maxVer = max(allrows);
                minVer = min(allrows);
                maxHori = max(allcolumns);
                minHori = min(allcolumns);
                
                % Grid the box
                [row1,col1] = meshgrid([minVer maxVer],[minHori:maxHori]);
                [row2,col2] = meshgrid([minVer:maxVer],[minHori maxHori]);
                
                rows = [row1(:)' row2(:)'];
                cols = [col1(:)' col2(:)'];
                
                % Pick random location within box
                rand_ind = randi(size(rows));
                r_border = rows(rand_ind);
                c_border = cols(rand_ind);
                
                % Move point until it hits a figure/ground border
                [r_border,c_border] = find_border(r_border,c_border,...
                    minHori,maxHori,minVer,maxVer,fgMap,labelM,objInd);
                
                % Check if valid point
                %%%% this doesn't do anything, need to put in while loop
                if isnan(currentRangeR(r_border,c_border,1))
                    invalid_count = invalid_count + 1;
                    continue
                end
                
                vergenceCheck = true;
                
            otherwise
                error('Not set up yet.');
        end
        
        % Make sure left eye image can fill 600x600
        % constraint on steropair: need to choose x coord that can
        % still produce N x N px image in left eye (i.e. location
        % must be 6 deg from right border of camera FoV)
        xf = currentRangeXR(r,c)-eye_2(1);
        yf = currentRangeYR(r,c);
        zf = currentRangeZR(r,c);
        
%         thisLELookDir = 90 - atan2d(zf,(xf + eye_1(1)));
        
%         if thisLELookDir > ((hFov/2)-6)
            lookDirCheck   = true;
%         end
        
        
        % AND make sure it doesn't have sky
        if ~isnan(currentRangeZR(r,c))
            skyCheck       = true;
        end
        
        % Finally, make sure this point is keeping with our desired
        % vergence distribution
        if strcmp(fpType,'vergence')
            
            thisDist = zf;
            
            % Which bins can accept more FPs this round?
            availableBins = nCounts < counts;
            
            if thisDist >= bounds(end)
                
                lowerEdge = thisDist > bounds;
                upperEdge = thisDist < bounds;
                
                thisBin   = abs(diff(lowerEdge)) & abs(diff(upperEdge));
                
                if availableBins(thisBin)
                    
                    % update bin counter
                    nCounts(thisBin) = nCounts(thisBin) + 1;
                    
                    % exit while loop
                    vergenceCheck   = true;
                    
                end
                
            end

        end
        
    end
    
    % We get out of the loop above once our checks have passed, so let's
    % assign this as a good fixation point
    r_fix = r;
    c_fix = c;
    
    % calculate horizontal disparity of all pixels (use shifted right eye depth map for cyclopean view)
    dsp = calc_disparity([xf zf yf], [currentRangeXR(:)-eye_2(1) currentRangeZR(:) currentRangeYR(:)]);
    
%     disp = flipud(reshape(disp,dims(1),dims(2)));
    dsp = reshape(dsp,dims(1),dims(2));
    
    %% Get axes for image in degrees eccentricity (not exact, due to tangent error)
    
    haxDeg = (hax-c)*degreePerPixel;
%     vaxDeg = fliplr((vax-r)*degreePerPixel);
    vaxDeg = (vax-r)*degreePerPixel;
    
    % Crop image pair/disparity image to radius set by vfPatchRad
    hmask = abs(haxDeg) <= vfPatchRad;
    vmask = abs(vaxDeg) <= vfPatchRad;
    
    haxDegCrop = haxDeg(hmask);
    vaxDegCrop = vaxDeg(vmask);
    dispCrop = dsp(vmask,hmask);
    
    % Axes flipped between image/figure map and disparity map?
    vaxDegImg = fliplr((vax-r)*degreePerPixel);
    vmaskImg = abs(vaxDegImg) <= vfPatchRad;
    
    fgMapCrop = fgMap(vmaskImg,hmask);
   
    % binary mask containing all objects
    fgtrue = nan(size(fgMapCrop));
    fgtrue(~isundefined(fgMapCrop)) = 1;
    fgtrue(isundefined(fgMapCrop)) = 0;
    
    fgtrueFull = nan(size(fgMap));
    fgtrueFull(~isundefined(fgMap)) = 1;
    fgtrueFull(isundefined(fgMap)) = 0;
    
%     fgtrue = flipud(fgtrue);
%     fgtrueFull = flipud(fgtrueFull);
%     fgtrue = fgtrue;
    
    
    %% Generate image pair (USING RIGHT EYE AS CYCLOPEAN)
    
    % Get index of closest pixel to looking direction in each eye
    hindL = round(dims(2)/2 + (dims(2)/2)*( (90-atan2d(zf,(xf+eye_1(1))))/(hFov/2) ));
%     hindR = round(dims(2)/2 + (dims(2)/2)*( (90-atan2d(zf,(xf+eye_2(1))))/(hFov/2) ));
    hindR = c;
    
    % Crop out square of each eye's image, centered on looking direction
    haxDegL = (hax-hindL)*degreePerPixel;
    hmaskL = abs(haxDegL) <= vfPatchRad;
    
    haxDegR = (hax-hindR)*degreePerPixel;
    hmaskR = abs(haxDegR) <= vfPatchRad;
    
    % Get indices of fixation point in left eye image
    imRCrop = single(imRgs(vmaskImg,hmaskR));
    imLCrop = single(imLgs(vmaskImg,hmaskL));
    
    % Stick in struct
    imdata.imageID{ns}  = num2str(dmap_indices(scene));
    imdata.fpPx(ns,:)   = [hax(c) vax(r)];
    imdata.fpType       = fpType;
%     imdata.imR{ns}      = imRCrop;
%     imdata.imL{ns}      = imLCrop;
    imdata.dispMap{ns}  = dispCrop;
    imdata.dispUnc{ns}  = dsp;
    imdata.fgMap{ns}    = logical(fgtrue);
    imdata.fgMapUnc{ns} = logical(fgtrueFull);
    imdata.subID{ns}    = subStr{subID};
    imdata.haxDeg{ns}   = haxDeg;
    imdata.vaxDeg{ns}   = vaxDeg;
    
    imR     = imRCrop;
    imL     = imLCrop;
        
    if saveOn
    save([neuronRootDir,'NaturalImageDB/imgPairs_',fpType,'/imR',num2str(ns)],'imR');
    save([neuronRootDir,'NaturalImageDB/imgPairs_',fpType,'/imL',num2str(ns)],'imL');
    end
    
    disp(['Completed image pair ',num2str(ns),' of ',num2str(numSamps)]);
    
end


%% Save image pair and disparity map for this image patch

if saveOn
    save([neuronRootDir,'NaturalImageDB/imgPairs_',fpType,'/natImgPairs'],'imdata','-v7.3');
%     save([rootDir,'natImgPairs'],'imdata');
end


%% Get coordinates in retinal eccentricity

if plotEx
    
    %%% make this a separate function to plot from saved structure
    
    % full Img
    fig1 = figure;
    fig1.Position = [2175 300 1400 715];
    hold on;
    
    imagesc(haxDeg,vaxDeg,dsp);
%     plot([-6 -6],[-6 6],'k','linewidth',3);
%     plot([-6 6],[6 6],'k','linewidth',3);
%     plot([6 6],[-6 6],'k','linewidth',3);
%     plot([-6 6],[-6 -6],'k','linewidth',3);
    
    set(gca,'plotboxaspectratio',[1 min(dims)/max(dims) 1],'fontsize',20,...
        'xlim',[haxDeg(1) haxDeg(end)],'ylim',[vaxDeg(1) vaxDeg(end)]);
    
    cb = colorbar('Ticks',[-1 -0.5 0 0.5 1]);
    colormap( redblue );
    cb.Label.String = 'Horizontal disparity (deg)';
    cb.FontSize = 20;
    caxis([-1 1]);

    
    % RF array patch
    fig2 = figure;
    fig2.Position = [3585 300 815 715];
    hold on;
    
    imagesc(haxDegCrop,vaxDegCrop,dispCrop);
    
    if plotArr
    [fig2] = drawPopRFLocations(NeuronMat.V1Neurons,[0 0 0],3,fig2);
    [fig2] = drawPopRFLocations(NeuronMat.MTNeurons,[0 0.8 0],6,fig2);
    end
    
    set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,...
        'xlim',[-6 6],'ylim',[-6 6]);
    
    cb = colorbar('Ticks',[-1 -0.5 0 0.5 1]);
    colormap( redblue );
    cb.Label.String = 'Horizontal disparity (deg)';
    cb.FontSize = 20;
    caxis([-1 1]);
    
    % Figure-ground patch
%     fig3 = figure;
%     fig3.Position = [3585 -345 815 715];
%     hold on;
%     
%     imagesc(haxDegCrop,vaxDegCrop,fgtrue);
%     
%     if plotArr
%     [fig3] = drawPopRFLocations(NeuronMat.V1Neurons,[1 1 1],3,fig3);
%     [fig3] = drawPopRFLocations(NeuronMat.MTNeurons,[0 0.8 0],6,fig3);
%     end
%     
%     set(gca,'plotboxaspectratio',[1 1 1],'fontsize',20,...
%         'xlim',[-6 6],'ylim',[-6 6]);
%     colormap copper;
    
     % Figure-ground patch
    fig3b = figure;
    fig3b.Position = [2175 300 1400 715];
    hold on;
    
    imagesc(haxDeg,vaxDeg,fgtrueFull);
    
    set(gca,'plotboxaspectratio',[1 min(dims)/max(dims) 1],'fontsize',20,...
        'xlim',[haxDeg(1) haxDeg(end)],'ylim',[vaxDeg(1) vaxDeg(end)]);
    colormap copper;
    
    % Stereo image patch
    imLgc = imread([scene_path '/lImage' sprintf('%03d',dmap_indices(scene)) 'V.png']);
    imRgc = imread([scene_path '/rImage' sprintf('%03d',dmap_indices(scene)) 'V.png']);
    
    imLgcgs = mean(imLgc,3);
    imRgcgs = mean(imRgc,3);
    
    imRCropgc = uint8(imRgcgs(vmaskImg,hmask));
    imLCropgc = uint8(imLgcgs(vmaskImg,hmask));
    
    fig4 = figure;
    fig4.Position = [2175 -345 1400 715];
    
    subplot(1,2,1);
    imshow(imRCropgc);
    subplot(1,2,2);
    imshow(imLCropgc);
    
    % Left eye full color img
    
    fig5 = figure;
    fig5.Position = [2175 300 1400 715];
    hold on;
    
    imagesc(haxDeg,fliplr(vaxDeg),imLgc);
    set(gca,'plotboxaspectratio',[1 min(dims)/max(dims) 1],'fontsize',20,...
        'xlim',[haxDeg(1) haxDeg(end)],'ylim',[vaxDeg(1) vaxDeg(end)],'ydir','normal');
    
end
