function [dispHist1,dispHist2,dispHist3,gHist,edges_disp,f1] = getDispStats_fgBorder(imdata,subset,plotOn)

% For image set what are the disparity statistics?

%% Collect disparities within masks and histogram

res        = 50;
% res        = 25;

numIms     = numel(imdata.dispUnc);
edges_disp = linspace(-1.5,1.5,res);
edges_diff = linspace(-1.5,1.5,res);
cntr_disp  = edges_disp(1:end-1) + diff(edges_disp(1:2))/2;

dispHist1 = nan(numIms,res-1);
dispHist2 = nan(numIms,res-1);

pmax = 5;

for ii = 1:numIms
    
    % Get image size
    [imxm,imym] = meshgrid(imdata.haxDeg{ii},imdata.vaxDeg{ii});
    ecc         = sqrt(imxm.^2 + imym.^2);
    
    switch subset
        case 'all'
            mask = ones(numel(imdata.haxDeg{ii}),numel(imdata.vaxDeg{ii}),'logical');
            
        case 'figure'
            mask = imdata.fgMapUnc{ii};
            
        case 'ecc'
            %%%% doesn't match median split in MT dataset - could decrease MT
            %%%% set split level? Though it doesn't have much small
            %%%% eccentricities...
%             medEcc = median(ecc(:));
            medEcc = 7;
            
            % Central
            mask = ecc<medEcc;
            
        case 'vertPos'
            mask = imym>0;
            
    end
    
    theseDisps   = imdata.dispUnc{ii};
    theseFigures = imdata.fgMapUnc{ii};
    
    [figureDisps,groundDisps,dispDiff] = getborderDisps(theseDisps,theseFigures);
    
    % Get a 2D histogram of disparity differences for a given figure
    % disparity
    disp2DHist(ii,:,:) = histcounts2(figureDisps,dispDiff,edges_disp,edges_diff);
    
    % Ground disparity histogram
    groundHist(ii,:) = histcounts(groundDisps,edges_disp);

end


%% Plotting

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

% Collapse across images
dispHist1 = squeeze(sum(disp2DHist,1));

gHist = sum(groundHist,1);

% Normalize to get 2D PDF
dispHist1 = dispHist1/(sum(dispHist1(:))*diff(edges_disp(1:2))*diff(edges_diff(1:2)));

gHist = gHist/(sum(gHist(:))*diff(edges_disp(1:2)));

% Collapse across pedestal disparities and normalize
dispHist2 = squeeze(sum(sum(disp2DHist,1),3));
dispHist2 = dispHist2/(sum(dispHist2(:))*diff(edges_disp(1:2)));

% Collapse across differences and normalize
dispHist3 = squeeze(sum(sum(disp2DHist,1),2));
dispHist3 = dispHist3/(sum(dispHist3(:))*diff(edges_disp(1:2)));

if plotOn

    % Plot 
    f1 = figure;
    f1.Position = [300 300 720 770];
    hold on;
    
    imagesc(cntr_disp,cntr_disp,log(dispHist1));
    set(gca,'fontsize',20,'plotboxaspectratio',[1 1 1],'xlim',[cntr_disp(1) cntr_disp(end)],...
        'ylim',[cntr_disp(1) cntr_disp(end)]);
    cb = colorbar;
    xlabel('Figure disparity');
    ylabel('F-G border disparity difference');
    cb.Label.String = 'log(probability)';
    
    % Plot difference marginal distribution
    f2 = figure;
    f2.Position = [600 300 380 612];
    hold on;
    
    plot(cntr_disp,dispHist2,'k','linewidth',4);
    xlabel('F-G border disparity difference');
    ylabel('Probability density');
    set(gca,'fontsize',20,'view',[90 -90]);
    
    % Plot figure disparity marginal distribution
    f3 = figure;
    f3.Position = [300 600 627 380];
    hold on;
    
    plot(cntr_disp,dispHist3,'k','linewidth',4);
    xlabel('Figure disparity');
    ylabel('Probability density');
    set(gca,'fontsize',20);
    
    % Plot ground disparity distribution
    f4 = figure;
    f4.Position = [300 600 627 380];
    hold on;
    
    plot(cntr_disp,gHist,'k','linewidth',4);
    xlabel('Ground disparity');
    ylabel('Probability density');
    set(gca,'fontsize',20);
    
end

end

%% Helper functions

function [figureDisps,groundDisps,dispDiff] = getborderDisps(dispMap,fgmap)

% For a given disparity map, find the figure-ground borders and get both
% the change in disparity across them and the figure pedestal

% dispMap = imdata.dispUnc{imInd};
% fgmap   = imdata.fgMapUnc{imInd};

figureMask = fgmap;

% Get vertical pixels
vbordDiff   = diff(figureMask,[],1);
topbords    = [zeros(1,size(figureMask,2),'logical'); (vbordDiff > 0)];
bottombords = [(vbordDiff < 0); zeros(1,size(figureMask,2),'logical')];

% Get horizontal border pixels
hbordDiff   = diff(figureMask,[],2);
leftbords   = [zeros(size(figureMask,1),1,'logical') (hbordDiff > 0)];
rightbords  = [(hbordDiff < 0) zeros(size(figureMask,1),1,'logical')];

borderPxMask = topbords | bottombords | leftbords | rightbords;

[r,c] = find(borderPxMask == 1);

n1 = nan(numel(r),2);
n2 = nan(numel(r),2);

% For each pixel along figure border, find one point within figure and one
% point outside it, x number of pixels out along normal
for ii = 1:numel(r)
    
    % Check to make sure query points aren't out of range - if so, don't count
    % them
    if (r(ii)-15 >= 1) && (r(ii)+15 <= size(borderPxMask,1)) && (c(ii)-15 >= 1) && (c(ii)+15 <= size(borderPxMask,2))
        
        % Query logical values in 5x5 grid around current border pixel
        thisPxMat = borderPxMask(r(ii)-2:r(ii)+2,c(ii)-2:c(ii)+2);
        
        % Make the longest line segment possible between two border pixels
        % within this box
        [thisR,thisC] = find(thisPxMat == 1);
        
        [maxPairRC] = findLongLine(thisR,thisC);
        
        % Find two points along normals
        dr = diff(maxPairRC(:,1));
        dc = diff(maxPairRC(:,2));
        
        n1(ii,:) = [maxPairRC(1,1)-5*dc maxPairRC(1,2)+5*dr] + [r(ii) c(ii)];
        n2(ii,:) = [maxPairRC(2,1)+5*dc maxPairRC(2,2)-5*dr] + [r(ii) c(ii)];
        
    end

end

% Cull any points we couldn't determine points along normal for
n1NanPass = ~isnan(n1(:,1));
n2NanPass = ~isnan(n2(:,1));

% Cull any normal points that are out of range of image
n1RangePass = (n1(:,1)>1) & (n1(:,2)>1) & (n1(:,1)<size(dispMap,1)) & (n1(:,2)<size(dispMap,2));
n2RangePass = (n2(:,1)>1) & (n2(:,2)>1) & (n2(:,1)<size(dispMap,1)) & (n2(:,2)<size(dispMap,2));

passingInds = n1NanPass & n2NanPass & n1RangePass & n2RangePass;

n1 = n1(passingInds,:);
n2 = n2(passingInds,:);

% Since direction conventions change for normal points, reassign query
% points based on whether they're inside or outside of the border
n1Lin = sub2ind([size(fgmap,1) size(fgmap,2)],n1(:,1),n1(:,2));
n2Lin = sub2ind([size(fgmap,1) size(fgmap,2)],n2(:,1),n2(:,2));

n1temp = n1Lin;
n2temp = n2Lin;

% Normal 1 (Light Blue) should always be ground
flipInds1        = fgmap(n1Lin) == 1;

n1Lin(flipInds1) = n2temp(flipInds1);

[n1(:,1),n1(:,2)] = ind2sub(size(fgmap),n1Lin);

% Normal 2 (Green) should always be figure
flipInds2        = fgmap(n2Lin) ~= 1;

n2Lin(flipInds2) = n1temp(flipInds2);

[n2(:,1),n2(:,2)] = ind2sub(size(fgmap),n2Lin);

% Use these indices to select figure and ground disparities along the
% figure border

groundDisps = dispMap(n1Lin);

figureDisps = dispMap(n2Lin);

% Get paired disparity difference
dispDiff    = figureDisps - groundDisps;

end

function [maxPairRC] = findLongLine(r,c)

    numPx    = numel(r);
    startInd = 1;
    
    % Make array of unique pixel pairings
    for ii = 1:numPx-1
        
        endInd = startInd + (numPx - ii - 1);
        
        pairs(startInd:endInd,1) = ii*ones(numPx-ii,1);
        pairs(startInd:endInd,2) = ii+1:numPx;
        
        startInd = endInd + 1;
        
    end
    
    numPairs = size(pairs,1);
    
    % Calculate lengths of each pair
    for ii = 1:numPairs
        lengths(ii) = sqrt( abs(r(pairs(ii,1)) - r(pairs(ii,2)))^2 + abs(c(pairs(ii,1)) - c(pairs(ii,2)))^2 );
    end
    
    % Find pair with longest length
    [~,maxInd] = max(lengths);
    maxPair    = pairs(maxInd,:);  
    
    maxPairRC  = [r(maxPair(1)) c(maxPair(1)); r(maxPair(2)) c(maxPair(2))];

end
