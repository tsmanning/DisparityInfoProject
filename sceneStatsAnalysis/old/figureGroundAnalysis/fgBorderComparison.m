
% load('/media/tyler/Data/MATLAB/cooperLab/2-Modeling_Simulations/disparityMTModel/NaturalImageDB/imgPairs_rand/natImgPairs')

imInd = 30; % WORKS
% imInd = 3;

dispMap = imdata.dispUnc{imInd};
fgmap   = imdata.fgMapUnc{imInd};

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
n1RangePass = (n1(:,1)>1) & (n1(:,2)>1) & (n1(:,1)<size(dispMap,1)) & (n2(:,2)<size(dispMap,2));
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


%% Plot
if 1
    
close all

% figure;
% imagesc(borderPxMask);

[r,c] = find(borderPxMask==1);

redblue = make_red_blue_colormap(1);

f = figure;
f.Position = [100 100 1500 750];
hold on;
imagesc(dispMap);
colormap(redblue);
scatter(c,r,6,[1 1 0]);
scatter(n1(:,2),n1(:,1),6,[0 1 1]);
scatter(n2(:,2),n2(:,1),6,[0 1 0]);
set(gca,'xlim',[0 1920],'ylim',[0 1080],'plotboxaspectratio',[1.78 1 1],'visible','off');

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