function [theseDisps,thisDispR] = grabRFDisps(MTNeurons,imdata)

% For each RF location, get the pixelwise disp range

numIms = numel(imdata.imageID);
numRFs = numel(MTNeurons.xpos);

%%% assume -6:6deg for now, fix later
xDeg = linspace(-6,6,601);
yDeg = linspace(-6,6,601);

[xd2,yd2] = meshgrid(xDeg,yDeg);

theseDisps = cell(numIms,numRFs);
thisDispR(1).disps = [];
thisDispR(2).disps = [];
thisDispR(3).disps = [];

for ii = 1:numIms

    thisIm = imdata.dispMap{ii};
    
    for jj = 1:numRFs
        
        thisXval = MTNeurons.xpos(jj);
        thisYval = MTNeurons.ypos(jj);
        thisrad  = MTNeurons.RFrad(jj);
        
        RFmask = sqrt((xd2 - thisXval).^2 + (yd2 - thisYval)).^2 <= thisrad;
        
        theseDisps{ii,jj} = thisIm(RFmask);
        
        thisDispR(jj).disps = [thisDispR(jj).disps thisIm(RFmask)];
        
    end

end


if 0 
    
   figure;
   subplot(3,1,1);
   histogram(thisDispR(1).disps);
   subplot(3,1,2);
   histogram(thisDispR(2).disps);
   subplot(3,1,3);
   histogram(thisDispR(3).disps);
   
   [r1,e1] = histcounts(thisDispR(1).disps,21);
   r1 = r1/sum(r1);
   [r2,e2] = histcounts(thisDispR(2).disps,21);
   r2 = r2/sum(r2);
   [r3,e3] = histcounts(thisDispR(3).disps,21);
   r3 = r3/sum(r3);
   
   x1 = e1(1:end-1) + diff(e1);
   x2 = e2(1:end-1) + diff(e2);
   x3 = e3(1:end-1) + diff(e3);
   
   figure;
   hold on;
   bar(x1,r1,'facecolor','r','facealpha',0.5);
   bar(x2,r2,'facecolor','g','facealpha',0.5);
   bar(x3,r3,'facecolor','b','facealpha',0.5);
   
end

end