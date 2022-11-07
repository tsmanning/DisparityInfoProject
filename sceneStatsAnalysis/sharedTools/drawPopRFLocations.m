function [fig] = drawPopRFLocations(neuronMat,linecolor,linewidth,varargin)

% Simple function to draw RF locations and extents (out to 2SD of gaussian envelope)
% for a neural population (V1 or MT)

numRFs = numel(neuronMat.RFrad);

xpos = neuronMat.xpos;
ypos = neuronMat.ypos;
rads = neuronMat.RFrad;

thisColor = linecolor;
thisWidth = linewidth;

xlim = ceil(neuronMat.vf(2))*[-1 1];
ylim = xlim;

if ~isempty(varargin)
    fig = varargin{1};
else
    fig = figure;
    fig.Position = [100 100 650 600];
end
hold on;

for ii = 1:numRFs

    drawEmptyCirc(xpos(ii),ypos(ii),rads(ii),thisColor,thisWidth);
    
end
set(gca,'xlim',xlim,'xtick',[xlim(1) 0 xlim(2)],'ylim',ylim,'ytick',[xlim(1) 0 xlim(2)],...
    'plotboxaspectratio',[1 1 1],'fontsize',20);
xlabel('Horizontal position (\circ)');
ylabel('Vertical position (\circ)');

end
