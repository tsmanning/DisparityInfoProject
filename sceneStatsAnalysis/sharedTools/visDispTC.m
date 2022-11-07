function [fig1] = visDispTC(neuronDataV1,neuronDataMT,neuronResps,tuningPars,neuronInd,cellInds)

% Plot datapoints and fitted tuning curve for a given neuron in pop
% Note: tuning curve indices may not line with responses if only a subset
% of cells were fit.

%% Organize input data & define some things

tuningFxn = @(disp,r0,A,prPosDisp,sig,prFreq,prPhase) ...
            r0 + A*exp(-0.5*((disp - prPosDisp)/sig).^2 ).*...
                   cos(2*pi*prFreq*(disp - prPosDisp) + prPhase);

[numDataDisps,numReps] = size(neuronResps);
numTCDisps  = size(tuningPars,2);

% Subset of V1 cell locations that were actually fit
v1subset    = cellInds{1};
v1subset    = unique(v1subset(:,4));

dataDisps   = linspace(neuronDataV1.PrefHDisp(1),neuronDataV1.PrefHDisp(end),numDataDisps);
TCDisps     = linspace(neuronDataV1.PrefHDisp(1),neuronDataV1.PrefHDisp(end),61);

xpos = repmat(dataDisps',[1 numReps]);

if numel(neuronInd) > 3
    
    respMatsz = size(neuronResps{1,1});
    respMatsz(4) = numel(v1subset);
    
    respData = cellfun(@(x) x(neuronInd(1),neuronInd(2),neuronInd(3),v1subset(neuronInd(4))),neuronResps);
    linInd = sub2ind(respMatsz,neuronInd(1),neuronInd(2),neuronInd(3),neuronInd(4));
else
    respData = cellfun(@(x) x(neuronInd(1),neuronInd(2),neuronInd(3)),neuronResps);
    linInd = sub2ind(size(neuronResps{1,1}),neuronInd(1),neuronInd(2),neuronInd(3));
end


%% Plot
fig1 = figure;
fig1.Position = [100 100 625 575];
hold on;

% Plot trialwise datapoints
scatter(xpos(:),respData(:),50,[0.7 0.7 0.7],'filled');

% Plot response means for each disparity
scatter(dataDisps,mean(respData,2),125,'k','filled');

% Plot TC
tp     = tuningPars(linInd,:);
thisTC = tuningFxn(TCDisps,tp(1),tp(2),tp(3),tp(4),tp(5),tp(6));
plot(TCDisps,thisTC,'b','linewidth',4);

% Plot line for preferred disparity
estPD = findPrefDisp(tp);
plot(estPD*[1 1],[0 1],'--r','linewidth',4);

% Labeling
if numel(neuronInd) == 4
    disp(['Pref HD = ',num2str(neuronDataV1.PrefHDisp(neuronInd(2))),...
           '; Ph Disp = ',num2str(round(neuronDataV1.PhaseDisp(neuronInd(3)),2)),'\circ',...
           '; SF = ',num2str(round(neuronDataV1.SF(neuronInd(4)),2)),'cyc/\circ',...
           '; Ori = ',num2str(round(neuronDataV1.Ori(neuronInd(1)),2)),'\circ']);
end

set(gca,'fontsize',13,'xtick',dataDisps,'ytick',[0 1],'ylim',[0 1]);
% set(gca,'fontsize',13,'xtick',dataDisps,'ytick',[-1 0 1],'ylim',[-1 1]);
ylabel('Response');
xlabel('Horizontal Disparity (\circ)');

% Draw which RF location cell was drawn from
fig2 = drawPopRFLocations(neuronDataV1,[0 0 0],1);
fig2.Position = [800 100 625 575];
hold on;
if numel(neuronInd) == 4
    drawEmptyCirc(neuronDataV1.xpos(v1subset(neuronInd(4))),...
                  neuronDataV1.ypos(v1subset(neuronInd(4))),...
                  neuronDataV1.RFrad(v1subset(neuronInd(4))),[1 0 0],2);
else
    drawEmptyCirc(neuronDataMT.xpos(neuronInd(3)),...
                  neuronDataMT.ypos(neuronInd(3)),...
                  neuronDataMT.RFrad(neuronInd(3)),[1 0 0],2);
end

end