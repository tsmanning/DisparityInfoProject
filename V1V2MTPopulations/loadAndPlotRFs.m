clear all; close all;

% load fitting results
V1 = load('resultsV1.mat');
V2 = load('resultsV2.mat');
MT = load('resultsMT.mat');

%% LOAD ADDITIONAL MT/V1 DATA AND APPEND TO DATA STRUCTURE

% load MT metadata
T_MT = readtable('metadataMT.xlsx');

% for each file name in experiments, find index from xls and get other data
for f = 1:length(MT.experiments)
    
    % get filename
    fname = MT.experiments{f}.fn;
    fname = fliplr(strtok(fliplr(fname),'_'));
    fname = fname(1:end-4);

    % find row in Greg's files
    this_row = find(strcmp(T_MT.FIleName,fname));
    
    % grab that from that row and assign it to this neuron
    % positive x is right in VF, positive y is up in VF
    MT.x_pos(f)             = T_MT.RF_x(this_row); % note this does not account for vergence, see Jenny's notes
    MT.y_pos(f)             = T_MT.RF_y(this_row);
    
end

% grab V1 RF location data, much simpler
for v = 1:length(V1.experiments)
    
    % positive x is to the animals left in their visual field so we flip to match MT data, positive y is up
    V1.x_pos(v) = -V1.experiments{v}.x_pos;
    V1.y_pos(v) = V1.experiments{v}.y_pos;
    
end

% grab V2 RF location data, much simpler
for v = 1:length(V2.experiments)
    
    % positive x is to the animals left in their visual field so we flip to match MT data, positive y is up
    V2.x_pos(v) = -V2.experiments{v}.x_pos;
    V2.y_pos(v) = V2.experiments{v}.y_pos;
    
end

%% RECEPTIVE FIELD LOCATIONS

figure; hold on; title('Receptive field locations');
scatter(V1.x_pos,V1.y_pos,[],ColorIt('b'),'filled');
scatter(V2.x_pos,V2.y_pos,[],ColorIt('g'),'filled');
scatter(MT.x_pos,MT.y_pos,[],ColorIt('r'),'filled');
plot([-25 25],[0 0],'k:');
plot([0 0],[-25 25],'k:');
axis([-25 25 -25 25]); axis equal tight;
legend('V1','MT-overlap','MT-other');
xlabel('horizontal eccentricity (deg)'); ylabel('vertical eccentricity (deg)'); box on;
saveas(gcf,['V1_V2_MT_RFlocations.png']); 

%% FISHER INFORMATION

% visualize population FI
x = [-2 : 0.01 : 2];

v1 = sum(V1.FI); v1 = v1 / sum(v1);
v2 = sum(V2.FI); v2 = v2 / sum(v2);
mt = sum(MT.FI); mt = mt / sum(mt);

figure; hold on; title('Raw FI');
plot(x, v1,'color',ColorIt('b'),'linewidth',2);
plot(x, v2,'color',ColorIt('g'),'linewidth',2);
plot(x, mt,'color',ColorIt('r'),'linewidth',2);
legend('V1','V2','MT');
axis square; box on;
xlabel('horizontal disparity (deg)'); ylabel('normalized FI');
saveas(gcf,['./plots/FI_Raw.png']); 



