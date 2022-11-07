clear all; 
close all;

% iterations for bootstrapping
iterations = 10;

%% load fitting results
load('resultsALL_final.mat');
%V1 = load('resultsV1_final.mat');
%V2 = load('resultsV2_final.mat');
%MT = load('resultsMT_final.mat');


%% force all FI values to be positive
V1neg = sum(V1.FI(:) < 0);
V2neg = sum(V2.FI(:) < 0);
MTneg = sum(MT.FI(:) < 0);

display(['V1 removing ' num2str(V1neg) ' negative FI values']);
display(['V2 removing ' num2str(V2neg) ' negative FI values']);
display(['MT removing ' num2str(MTneg) ' negative FI values']);

V1.FI(V1.FI < 0) = 0;
V2.FI(V2.FI < 0) = 0;
MT.FI(MT.FI < 0) = 0;


%% 1) plot and analyze full dataset
plot_RF_locations(V1,V2,MT,'');
plot_population_FI(V1,V2,MT,'',iterations);


%% 2) only keep neurons with R2 > 0.75
flag = '_R2morethanPt75';

V1R2 = V1;
V2R2 = V2;
MTR2 = MT;

V1_inds = find(V1.r2 > 0.75);
V2_inds = find(V2.r2 > 0.75);
MT_inds = find(MT.r2 > 0.75);

display(['R2 .75 V1 =  ' num2str(numel(V1_inds))]);
display(['R2 .75 V2 =  ' num2str(numel(V2_inds))]);
display(['R2 .75 MT =  ' num2str(numel(MT_inds))]);

V1R2.x_pos  = V1.x_pos(V1_inds);
V1R2.y_pos  = V1.y_pos(V1_inds);
V1R2.FI     = V1.FI(V1_inds,:);
V1R2.allresp = V1.allresp(V1_inds,:);

V2R2.x_pos = V2.x_pos(V2_inds);
V2R2.y_pos = V2.y_pos(V2_inds);
V2R2.FI = V2.FI(V2_inds,:);
V2R2.allresp = V2.allresp(V2_inds,:);

MTR2.x_pos = MT.x_pos(MT_inds);
MTR2.y_pos = MT.y_pos(MT_inds);
MTR2.FI = MT.FI(MT_inds,:);
MTR2.allresp = MT.allresp(MT_inds,:);

plot_RF_locations(V1R2,V2R2,MTR2,flag);
plot_population_FI(V1R2,V2R2,MTR2,flag,iterations);


%% 3) crop everything to 10 deg and resample MT to same elevation/abs azimuth as V1/V2
flag = '_resampledEccentricities';

% grab VF subset of MT
% get indices of MT RF's that overlap with V1/V2; allow points to be either
% left or right of fixation, at same eccentricities
% in practice, this is just the V1 range because all V2 is within V2
V1V2_x_range = max(abs([V1.x_pos V2.x_pos])); % don't care about left/right, just the max x-coordinate
V1V2_y_range = [min([V1.y_pos V2.y_pos]) max([V1.y_pos V2.y_pos])];
MT.overlapping_RF_inds = abs(MT.x_pos) <= V1V2_x_range ...
    & MT.y_pos >= V1V2_y_range(1) & MT.y_pos <= V1V2_y_range(2);

MT.FI_overlap = MT.FI(MT.overlapping_RF_inds,:);

MT_resam.x_pos = MT.x_pos(MT.overlapping_RF_inds);
MT_resam.y_pos = MT.y_pos(MT.overlapping_RF_inds);
MT_resam.FI    = MT.FI(MT.overlapping_RF_inds,:);
MT_resam.allresp    = MT.allresp(MT.overlapping_RF_inds,:);

% placeholder
MT_resam.FI_overlap = MT.FI_overlap;

% resampled V1 and MT
% plot_RF_locations(V1_resam,V2_resam,MT_resam,flag); % Not sure where
% V1/V2_resam come from or why they're necessary, so just changed to last
% subsampling
% plot_population_FI(V1_resam,V2_resam,MT_resam,flag);
plot_RF_locations(V1R2,V2R2,MT_resam,flag);
plot_population_FI(V1R2,V2R2,MT_resam,flag,iterations);


%% 4) subsample neurons so VF sampling is uniform *****[TO DO]

% keyboard


%% 5) replot FI with 10 most informative neurons removed
% info_V1 = sum(V1.FI,2);
% info_V2 = sum(V2.FI,2);
% info_MT = sum(MT.FI,2);
% info_MT_overlap = sum(MT.FI_overlap,2);
% 
% [B,I] = sort(info_V1,'descend');
% most_info_V1 = I(1:10);
% 
% [B,I] = sort(info_V2,'descend');
% most_info_V2 = I(1:10);
% 
% [B,I] = sort(info_MT,'descend');
% most_info_MT = I(1:10);
% 
% [B,I] = sort(info_MT_overlap,'descend');
% most_info_MT_overlap = I(1:10);
% 
% V1.FI(most_info_V1,:) = [];
% V2.FI(most_info_V2,:) = [];
% MT.FI(most_info_MT,:) = [];
% MT.FI_overlap(most_info_MT_overlap,:) = [];
% 
% plot_population_FI(V1,V2,MT,'minus10');
