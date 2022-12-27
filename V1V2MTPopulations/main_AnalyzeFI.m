% This script analyzes (plots FI distribution and RF centers for) the FI distributions 
% from each of the cortical areas

% NOTE: each cell subsampling step below includes the previous steps as well

clear all; 
close all;

% Iterations for bootstrapping
iterations = 250;

% Define path to saved tuning curve data
splPath = regexp(which('main_AnalyzeFI'),filesep,'split');
topDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
dataDir = [topDir,filesep,'analysisFiles',filesep];


%% Load fitting results
load([dataDir,'resultsALL_final.mat']);


%% Force all FI values to be positive
V1neg = sum(V1.FI(:) < 0);
V2neg = sum(V2.FI(:) < 0);
MTneg = sum(MT.FI(:) < 0);

display(['V1 removing ' num2str(V1neg) ' negative FI values']);
display(['V2 removing ' num2str(V2neg) ' negative FI values']);
display(['MT removing ' num2str(MTneg) ' negative FI values']);

V1.FI(V1.FI < 0) = 0;
V2.FI(V2.FI < 0) = 0;
MT.FI(MT.FI < 0) = 0;


%% 1) Plot and analyze full dataset
% plot_RF_locations(V1,V2,MT,'',topDir);
% plot_population_FI(V1,V2,MT,'',iterations,topDir);


%% 2) Only keep neurons with R2 > 0.75
[r2Inds,eccInds] = subsampleCells();

flag = '_R2morethanPt75';

V1R2 = V1;
V2R2 = V2;
MTR2 = MT;

V1R2.x_pos  = V1.x_pos(r2Inds.V1);
V1R2.y_pos  = V1.y_pos(r2Inds.V1);
V1R2.FI     = V1.FI(r2Inds.V1,:);
V1R2.allresp = V1.allresp(r2Inds.V1,:);

V2R2.x_pos = V2.x_pos(r2Inds.V2);
V2R2.y_pos = V2.y_pos(r2Inds.V2);
V2R2.FI = V2.FI(r2Inds.V2,:);
V2R2.allresp = V2.allresp(r2Inds.V2,:);

MTR2.x_pos = MT.x_pos(r2Inds.MT);
MTR2.y_pos = MT.y_pos(r2Inds.MT);
MTR2.FI = MT.FI(r2Inds.MT,:);
MTR2.allresp = MT.allresp(r2Inds.MT,:);

% plot_RF_locations(V1R2,V2R2,MTR2,flag,topDir);
% plot_population_FI(V1R2,V2R2,MTR2,flag,iterations,topDir);


%% 3) Crop everything to 10 deg and resample MT to same elevation/abs azimuth as V1/V2
flag = '_resampledEccentricities';

V1ecc = V1;
V2ecc = V2;
MTecc = MT;

V1ecc.x_pos   = V1.x_pos(eccInds.V1);
V1ecc.y_pos   = V1.y_pos(eccInds.V1);
V1ecc.FI      = V1.FI(eccInds.V1,:);
V1ecc.allresp = V1.allresp(eccInds.V1,:);

V2ecc.x_pos   = V2.x_pos(eccInds.V2);
V2ecc.y_pos   = V2.y_pos(eccInds.V2);
V2ecc.FI      = V2.FI(eccInds.V2,:);
V2ecc.allresp = V2.allresp(eccInds.V2,:);

MTecc.x_pos   = MT.x_pos(eccInds.MT);
MTecc.y_pos   = MT.y_pos(eccInds.MT);
MTecc.FI      = MT.FI(eccInds.MT,:);
MTecc.allresp = MT.allresp(eccInds.MT,:);

% Plot resampled population
plot_RF_locations(V1ecc,V2ecc,MTecc,flag,topDir);
plot_population_FI(V1ecc,V2ecc,MTecc,flag,iterations,topDir);


