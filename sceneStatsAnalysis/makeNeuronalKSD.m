%% Make kernel-smoothed densities

clear all
close all

% Define dirs
splPath = regexp(which('makeNeuronalKSD'),filesep,'split');
topDir  = [filesep,fullfile(splPath{1:numel(splPath)-2}),filesep];
saveDir = [topDir,filesep,'SceneStatsAnalysis/savedKSDmatFiles_BORISdataset/'];
analyDir = [topDir,filesep,'V1V2MTPopulations/analysisFiles/'];

% Toggle saving KSD mat files
saveOn = 1;


%% Grab all RF x/y locations
% load fitting results
load([analyDir,'resultsALL_final.mat']);
load([analyDir,'eccInds.mat']);
load([analyDir,'r2Inds.mat']);

% Eccentricity + fit limited subsample
V1a.x_pos = V1.x_pos(eccInds.V1);
V1a.y_pos = V1.y_pos(eccInds.V1);
V2a.x_pos = V2.x_pos(eccInds.V2);
V2a.y_pos = V2.y_pos(eccInds.V2);

% Eccentricity + fit limited + V1/V2-overlapping subsample
MTa.x_posBnd = MT.x_pos(eccInds.MT);
MTa.y_posBnd = MT.y_pos(eccInds.MT);

% Eccentricity + fit limited subsample
MTa.x_pos = MT.x_pos(r2Inds.MT);
MTa.y_pos = MT.y_pos(r2Inds.MT);

incInds = sqrt(MTa.x_pos.^2 + MTa.y_pos.^2) <= 10;

MTa.x_pos = MTa.x_pos(incInds);
MTa.y_pos = MTa.y_pos(incInds);

% Overwrite structs to keep only x/y pos to save with the KSDs
V1 = V1a;

% Eccentricity limited + fit limited cell samples
V1 = V1a;
V2 = V2a;
MT = MTa;


%% Make KSD plots from each area

rangeMax = 10;         % Limit sampling window to +/-10deg around fixation point
dx       = 10/103;     % Deg (from BORIS dataset)

supp1D   = -rangeMax:dx:rangeMax;
[gx,gy]  = meshgrid(supp1D,supp1D);
gxL      = gx(:);
gyL      = gy(:);
suppSz   = size(gx,1);

% V1 density plot
V1density = ksdensity([V1a.x_pos' V1a.y_pos'],[gxL gyL]);

V1densityMat = reshape(V1density,[suppSz suppSz]);

% V2 density plot
V2density = ksdensity([V2a.x_pos' V2a.y_pos'],[gxL gyL]);

V2densityMat = reshape(V2density,[suppSz suppSz]);

% MT density plot
MTdensity = ksdensity([MTa.x_pos' MTa.y_pos'],[gxL gyL]);

MTdensityMat = reshape(MTdensity,[suppSz suppSz]);

% Resampled MT density plot
V1V2Rectdensity = ksdensity([MTa.x_posBnd' MTa.y_posBnd'],[gxL gyL]);

V1V2RectMat = reshape(V1V2Rectdensity,[suppSz suppSz]);

% Circular center 10 deg
CircDensity = sqrt(gx.^2 + gy.^2)<=10;

circDensityMat = CircDensity/sum(CircDensity(:));


%% Save these to plug into image stats script
if saveOn
    save([saveDir,'V1densityMat_BORIS.mat'],'V1densityMat','V1');
    save([saveDir,'V2densityMat_BORIS.mat'],'V2densityMat','V2');
    save([saveDir,'MTdensityMat_BORIS.mat'],'MTdensityMat','MT');
    save([saveDir,'V1V2Rect_BORIS.mat'],'V1V2RectMat');
    save([saveDir,'circDensityMat_BORIS.mat'],'circDensityMat');
end

