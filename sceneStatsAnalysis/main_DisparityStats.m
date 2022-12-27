% This script generates KSDs from the V1, V2, and MT neuronal datasets,
% samples from the true and bootstrapped BORIS image sets, 

clear all;
close all;

% Define path to saved distribution data
splPath    = regexp(which('main_DisparityStats'),filesep,'split');
imStatsDir = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
saveDir    = [imStatsDir,'SceneStatsAnalysis/savedImageStats_BORISdataset/'];

addpath([imStatsDir,'sharedTools']);

% Toggle save figures on/off
saveOn = 0;


%% Generate kernel-smoothed densities from the neuronal RF centers + plot

% Generate sampling densities
makeNeuronalKSD

% Plot
KSDplots


%% Sample from true BORIS image set and bootstrapped sets
resampIter = 0;

% Get stats from V1/V2/MT from original KSDs + V1/V2 limited MT
getDispStats_KSD_BORIS('sando','all',resampIter);
getDispStats_KSD_BORIS('walking','all',resampIter);

getDispStats_KSD_BORIS('sando','ecc',resampIter);
getDispStats_KSD_BORIS('walking','ecc',resampIter);

getDispStats_KSD_BORIS('sando','vertPos',resampIter);
getDispStats_KSD_BORIS('walking','vertPos',resampIter);

% Bootstrap sampling
numReps = 100;

parfor ii = 1:numReps
    resampIter = ii;

    getDispStats_KSD_BORIS('sando','all',resampIter);
    getDispStats_KSD_BORIS('walking','all',resampIter);
    
    getDispStats_KSD_BORIS('sando','ecc',resampIter);
    getDispStats_KSD_BORIS('walking','ecc',resampIter);
    
    getDispStats_KSD_BORIS('sando','vertPos',resampIter);
    getDispStats_KSD_BORIS('walking','vertPos',resampIter);
end


%% Collect image stats

BORISbootstrap

dispDistsPlots

