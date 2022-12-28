% This script generates KSDs from the V1, V2, and MT neuronal datasets,
% samples from the true and bootstrapped BORIS image sets, 

% NOTE: the somewhat weird design of the pipeline here is due to running it
% on a high-performance computing cluster with lots of nodes where one of
% the parallel threads might fail. There's a lot of number crunching here
% that required ~a day on the cluster: keep in mind if trying to bootstrap
% on your own machine.

clear all;
close all;

% Define path to saved distribution data
splPath    = regexp(which('main_DisparityStats'),filesep,'split');
rootDir    = [filesep,fullfile(splPath{1:numel(splPath)-2}),filesep];
imStatsDir = [rootDir,'SceneStatsAnalysis',filesep];
dataDir    = [rootDir,'V1V2MTPopulations',filesep];
saveDir    = [imStatsDir,'SceneStatsAnalysis',filesep,'savedImageStats_BORISdataset',filesep];

addpath([imStatsDir,'sharedTools']);
addpath([imStatsDir,'BORISimageSet']);
addpath([dataDir,'analysisFiles']);

% Toggle save figures on/off
saveOn = 0;

% Bootstrap sampling
numReps = 1;
% numReps = 100;

% Parallel computing on (not recommended to turn off unless testing)
parallelOn = 1;

% Define number of histogram bin edges
res = 52;

% Define disparity boundaries of histogram (deg)
ub = 2;
lb = -2;

% Define bin edges and centers
edges_disp = linspace(lb,ub,res);
cntr_disp  = edges_disp(1:end-1) + diff(edges_disp(1:2))/2;


%% Generate kernel-smoothed densities from the neuronal RF centers + plot

% Generate sampling densities
makeNeuronalKSD

% Plot
KSDplots


%% Sample from true BORIS image set and bootstrapped sets

% NOTE: although the functionality is there, the comparisons of different
% VF regions (parafoveal vs eccentric, upper vs lower) were not included
% in the paper and are commented out here. See Sprague & Cooper et al.
% 2015 for some of those comparisons using 10deg uniform sampling.

% Get stats from V1/V2/MT from original KSDs + V1/V2 limited MT
resampIter = 0;

getDispStats('sando','all',res,resampIter);
getDispStats('walking','all',res,resampIter);

% getDispStats('sando','ecc',res,resampIter);
% getDispStats('walking','ecc',res,resampIter);
% 
% getDispStats('sando','vertPos',res,resampIter);
% getDispStats('walking','vertPos',res,resampIter);

% Setup parallel processing workers
if parallelOn
    numCores = feature('numcores');
    nCores   = numCores;
    maxNumCompThreads(nCores);
    parpool(nCores);

    parfor ii = 1:numReps
        resampIter = ii;

        getDispStats('sando','all',res,resampIter);
        getDispStats('walking','all',res,resampIter);

        %     getDispStats('sando','ecc',res,resampIter);
        %     getDispStats('walking','ecc',res,resampIter);
        %
        %     getDispStats('sando','vertPos',res,resampIter);
        %     getDispStats('walking','vertPos',res,resampIter);
    end
else
    for ii = 1:numReps
        resampIter = ii;

        getDispStats('sando','all',res,resampIter);
        getDispStats('walking','all',res,resampIter);

        %     getDispStats('sando','ecc',res,resampIter);
        %     getDispStats('walking','ecc',res,resampIter);
        %
        %     getDispStats('sando','vertPos',res,resampIter);
        %     getDispStats('walking','vertPos',res,resampIter);
    end
end


%% Collect image stats

numHistBins = res-1;

% Collect individual bootstrap stat runs from HPC into a single structure
collectBootstrapRuns(numReps,numHistBins);

% Generate plots from the collected stats
dispDistsPlots

