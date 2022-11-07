% McCann analyses to run
% load('/Users/tylermanning/Documents/MATLAB/cooperlab/2-Modeling_simulations/disparityMTModel/NaturalImageDB/imgPairs_rand/natImgPairs.mat')
% [edges_disp,f1,dispHistV1,dispHistV2,dispHistMT,dispHistCirc]  = getDispStats_KSD(imdata,'all',1);
% [edges_disp,f1,dispHistV1_m1,dispHistV2_m1,dispHistMT_m1,dispHistCirc_m1,...
%     dispHistV1_m2,dispHistV2_m2,dispHistMT_m2,dispHistCirc_m2] = getDispStats_KSD(imdata,'ecc',1);
% [edges_disp,f1,dispHistV1_m1,dispHistV2_m1,dispHistMT_m1,dispHistCirc_m1,...
%     dispHistV1_m2,dispHistV2_m2,dispHistMT_m2,dispHistCirc_m2] = getDispStats_KSD(imdata,'vertPos',1);

numCores = feature('numcores');
nCores = numCores;
maxNumCompThreads(nCores);
parpool(nCores);

splPath  = regexp(which('imStats9Mar2022'),filesep,'split');
rootDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];

addpath([rootDir,'1-baseCode']);
addpath([rootDir,'3-sharedTools']);
addpath([rootDir,'5-sceneStatsAnalysis']);
    
% BORIS analyses to run
numReps = 100;

parfor ii = 1:numReps
    resampIter = ii;
    
    getDispStats_KSD_BORIS('sando','all',0,resampIter);
    getDispStats_KSD_BORIS('walking','all',0,resampIter);
    
    getDispStats_KSD_BORIS('sando','ecc',0,resampIter);
    getDispStats_KSD_BORIS('walking','ecc',0,resampIter);
    
    getDispStats_KSD_BORIS('sando','vertPos',0,resampIter);
    getDispStats_KSD_BORIS('walking','vertPos',0,resampIter);
end