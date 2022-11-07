% McCann analyses to run

numCores = feature('numcores');
nCores = numCores;
maxNumCompThreads(nCores);
parpool(nCores);

splPath  = regexp(which('imStats11Apr2022A'),filesep,'split');
rootDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];

addpath([rootDir,'1-baseCode']);
addpath([rootDir,'3-sharedTools']);
addpath([rootDir,'5-sceneStatsAnalysis']);
    
% BORIS analyses to run
numReps = 50;

parfor ii = 1:numReps
    resampIter = ii;
    
    getDispStats_KSD_BORIS_rect('sando','all',0,resampIter);
    getDispStats_KSD_BORIS_rect('walking','all',0,resampIter);
    
    getDispStats_KSD_BORIS_rect('sando','ecc',0,resampIter);
    getDispStats_KSD_BORIS_rect('walking','ecc',0,resampIter);
    
    getDispStats_KSD_BORIS_rect('sando','vertPos',0,resampIter);
    getDispStats_KSD_BORIS_rect('walking','vertPos',0,resampIter);
end