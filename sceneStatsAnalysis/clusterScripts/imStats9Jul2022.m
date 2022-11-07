% McCann analyses to run
splPath  = regexp(which('imStats9Jul2022'),filesep,'split');
rootDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];

addpath([rootDir,'1-baseCode']);
addpath([rootDir,'3-sharedTools']);
addpath([rootDir,'5-sceneStatsAnalysis']);

% BORIS analyses to run
resampIter = 0;

getDispStats_KSD_BORIS_rect('sando','all',0,resampIter);
getDispStats_KSD_BORIS_rect('walking','all',0,resampIter);

getDispStats_KSD_BORIS_rect('sando','ecc',0,resampIter);
getDispStats_KSD_BORIS_rect('walking','ecc',0,resampIter);

getDispStats_KSD_BORIS_rect('sando','vertPos',0,resampIter);
getDispStats_KSD_BORIS_rect('walking','vertPos',0,resampIter);
