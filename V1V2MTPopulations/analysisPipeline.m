% NOTE: While this script could be run all the way through, it's designed
% to be a description of the analysis pipeline for the paper. If you'd like
% to rerun, it's best to run the scripts one at a time, since the tuning
% curve fitting and disparity statistics scripts take a while!

splPath = regexp(which('analysisPipeline'),filesep,'split');
topDir  = [filesep,fullfile(splPath{1:numel(splPath)-2}),filesep];
addpath(genpath([topDir,'SceneStatsAnalysis/helper_functions']));


%----------------------------------------------------------------%
% This script reads in the spiking data, processes, and fits tuning
% curves The fits used for the paper are stored with
% results[area]_final.mat, if it's re-run the outputs will be results[area].mat
main_FitTuningCurves

% Calls following custom functions:
% - loadDataV1V2.m
% - fit1DGabor.m
% - loadDataMT.m


%----------------------------------------------------------------%
% This script enables refitting of individual neurons in cases of
% catastophic failure of the main fitting routine
main_FitTuningCurvesCustom


%----------------------------------------------------------------%
% This script plots the fits and assesses quality in various ways, it saves
% out a file called resultsALL_final.mat, which has all the fits plus some
% extra characterizations

% Note: includes options for subsampling cells based on fit quality (r2) and
% eccentricity as controls to make sure we are not selecting against certain
% preferences in final paper
main_AssessFits


%----------------------------------------------------------------%
% This script does some basic summaries of tuning characteristics and
% computes/plots FI
main_AnalyzeFI

% Calls following custom functions:
% - plot_RF_locations.m
% - plot_population_FI.m


%----------------------------------------------------------------%
% IN PARALLEL: This script was run on a high-performance computing cluster
% (it takes a while!) and samples disparities from the BORIS image set
% according to sampling distributions derived from neuronal RF center KSDs
main_DisparityStats


%----------------------------------------------------------------%
% this script loads the FI results and compares to the natural disparity histograms
main_CompareToNDS

% Calls following custom function:
% - fitGeneralizedLaplacian


%----------------------------------------------------------------%
% See which parameters of the Gabor fit can best recreate the MT FI
% distribution using the rest of the parameters drawn from V1
main_reparameterizeV1

% Calls following custom functions:
% - permuteFits
% - fitGeneralizedLaplacian
% - getJSDiv
