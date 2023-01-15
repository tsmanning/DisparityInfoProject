% NOTE: While this script could be run all the way through, it's designed
% to be a description of the analysis pipeline for the paper. If you'd like
% to rerun, it's best to run the scripts one at a time, since the tuning
% curve fitting and disparity statistics scripts take a while!

clear all
close all

splPath = regexp(which('analysisPipeline'),filesep,'split');
topDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
addpath(genpath([topDir,'sceneStatsAnalysis']));
addpath(genpath([topDir,'V1V2MTPopulations']));
addpath(genpath([topDir,'helper_functions']));


%------------------------------(1)-------------------------------%
% This script reads in the spiking data, processes, and fits tuning
% curves The fits used for the paper are stored with
% results[area]_final.mat, if it's re-run the outputs will be results[area].mat
% - areas: flag which cell population sample we want to fit
% areas = {'V1'};
% areas = {'V2'};
% areas = {'MT'};
areas = {'V1','V2','MT'};

main_FitTuningCurves(area);

% Calls following custom functions:
% - loadDataV1V2.m
% - fit1DGabor.m
% - loadDataMT.m
% - errfun1DGabor.m
% - screen2retDisp.m
% - uniform_sample_in_range.m


%------------------------------(2)-------------------------------%
% This script enables refitting of individual neurons in cases of
% catastophic failure of the main fitting routine
main_FitTuningCurvesCustom

% Calls same functions as main_FitTuningCurves


%------------------------------(3)-------------------------------%
% This script plots the fits and assesses quality in various ways, it saves
% out a file called resultsALL_final.mat, which has all the fits plus some
% extra characterizations

%----------------%
% Options for Subsampling cells based on fit quality (r2) and
% eccentricity as controls (to make sure we are not selecting against certain
% preferences in final paper):
% - Toggle whether we want to apply all subsampling criteria applied in main
% - figures of paper
doAllSubsampling = 1;

% Define which subsample of the population we want to analyze:
subsample = 'eccentricity';
% subsample = 'fitQuality';

% Define which subsample we want: greater than or less than the threshold
% for eccentricity or fit quality:
direc  = 'lt';
% direc  = 'gt';

subsampleFlags.subsample = subsample;
subsampleFlags.direc     = direc;
%----------------%

main_AssessFits(doAllSubsampling,subsampleFlags)

% Calls following custom functions:
% - distributionPlot.m (Jonas Dorn, 2017)
% - ColorIt.m
% - subsampleCells.m


%------------------------------(4)-------------------------------%
% This script does some basic summaries of tuning characteristics and
% computes/plots FI
main_AnalyzeFI

% Calls following custom functions:
% - plot_RF_locations.m
% - plot_population_FI.m
% - subsampleCells.m


%------------------------------(5)-------------------------------%
% IN PARALLEL: This script was run on a high-performance computing cluster
% (it takes a while!) and samples disparities from the BORIS image set
% according to sampling distributions derived from neuronal RF center KSDs
main_DisparityStats

% Calls following custom functions:
% - makeNeuronalKSD.m
% - KSDplots.m
% - getDispStats.m
% - collectBootstrapRuns.m
% - dispDistsPlots.m
% - pinky.m (Tristan Ursell, 2012)


%------------------------------(6)-------------------------------%
% this script loads the FI results and compares to the natural disparity histograms
main_CompareToNDS

% Calls following custom function:
% - fitGeneralizedLaplacian


%------------------------------(7)-------------------------------%
% See which parameters of the Gabor fit can best recreate the MT FI
% distribution using the rest of the parameters drawn from V1
main_reparameterizeV1

% Calls following custom functions:
% - permuteFits.m
% - fitGeneralizedLaplacian.m
% - getJSDiv.m
% - shadeplot.m
% - pinky.m (Tristan Ursell, 2012)


% Supplementary analyses
%------------------------------(8)-------------------------------%
% Plot the distribution of spike count means/variances between areas 
% to compare "Poisson-ness"
supp_runFanoFacCheck

% Calls following custom functions:
% - FanoFacCheck.m


% Figure components
%------------------------------(9)-------------------------------%
% Plot example tuning curves used in Figures
plotFigure1
plotFigure4A
plotFigure6A

