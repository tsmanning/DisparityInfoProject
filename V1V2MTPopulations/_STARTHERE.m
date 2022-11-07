% This script reads in the spiking data, processes, and fits tuning
% curves. No need to re-run. The fits used for the paper are stored with
% results[area]_final.mat, if it's re-run the outputs will be results[area].mat
main_FitTuningCurves

% This script enables refitting of individual neurons in cases of
% catastophic failure of the main fitting routine
main_FitTuningCurvesCustom

% This script plots the fits and assesses quality in various ways, it saves
% out a file called resultsALL_final.mat, which has all the fits plus some
% extra characterizations
main_AssessFits

% This script ddoes some basic summaries of tuning characteristics and
% computes/plots FI
main_AnalyzeFI

% this script loads the FI results and compares to the natural disparity histograms
main_CompareToNDS