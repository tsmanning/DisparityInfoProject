# DisparityInfoProject

This repository contains the analysis code for the Disparity Information Coding Project.
All analyses can be run from ./analysisPipeline.m, but this is not recommended since the
tuning curve fitting and disparity statistics estimation are best run in parallel on a
high performance computing cluster (even then it takes ~19hrs to run!).

The neural data and image dataset can be downloaded at the following location: [LINK]

Running the script included with the data will place the datasets in their expected locations
for the analysis scripts.

The various subfunctions and scripts are split into two directories: 1) sceneStatsAnalysis
(which collects the disparity statistics from the BORIS image set), and 2) V1V2MTPopulations
(which performs the analyses of the neural datasets and the comparison with the natural scene
statistics.)

## sceneStatsAnalysis

./BORISimageSet - contains the food prep and navigation disparity image sets
./savedImageStats_BORISdataset - contains saved disparity stats
./savedKSDmatFiles_BORISdataset - contains saved KSDs generated from neuronal data
./plots - holds all the scene stats plots used in the paper

main_DisparityStats.m - top level script for running the disparity statistics

makeNeuronalKSD.m - makes kernel-smoothed densities to guide disparity sampling
KSDplots.m - visualizes resulting KSDs
getDispStats.m - script run on high performance cluster to sample disparity histograms from image sets
collectBootstrapRuns.m - used to combine/sort output files from HPC
dispDistsPlots.m - plot resulting disparity histograms from analysis

## V1V2MTPopulations

./analysisFiles - holds stats/calculated tuning curves for neuronal data sets
./dataMT - raw MT neuronal data
./dataV1V2Combo - raw V1 and V2 neuronal data
./plots - holds all the neuronal plots from paper
./screen_disparity_correction - scripts used to put disparity supports from neuronal data into true disparity units

main_FitTuningCurves.m - Used to fit all neuron responses with Gabor functions
main_FitTuningCurvesCustom.m - Used to fit some edge cases
main_AssessFits.m - Calculate histograms, etc summarizing Gabor parameters
main_AnalyzeFI.m - Top level script to calculate single cell and population Fisher information
main_CompareToNDS.m - Compare population FI and disparity image stats
main_reparameterizeV1.m - Resample V1 Gabor fits using V1 or MT distributions
supp_runFanoFacCheck.m - Top level script for running supplementary spike count stats

errfun1DGabor.m - objective function for Gabor fitting
FanoFacCheck.m - Calculate spike count statistics from neuronal data
plotFigure1 - plot figure in paper
plotFigure4A - plot figure in paper
plotFigure6A - plot figure in paper
fit1DGabor.m - fits Gabor function to neuronal data
fitGeneralizedLaplacian.m - fits Laplace function to disparity histograms
loadDataMT.m - loads MT neuronal data, cleans it, corrects disparity
loadDataV1V2.m - loads V1 and V2 neuronal data, cleans it, corrects disparity
permuteFits.m - Used to resample best fit parameters in V1 data set
plot_population_FI.m - Used to calculate population FI
plot_RF_locations.m - Plot scatterplots of RF centers for each neuronal data set
subsampleCells.m - Used to define subsamples of cells depending on inclusion criteria

## sharedTools

./distributionPlot - generate fancy smoothed histograms for Figures 4/5
./pinky - script used for sampling disparities/parameters from a KSD
ColorIt.m - Define an RGB color pallet that's a little easier on the eyes
getJSDiv.m - Calculated Jensen-Shannon divergence between two distributions
shadeplot.m - Plots shaded/interpolated errorbars
uniform_sample_in_range.m - Generates a support for doing stats

## Additional notes:
This repo contains a few tools developed by others: "pinky"[^1], "distributionPlot"[^2]. See reference for a link to download file from MATLAB Central File Exchange.

## References:
[^1]: Tristan Ursell (2012). pinky (https://www.mathworks.com/matlabcentral/fileexchange/35797-generate-random-numbers-from-a-2d-discrete-distribution)
[^2]: Jonas Dorn (2017). distributionPlot (https://www.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions-distributionplot-m)
