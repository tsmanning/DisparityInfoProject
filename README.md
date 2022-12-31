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

### ./BORISimageSet -
### ./savedImageStats_BORISdataset -
### ./savedKSDmatFiles_BORIDdataset -
### ./plots -

### main_DisparityStats.m -

### makeNeuronalKSD.m -
### KSDplots.m -
### getDispStats.m -
### collectBootstrapRuns.m -
### dispDistsPlots.m -

## V1V2MTPopulations

### ./analysisFiles -
### ./dataMT -
### ./dataV1V2Combo -
### ./plots -
### ./screen_disparity_correction -

### main_FitTuningCurves.m -
### main_FitTuningCurvesCustom.m -
### main_AssessFits.m -
### main_AnalyzeFI.m -
### main_CompareToNDS.m -
### main_reparameterizeV1.m -
### supp_runFanoFacCheck.m -

### errfun1DGabor.m -
### FanoFacCheck.m -
### fig1c_demo_infodiscrimax.m -
### fit1DGabor.m -
### fitGeneralizedLaplacian.m -
### loadDataMT.m -
### loadDataV1V2.m -
### permuteFits.m -
### plotExampleTC.m -
### plot_population_FI.m -
### plot_RF_locations.m -
### subsampleCells.m -

## sharedTools

### ./distributionPlot - 
### ./pinky -
### ColorIt.m -
### getJSDiv.m -
### shadeplot.m -
### uniform_sample_in_range.m -
