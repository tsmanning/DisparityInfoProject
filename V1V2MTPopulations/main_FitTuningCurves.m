function [] = main_FitTuningCurves(areas)

close all;

splPath  = regexp(which('main_FitTuningCurves.m'),filesep,'split');
topDir   = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
saveDir  = [topDir,'analysisFiles',filesep];

addpath([topDir,filesep,'helper_functions']);


% flag indicating whether to correct the stimulus disparities to account
% for planar screen (horopter deviation and foreshortening)
correct_screen_disparity = 1;

% V1
if any(contains(areas,'V1'))

    % load and preprocess data
    experiments = loadDataV1V2(1,correct_screen_disparity);

    FI    = []; % Fisher information (for fit1DGabor)
    P     = []; % parameters
    S     = []; % interpolated spike rates
    X     = []; % horizontal disparity
    E     = []; % fitting error

   
    for n = 1 : length(experiments)
        dat = experiments{n}.dat; % pre-processed data

        figure(1); clf;
        [FI,P,S,X,E] = fit1DGabor( 'V1', dat, FI, P, S, X, E, n, length(experiments) );
    end

    if correct_screen_disparity
        save( [saveDir,'resultsV1.mat'], 'FI', 'P', 'S', 'X', 'E', 'experiments' );
    else
        save( [saveDir,'resultsV1_no correction.mat'], 'FI', 'P', 'S', 'X', 'E', 'experiments' );
    end
end

% V2
if any(contains(areas,'V2'))

    % load and preprocess data
    experiments = loadDataV1V2(2,correct_screen_disparity);

    FI    = []; % Fisher information (for fit1DGabor)
    P     = []; % parameters
    S     = []; % interpolated spike rates
    X     = []; % horizontal disparity
    E     = []; % fitting error

    for n = 1 : length(experiments)
        dat = experiments{n}.dat; % pre-processed data

        figure(1); clf;
        [FI,P,S,X,E] = fit1DGabor( 'V2', dat, FI, P, S, X, E, n, length(experiments) );
    end

    if correct_screen_disparity
        save( [saveDir,'resultsV2.mat'], 'FI', 'P', 'S', 'X', 'E', 'experiments' );
    else
        save( [saveDir,'resultsV2_no correction.mat'], 'FI', 'P', 'S', 'X', 'E', 'experiments' );
    end
end

% MT
if any(contains(areas,'MT'))
    % load and preprocess data
    experiments = loadDataMT(correct_screen_disparity);

    % fit models
    FI    = []; % Fisher information (for fit1DGabor)
    P     = []; % parameters
    S     = []; % interpolated spike rates
    X     = []; % horizontal disparity
    E     = []; % fitting error

    for n = 1 : length(experiments)
        dat = experiments{n}.dat;

        figure(1); clf;
        [FI,P,S,X,E] = fit1DGabor( 'MT', dat, FI, P, S, X, E, n, length(experiments) );
    end

    if correct_screen_disparity
        save( [saveDir,'resultsMT.mat'], 'FI', 'P', 'S', 'X', 'E', 'experiments' );
    else
        save( [saveDir,'resultsMT_no correction.mat'], 'FI', 'P', 'S', 'X', 'E', 'experiments' );
    end
end

end