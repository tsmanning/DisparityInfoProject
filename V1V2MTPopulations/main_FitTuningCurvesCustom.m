% Handle a couple of edge cases for tuning curve fitting

clear;

splPath  = regexp(which('main_FitTuningCurves.m'),filesep,'split');
topDir   = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
saveDir  = [topDir,'analysisFiles',filesep];

addpath([topDir,filesep,'helper_functions']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%area = 'V1';  % max fitted spike rate exceeds measured max by a factor of
%10, better fits acheived penalizing fits where the max fitted rate is more than
%2x the observed max

%n = 70;
%n = 74;
%n = 195;

%area = 'V2'; % max fitted spike rate exceeds measured max by a factor of
%10, better fits acheived penalizing fits where the max fitted rate is more than
%2x the observed max

%n = 485; 
%n = 513; 


%area = 'MT'; % gaps between disparity samples are too broad so we got fits
% with aliasing (high frequency that goes up and down through all points.
% Re-running with max freq as 3 cpd

%n = 18;
%n = 20;
%n = 24;
%n = 61;
%n = 81;
%n = 109;
%n = 139;
%n = 144;
%n = 145;
%n = 187;
%n = 197;
%n = 198;
%n = 204;
%n = 216;
%n = 220;
%n = 247;
%n = 271;
%n = 280;
%n = 293;
%n = 343;
%n = 344;
%n = 383;
%n = 442;

% Re-running with max freq as 1 cpd
%n = 285;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% flag indicating whether to correct the stimulus disparities to account
% for planar screen (horopter deviation and foreshortening)

correct_screen_disparity = 1;

% load and preprocess data
if strcmp(area,'V1')
    experiments = loadDataV1V2(1,correct_screen_disparity);
elseif strcmp(area,'V2')
    experiments = loadDataV1V2(2,correct_screen_disparity);
elseif strcmp(area,'MT')
    experiments = loadDataMT(correct_screen_disparity);
end


FI    = []; % Fisher information (for fit1DGabor)
P     = []; % parameters
S     = []; % interpolated spike rates
X     = []; % horizontal disparity
E     = []; % fitting error


dat = experiments{n}.dat; % pre-processed data

figure(1); clf;

[FI,P,S,X,E] = fit1DGabor( area, dat, FI, P, S, X, E, 1, 1 );

FIind = FI; Pind = P; Sind = S; Xind = X; Eind = E;

if strcmp(area,'V1')
    load resultsV1.mat;
elseif strcmp(area,'V2')
    load resultsV2.mat;
elseif strcmp(area,'MT')
    load resultsMT.mat;
end

FI(n,:)  = FIind;
P(n,:)   = Pind;
S(n).s   = Sind.s;
X(n).x   = Xind.x;
E(n)     = Eind;

if correct_screen_disparity
    suffix = '_no correction';
else
    suffix = '';
end

if strcmp(area,'V1')
    save( [saveDir,'resultsV1',suffix,'.mat'], 'FI', 'P', 'S', 'X', 'E', 'experiments' );
elseif strcmp(area,'V2')
    save( [saveDir,'resultsV2',suffix,'.mat'], 'FI', 'P', 'S', 'X', 'E', 'experiments' );
elseif strcmp(area,'MT')
    save( [saveDir,'resultsMT',suffix,'.mat'], 'FI', 'P', 'S', 'X', 'E', 'experiments' );
end
