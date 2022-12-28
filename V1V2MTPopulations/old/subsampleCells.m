function [r2Inds,eccInds] = subsampleCells(varargin)

% Take full sample of V1/V2/MT cells and resample them so we're only
% looking at those well-described by Gabor tuning, those with RF centers
% within 10deg of fixation, and roughly all in the same visual field
% regions

% Define path to saved tuning curve data
splPath = regexp(which('subsampleCells'),filesep,'split');
topDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
dataDir = [topDir,filesep,'analysisFiles',filesep];

saveOn  = 1;


%% Load fitting results
if nargin == 0
    load([dataDir,'resultsALL_final.mat']);
else
    V1 = varargin{1};
    V2 = varargin{2};
    MT = varargin{3};
end


%% Only keep neurons with R2 > 0.75
r2Inds.V1 = V1.r2 > 0.75;
r2Inds.V2 = V2.r2 > 0.75;
r2Inds.MT = MT.r2 > 0.75;

display(['R2 .75 V1 =  ' num2str(sum(r2Inds.V1))]);
display(['R2 .75 V2 =  ' num2str(sum(r2Inds.V2))]);
display(['R2 .75 MT =  ' num2str(sum(r2Inds.MT))]);


%% Crop everything to 10 deg and resample MT to same elevation/abs azimuth as V1/V2

% Grab VF subset of MT iwthin V1/V2 bounds
% get indices of MT RF's that overlap with V1/V2; allow points to be either
% left or right of fixation, at same eccentricities
% in practice, this is just the V1 range because all V2 is within V2
V1V2_x_range = max(abs([V1.x_pos V2.x_pos])); % don't care about left/right, just the max x-coordinate
V1V2_y_range = [min([V1.y_pos V2.y_pos]) max([V1.y_pos V2.y_pos])];

overlapping_RF_inds = abs(MT.x_pos) <= V1V2_x_range ...
                         & MT.y_pos >= V1V2_y_range(1) ...
                         & MT.y_pos <= V1V2_y_range(2);

% Cull cells outside of 10deg of eccentricities
eccInds.V1 = (sqrt(V1.x_pos.^2 + V1.y_pos.^2) <= 10) & r2Inds.V1;
eccInds.V2 = (sqrt(V2.x_pos.^2 + V2.y_pos.^2) <= 10) & r2Inds.V2;
eccInds.MT = (sqrt(MT.x_pos.^2 + MT.y_pos.^2) <= 10) & r2Inds.MT & overlapping_RF_inds;

% Report final number of included cells
display(['V1 <= 10deg:  ' num2str(sum(eccInds.V1))]);
display(['V2 <= 10deg:  ' num2str(sum(eccInds.V2))]);
display(['MT <= 10deg:  ' num2str(sum(eccInds.MT))]);

%% Save indices

if saveOn
save([dataDir,'r2Inds'],'r2Inds');
save([dataDir,'eccInds'],'eccInds');
end

end