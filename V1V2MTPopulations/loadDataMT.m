function [experiments] = loadDataMT(correct_screen_disparity)

% Define path to saved distribution data
splPath  = regexp(which('FanoFacCheck'),filesep,'split');
rootDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
analyDir = [rootDir,'analysisFiles',filesep];
baseDir  = [rootDir,'dataMT',filesep];

addpath(genpath([rootDir,'helper_functions']));
addpath(genpath([rootDir,'screen_disparity_correction']));

D = dir([baseDir,'DispTunRaw*.mat']);
experiments = {}; % initialize

% load meta-data to grab RF eccentricities
T_MT = readtable([baseDir,'metadataMT.xlsx']);

cnt = 1;

for k = 1 : numel(D)

    load([baseDir,D(k).name]);

    % set valid indices for disparity trials
    dispInds = ddat.disparity~=-9999 & ddat.disparity~=-99 & ddat.disparity~=99 & ddat.disparity~=98;

    % horizontal disparity levels and counts for all stimulus presentations
    dx = ddat.disparity(dispInds);
    count = ddat.firing_rate(dispInds);

    % calculate mean spike rate and variance for each unique disparity
    [disparities,~,ic] = unique(dx); % get unique disparities tested and indices for these disparities (indices used below for anova)
    mean_count = zeros(size(disparities));
    var_count = zeros(size(disparities));
    trials = zeros(size(disparities));

    for d = 1:length(disparities)
        trials(d) = sum(dx==disparities(d));
        mean_count(d) = mean(count(dx==disparities(d))); % count number of trials and mean spike count
        var_count(d) = var(count(dx==disparities(d))); % count number of trials and mean spike count
    end

    % store mean repeats
    mean_repeats = mean(trials);

    % compare filename to Greg's metadata in order to get RF eccentricity
    % get filename
    fname = D(k).name;
    fname = fliplr(strtok(fliplr(fname),'_'));
    fname = fname(1:end-4);

    % find row in Greg's metadata files
    this_row = find(strcmp(T_MT.FIleName,fname));

    % Grab RF location
    % positive x is right in VF, positive y is up in VF
    x_pos  = T_MT.RF_x(this_row); % note this does not account for vergence, see Jenny's notes
    y_pos  = T_MT.RF_y(this_row);


    % apply screen disparity distortion correction to get disparity in
    % units where 0 = horopter rather than 0 = screen (also, corrects
    % for foreshortening on flat screen)
    if correct_screen_disparity
        horiz_ecc_deg = x_pos;
        disparities = screen2retDisp(disparities,horiz_ecc_deg,['MT']);
    end


    % to match V1/V2 data, we run an ANOVA and only keep cells with p <
    % 0.01 (threshold Bruce used to get samples for us)
    groups = ic;
    p  = anova1(count,groups,'off');

    if p >= 0.01
        display(['p >= 0.01 - skipping ' num2str(k)]);
        continue
    end

    % also only keep neurons with average repeats 3 or more
    if mean_repeats < 3
        display(['fewer than 3 repeats on average - skipping ' num2str(k)]);
        continue
    end

    %     % omit cells with a max spike rate < 0.5 sps
    %     if max(mean_count) < 0.5
    %         display(['max spike rate less than 0.5 sps - skipping ' num2str(k)]);
    %         continue;
    %     end

    % store into data structure
    experiments{cnt}.dat            = [disparities' mean_count var_count trials];
    experiments{cnt}.fn             = D(k).name;
    experiments{cnt}.mean_repeats   = mean_repeats;

    % grab RF location
    % atan(x/d) where x is displacement on screen, d is distance to nodal point from the screen
    experiments{cnt}.x_pos = x_pos;
    experiments{cnt}.y_pos = y_pos;

    cnt = cnt + 1;

end



