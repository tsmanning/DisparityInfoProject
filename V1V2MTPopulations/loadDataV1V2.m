function [experiments] = loadDataV1V2(area,correct_screen_disparity)

% Define path to saved distribution data
splPath  = regexp(which('loadDataV1V2'),filesep,'split');
rootDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep];
basePath = [rootDir,'dataV1V2Combo',filesep];

addpath(genpath([rootDir,'helper_functions']));
addpath(genpath([rootDir,'screen_disparity_correction']));

% Load in raw data
load([basePath,'AllDTData.mat']);

experiments = {};
cnt = 1;
for x = 1:length(DTdata)

    if DTdata(x).visualarea == area

        % load data
        disparities = DTdata(x).dx;
        mean_count  = DTdata(x).resp;
        var_count   = cellfun(@(z) var(z),DTdata(x).counts);

        % some data has more than 2 dimensions
        var_count = squeeze(var_count);

        % count number of repeats per disparity and mean repeats across all disparities
        counts = DTdata(x).counts;
        repeats = [];
        for c = 1:length(counts)
            repeats(c) = length(DTdata(x).counts{c});
        end
        mean_repeats = mean(repeats);

        % grab RF location
        % atan(x/d) where x is displacement on screen, d is distance to nodal point from the screen
        x_pos = DTdata(x).rf(1);
        y_pos = DTdata(x).rf(2);


        % Note that Bruce already ran an ANOVA and only kept cells with p < 0.01

        % remove non disparity conditions
        remove      = disparities < -500;
        disparities = disparities(~remove);
        mean_count   = mean_count(~remove)';
        var_count   = var_count(~remove)';
        repeats = repeats(~remove);

        % ignore cases where there are repeated mean responses given for the same disparity
        if numel(unique(disparities)) < numel(disparities)
            display(['file abnormality - skipping ' num2str(x)]);
            continue
        end

        % apply screen disparity distortion correction to get disparity in
        % units where 0 = horopter rather than 0 = screen (also, corrects
        % for foreshortening on flat screen)
        if correct_screen_disparity
            horiz_ecc_deg = x_pos;
            disparities = screen2retDisp(disparities,horiz_ecc_deg,['V' num2str(area)]);
            disparities = disparities';
        end

        % ignore cases where orientation of disparity is not horizontal
        if ~isnan(DTdata(x).or) && abs(DTdata(x).or) ~= 90
            display(['disparity not horizontal - skipping ' num2str(x)]);
            continue;
        end

        % when RF orientation is -90 the disparity is still horizontal, but the direction is reversed.  So now negative values are uncrossed.
        if DTdata(x).or == -90
            display(['disparity axis reversed, fixing valuees ' num2str(x)]);
            disparities = -disparities;
        end


        % also only keep neurons with average repeats 3 or more
        if mean_repeats < 3
            display(['fewer than 3 repeats on average - skipping ' num2str(x)]);
            continue;
        end


        %         % omit cells with a max spike rate < 0.5 sps
        %         if max(mean_count) < 0.5
        %             display(['max spike rate less than 0.5 sps - skipping ' num2str(x)]);
        %             continue;
        %         end

        % omit cells fit to a very small range of disparities,
        % because this is hard to compare to MT data
        if range(disparities) < 1
            display(['disparity range too small - skipping ' num2str(x)]);
            continue;
        end


        % convert responses to spikes per second
        % .resp is in units of mean spike count, not divided by duration. Duration was 0.5sec  0.47 or 0.42 sec (different studies).  ]
        % There is a field in the dxresp that records a  number related to the stimulus duration.
        % (This gets complicated, on some trials the graphics card fails to keep up, so it takes 43 video frames to paint 42 actual frames.
        % The duration number in these files depends on the distribution of real durations in the expt.)
        % But I always just count spikes over either 420ms o470, or  500ms, whichever is the nominal (and so shortest) duration.
        % SO if "duration" is <  0.47, its 420ms.  if its > 0.47 and < 0.49, its 470ms (names all begin 'jbeG')  and if its > 0.49, its 500ms (names all begin 'jbeM' or 'lemM')
        if DTdata(x).duration < 0.47
            duration = 0.42;
        elseif DTdata(x).duration >= 0.47 && DTdata(x).duration < 0.49
            duration = .47;
        elseif DTdata(x).duration >= 0.49
            duration = .5;
        else
            error('invalid duration');
        end

        mean_count = mean_count/duration;
        var_count = var_count/duration;

        % store disparities and responses
        experiments{cnt}.dat = [disparities mean_count' var_count' repeats'];
        % mean repeats
        experiments{cnt}.mean_repeats = mean_repeats;

        % grab RF location
        % atan(x/d) where x is displacement on screen, d is distance to nodal point from the screen
        experiments{cnt}.x_pos = x_pos;
        experiments{cnt}.y_pos = y_pos;


        % store all other data about this cell
        experiments{cnt}.alldata = DTdata(x);

        cnt = cnt + 1;

    end

end


