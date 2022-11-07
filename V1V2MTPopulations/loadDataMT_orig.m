function[experiments] = loadDataMT()

D = dir( 'dataMT/DispTunRaw*.mat' );
experiments = {}; % initialize

cnt = 1;

for k = 1 : numel(D)
    
    load( ['dataMT/' D(k).name] );
    
    % set valid indices for disparity trials
    dispInds = ddat.disparity~=-9999 & ddat.disparity~=-99 & ddat.disparity~=99 & ddat.disparity~=98;
    
    % horizontal disparity levels and counts
    dx = ddat.disparity(dispInds);
    count = ddat.firing_rate(dispInds);
    
    % calculate mean spike rate for each disparity
    [udx,~,ic] = unique(dx); % get unique disparities tested and indices for these disparities (indices used below for anova)
    mean_count = zeros(size(udx));
    trials = [];
    for d = 1:length(udx)
        trials(d) = sum(dx==udx(d));
        mean_count(d) = mean(count(dx==udx(d))); % count number of trials and mean spike count
    end
    mean_repeats = mean(trials);
    
    % to match V1/V2 data, we run an ANOVA and only keep cells with p <
    % 0.01 (threshold Bruce used to get samples for us
    %for ii = 1:size(neuronFitPars,1)
    groups = ic;
    p  = anova1(count,groups,'off');
    
    % also only keep neurons with average repeats 3 or more
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
    
    experiments{cnt}.dat            = [udx mean_count];
    experiments{cnt}.fn             = D(k).name;
    experiments{cnt}.mean_repeats   = mean_repeats;
    cnt = cnt + 1;
    
end



