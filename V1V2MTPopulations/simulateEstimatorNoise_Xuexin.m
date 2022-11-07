clear all; close all;

% determine the noise of an estimator for disparity based on these
% populations

% Wei & Stocker point out that FI/MI relationship holds best if the
% estimator noise is Gaussian
% https://www.sas.upenn.edu/~astocker/lab/publications-files/journals/NC2016/Wei_Stocker2016.pdf
% https://hal.archives-ouvertes.fr/hal-00143781/document

% disparities to test estimator 
disps = -2:.2:2;

% a test population of identical shifted Gaussians with mean at each disp level and sigma = 0.5
Gauss.P = [disps' repmat(0.5,numel(disps),1)];

% shifted sine waves (freq and phase)
Sine.P = [repmat(4,numel(disps)*3,1) linspace(-2*pi,2*pi,numel(disps)*3)'];

% empirical populations
V1 = load('resultsV1.mat');
V2 = load('resultsV2.mat');
MT = load('resultsMT.mat');

% compute preferred disparity for empirical tuning curves
disps_fine = -4:.0001:4;

for v = 1:size(V1.P,1)
    resps_tmp = V1.P(v,1) + V1.P(v,2)*exp( -(disps_fine-V1.P(v,3)).^2 / (2*V1.P(v,4)^2) ) .* cos(2*pi*V1.P(v,5)*(disps_fine-V1.P(v,3))+V1.P(v,6));
    pref_disp = disps_fine(resps_tmp == max(resps_tmp));
    if numel(pref_disp) > 1; display('V1 cell with 2 peaks'); end
    pref_disp = pref_disp(abs(pref_disp) == min(abs(pref_disp))); % when multiple, take the one with the smallest disparity
    v1_pref(v) = pref_disp;
end


for v = 1:size(V2.P,1)
    resps_tmp = V2.P(v,1) + V2.P(v,2)*exp( -(disps_fine-V2.P(v,3)).^2 / (2*V2.P(v,4)^2) ) .* cos(2*pi*V2.P(v,5)*(disps_fine-V2.P(v,3))+V2.P(v,6));
    pref_disp = disps_fine(resps_tmp == max(resps_tmp));
    if numel(pref_disp) > 1; display('V2 cell with 2 peaks'); end
    pref_disp = pref_disp(abs(pref_disp) == min(abs(pref_disp))); % when multiple, take the one with the smallest disparity
    v2_pref(v) = pref_disp;
end

for v = 1:size(MT.P,1)
    resps_tmp = MT.P(v,1) + MT.P(v,2)*exp( -(disps_fine-MT.P(v,3)).^2 / (2*MT.P(v,4)^2) ) .* cos(2*pi*MT.P(v,5)*(disps_fine-MT.P(v,3))+MT.P(v,6));
    pref_disp = disps_fine(resps_tmp == max(resps_tmp));
    if numel(pref_disp) > 1; display('MT cell with 2 peaks'); end
    pref_disp = pref_disp(abs(pref_disp) == min(abs(pref_disp))); % when multiple, take the one with the smallest disparity
    mt_pref(v) = pref_disp;
end


% record mean population response for all possible stimuli

% for each possible disparity
for s = 1:length(disps)
    
    % stimulus value
    stim_disp = disps(s);
    
    % measure response of each neuron
    
    % identical Gaussians
    for g = 1:size(Gauss.P,1)
        g_resps(s,g) = exp( (-(stim_disp - Gauss.P(g,1))^2) / (2*Gauss.P(g,2)^2) );
    end
    
    % sine waves
    for n = 1:size(Sine.P,1)
        s_resps(s,n) = sin(stim_disp*Sine.P(n,1) + Sine.P(n,2));
    end
    
    % V1 Gabors
    for v = 1:size(V1.P,1)
        v1_resps(s,v) = V1.P(v,1) + V1.P(v,2)*exp( -(stim_disp-V1.P(v,3)).^2 / (2*V1.P(v,4)^2) ) .* cos(2*pi*V1.P(v,5)*(stim_disp-V1.P(v,3))+V1.P(v,6));
    end
    
    % V2 Gabors
    for v = 1:size(V2.P,1)
        v2_resps(s,v) = V2.P(v,1) + V2.P(v,2)*exp( -(stim_disp-V2.P(v,3)).^2 / (2*V2.P(v,4)^2) ) .* cos(2*pi*V2.P(v,5)*(stim_disp-V2.P(v,3))+V2.P(v,6));
    end
    
    % MT Gabors
    for v = 1:size(MT.P,1)
        mt_resps(s,v) = MT.P(v,1) + MT.P(v,2)*exp( -(stim_disp-MT.P(v,3)).^2 / (2*MT.P(v,4)^2) ) .* cos(2*pi*MT.P(v,5)*(stim_disp-MT.P(v,3))+MT.P(v,6));
    end
end

% normalize responses
v1_resps_norm = (v1_resps - min(v1_resps))./range(v1_resps);
v2_resps_norm = (v2_resps - min(v2_resps))./range(v2_resps);
mt_resps_norm = (mt_resps - min(mt_resps))./range(mt_resps);

if(0)
    % for each unit, plot the response to each stimulus disparity
    fig_g       = figure; hold on; sgtitle('Gaussians - all responses');
    fig_sine    = figure; hold on; sgtitle('Sine Waves - all responses');
    fig_v1      = figure; hold on; sgtitle('V1 - all responses');
    fig_v2      = figure; hold on; sgtitle('V2 - all responses');
    fig_mt      = figure; hold on; sgtitle('MT - all responses');
    
    % for each possible disparity
    for s = 1:length(disps)
        
        % GAUSSIANS
        this_resp = g_resps(s,:);
        
        figure(fig_g); hold on;
        subplot(5,5,s); hold on;
        plot(this_resp,'ko');
        xlabel('neuron index (deg)'); ylabel('response');
        
        % SINE WAVES
        this_resp = s_resps(s,:);
        
        figure(fig_sine); hold on;
        subplot(5,5,s); hold on;
        plot(this_resp,'ko');
        xlabel('neuron index (deg)'); ylabel('response');
        
        % V1
        this_resp = v1_resps(s,:);
        
        figure(fig_v1); hold on;
        subplot(5,5,s); hold on;
        plot(this_resp,'ko');
        xlabel('neuron index (deg)'); ylabel('response');
        
        % V2
        this_resp = v2_resps(s,:);
        
        figure(fig_v2); hold on;
        subplot(5,5,s); hold on;
        plot(this_resp,'ko');
        xlabel('neuron index (deg)'); ylabel('response');
        
        % MT
        this_resp = mt_resps(s,:);
        
        figure(fig_mt); hold on;
        subplot(5,5,s); hold on;
        plot(this_resp,'ko');
        xlabel('neuron index (deg)'); ylabel('response');
        
    end
end


% now, run estimation based on population response

% for everything but the sine waves, can just estimate using the preferred
% disparity, but this is messy since some neurons have multiple peaks
fig_g_pop = figure; hold on; sgtitle('Gaussians - labeled line code - responses to given disparity');
fig_v1_pop = figure; hold on; sgtitle('V1 - labeled line code - responses to given disparity');
fig_v2_pop = figure; hold on; sgtitle('V2 - labeled line code - responses to given disparity');
fig_mt_pop = figure; hold on; sgtitle('MT - labeled line code - responses to given disparity');

fig_g_sampling = figure; hold on; sgtitle('Gaussians - labeled line code - decoder sampling distribution');
fig_v1_sampling = figure; hold on; sgtitle('V1 - labeled line code - decoder sampling distribution');
fig_v2_sampling = figure; hold on; sgtitle('V2 - labeled line code - decoder sampling distribution');
fig_mt_sampling = figure; hold on; sgtitle('MT - labeled line code - decoder sampling distribution');

ntrials = 1000; % number of trials to simulate
T       = 1; % time constant

% for each possible disparity
for s = 1:length(disps)
    
    % stimulus value
    stim_disp = disps(s);
    
    % Gaussians
    
    % average response of each unit, given this disparity
    
    % compute D0
    % D0 = entropy - 0.5*ln(2*pi*e / FI);
    input_pdf               = g_resps(s,:)/sum(g_resps(s,:)); % convert to pdf
    bin_val                 = (input_pdf .* log(input_pdf./(Gauss.P(2,1) - Gauss.P(1,1)))); % Compute value at each bin
    bin_val(input_pdf == 0) = 0; % Make sure values = 0 for events with 0 probability
    entropy_num             = -sum(bin_val); % Compute H(X) entropy
    FI                      = sum( ( (diff(log(input_pdf)) ./ (Gauss.P(2,1) - Gauss.P(1,1)) ).^2).*input_pdf(1:end-1) );
    D0                      = entropy_num - 0.5*log( (2*pi*exp(1)) / FI );
    
    figure(fig_g_pop); hold on;
    subplot(5,5,s); hold on; title(['D0 = ' num2str(D0,2)]);
    plot(Gauss.P(:,1),g_resps(s,:),'ko');
    plot([stim_disp stim_disp],[0 max(g_resps(s,:))],'r-');
    xlabel('preferred disparity (deg)'); ylabel('response');
    
    % compute probability of disparity, given the observed responses
    for t = 1:ntrials
        these_mean_resps    = T*g_resps(s,:); % mean resps for this disparity based on tuning curve
        these_noisy_resps   = poissrnd(these_mean_resps); % noisy responses on one trial with Poisson noise
        ind                 = find(these_noisy_resps == max(these_noisy_resps)); % find neuron with max response
        if numel(ind) > 1; ind = randsample(ind,1); end % if multiple, pick one at random
        g_disparity_estimate(t) = disps(ind);
    end
    
    figure(fig_g_sampling); hold on;
    subplot(5,5,s); hold on;
    histogram(g_disparity_estimate,disps)
    plot(stim_disp,0,'ro');
    xlabel('decoded disparity (deg)'); ylabel('probability');
    
    
    % V1
    
    % average response of each unit, given this disparity
    
    % sort array, average cells with same pref disparity, and interpolate to a regular grid
    [B,I]           = sort(v1_pref);
    v1_pref_sort    = B;
    v1_resps_smooth = smooth(v1_pref,v1_resps_norm(s,:),0.1,'rloess');
    v1_resps_sort   = v1_resps_smooth(I);
    [C,ia,idx]      = unique(v1_pref_sort,'stable');
    v1_pref_sort    = accumarray(idx,v1_pref_sort,[],@mean); 
    v1_resps_sort   = accumarray(idx,v1_resps_sort,[],@mean); 
    v1_resps_interp = interp1(v1_pref_sort,v1_resps_sort,Gauss.P(:,1),'linear',eps); % interpolate population response to a linear lattice
    
    input_pdf       = v1_resps_interp/sum(v1_resps_interp); % convert to pdf
    bin_val         = (input_pdf .* log(input_pdf./(Gauss.P(2,1) - Gauss.P(1,1)))); % Compute value at each bin
    bin_val(input_pdf == 0) = 0; % Make sure values = 0 for events with 0 probability
    entropy_num     = -sum(bin_val); % Compute H(X) entropy
    FI              = sum( ( (diff(log(input_pdf)) ./ (Gauss.P(2,1) - Gauss.P(1,1)) ).^2).*input_pdf(1:end-1) );
    D0              = entropy_num - 0.5*log( (2*pi*exp(1)) / FI );
    
    figure(fig_v1_pop); hold on;
    subplot(5,5,s); hold on; title(['D0 = ' num2str(D0,2)]);
    plot(v1_pref,v1_resps_norm(s,:),'ko');
    plot(v1_pref,smooth(v1_pref,v1_resps_norm(s,:),0.1,'rloess'),'r.');
    plot(Gauss.P(:,1),v1_resps_interp,'co','markerfacecolor','c');
    plot([stim_disp stim_disp],[0 max(v1_resps_norm(s,:))],'r-');
    xlabel('preferred disparity (deg)'); ylabel('response');
    
    % compute probability of disparity, given the observed responses
    % NOTE: using normed responses scaled by T because otherwise high spiking neurons dominate
    for t = 1:ntrials
        these_mean_resps        = T*v1_resps_norm(s,:); % mean resps for this disparity based on tuning curve
        these_noisy_resps       = poissrnd(these_mean_resps); % noisy responses on one trial with Poisson noise
        ind                     = find(these_noisy_resps == max(these_noisy_resps)); % find neuron with max response
        if numel(ind) > 1; ind = randsample(ind,1); end % if multiple, pick one at random
        v1_disparity_estimate(t) = v1_pref(ind);
    end
    
    figure(fig_v1_sampling); hold on;
    subplot(5,5,s); hold on;
    histogram(v1_disparity_estimate,disps)
    plot(stim_disp,0,'ro');
    xlabel('decoded disparity (deg)'); ylabel('probability');
    
    % V2
    
    % sort array, average cells with same pref disparity, and interpolate
    % to a regular grid
    [B,I] = sort(v2_pref);
    v2_pref_sort = B;
    v2_resps_smooth = smooth(v2_pref,v2_resps_norm(s,:),0.1,'rloess');
    v2_resps_sort = v2_resps_smooth(I);
    [C,ia,idx] = unique(v2_pref_sort,'stable');
    v2_pref_sort = accumarray(idx,v2_pref_sort,[],@mean); 
    v2_resps_sort = accumarray(idx,v2_resps_sort,[],@mean); 
    v2_resps_interp = interp1(v2_pref_sort,v2_resps_sort,Gauss.P(:,1),'linear',eps); % interpolate population response to a linear lattice
    
    input_pdf = v2_resps_interp/sum(v2_resps_interp); % convert to pdf
    bin_val = (input_pdf .* log(input_pdf./(Gauss.P(2,1) - Gauss.P(1,1)))); % Compute value at each bin
    bin_val(input_pdf == 0) = 0; % Make sure values = 0 for events with 0 probability
    entropy_num = -sum(bin_val); % Compute H(X) entropy
    FI = sum( ( (diff(log(input_pdf)) ./ (Gauss.P(2,1) - Gauss.P(1,1)) ).^2).*input_pdf(1:end-1) );
    D0 = entropy_num - 0.5*log( (2*pi*exp(1)) / FI );
    
    figure(fig_v2_pop); hold on;
    subplot(5,5,s); hold on; title(['D0 = ' num2str(D0,2)]);
    plot(v2_pref,v2_resps_norm(s,:),'ko');
    plot(v2_pref,smooth(v2_pref,v2_resps_norm(s,:),0.1,'rloess'),'r.');
    plot(Gauss.P(:,1),v2_resps_interp,'co','markerfacecolor','c');
    plot([stim_disp stim_disp],[0 max(v2_resps_norm(s,:))],'r-');
    xlabel('preferred disparity (deg)'); ylabel('response');
    
    % compute probability of disparity, given the observed responses
    % NOTE: using normed responses scaled by T because otherwise high spiking neurons dominate
    for t = 1:ntrials
        these_mean_resps        = T*v2_resps_norm(s,:); % mean resps for this disparity based on tuning curve
        these_noisy_resps       = poissrnd(these_mean_resps); % noisy responses on one trial with Poisson noise
        ind                     = find(these_noisy_resps == max(these_noisy_resps)); % find neuron with max response
        if numel(ind) > 1; ind = randsample(ind,1); end % if multiple, pick one at random
        v2_disparity_estimate(t) = v2_pref(ind);
    end
    
    figure(fig_v2_sampling); hold on;
    subplot(5,5,s); hold on;
    histogram(v2_disparity_estimate,disps)
    plot(stim_disp,0,'ro');
    xlabel('decoded disparity (deg)'); ylabel('probability');
    
    % MT
    
    % sort array, average cells with same pref disparity, and interpolate
    % to a regular grid
    [B,I] = sort(mt_pref);
    mt_pref_sort = B;
    mt_resps_smooth = smooth(mt_pref,mt_resps_norm(s,:),0.1,'rloess');
    mt_resps_sort = mt_resps_smooth(I);
    [C,ia,idx] = unique(mt_pref_sort,'stable');
    mt_pref_sort = accumarray(idx,mt_pref_sort,[],@mean); 
    mt_resps_sort = accumarray(idx,mt_resps_sort,[],@mean); 
    mt_resps_interp = interp1(mt_pref_sort,mt_resps_sort,Gauss.P(:,1),'linear',eps); % interpolate population response to a linear lattice
    
    input_pdf = mt_resps_interp/sum(mt_resps_interp); % convert to pdf
    bin_val = (input_pdf .* log(input_pdf./(Gauss.P(2,1) - Gauss.P(1,1)))); % Compute value at each bin
    bin_val(input_pdf == 0) = 0; % Make sure values = 0 for events with 0 probability
    entropy_num = -sum(bin_val); % Compute H(X) entropy
    FI = sum( ( (diff(log(input_pdf)) ./ (Gauss.P(2,1) - Gauss.P(1,1)) ).^2).*input_pdf(1:end-1) );
    D0 = entropy_num - 0.5*log( (2*pi*exp(1)) / FI );
    
    figure(fig_mt_pop); hold on;
    subplot(5,5,s); hold on; title(['D0 = ' num2str(D0,2)]);
    plot(mt_pref,mt_resps_norm(s,:),'ko');
    plot(mt_pref,smooth(mt_pref,mt_resps_norm(s,:),0.1,'rloess'),'r.');
    plot(Gauss.P(:,1),mt_resps_interp,'co','markerfacecolor','c');
    plot([stim_disp stim_disp],[0 max(mt_resps_norm(s,:))],'r-');
    xlabel('preferred disparity (deg)'); ylabel('response');
    
    % compute probability of disparity, given the observed responses
    % NOTE: using normed responses scaled by T because otherwise high spiking neurons dominate
    for t = 1:ntrials
        these_mean_resps        = T*mt_resps_norm(s,:); % mean resps for this disparity based on tuning curve
        these_noisy_resps       = poissrnd(these_mean_resps); % noisy responses on one trial with Poisson noise
        ind                     = find(these_noisy_resps == max(these_noisy_resps)); % find neuron with max response
        if numel(ind) > 1; ind = randsample(ind,1); end % if multiple, pick one at random
        mt_disparity_estimate(t) = mt_pref(ind);
    end
    
    figure(fig_mt_sampling); hold on;
    subplot(5,5,s); hold on;
    histogram(mt_disparity_estimate,disps)
    plot(stim_disp,0,'ro');
    xlabel('decoded disparity (deg)'); ylabel('probability');
    
end

keyboard

% now lets try Euclidian distance
fig_g_euc = figure; hold on; sgtitle('Gaussians - Euclidian distance');
fig_sine_euc = figure; hold on; sgtitle('Sine Waves - Euclidian distance');
fig_v1_euc = figure; hold on; sgtitle('V1 - Euclidian distance');
fig_v2_euc = figure; hold on; sgtitle('V2 - Euclidian distance');
fig_mt_euc = figure; hold on; sgtitle('MT - Euclidian distance');

% for each possible disparity
for s = 1:length(disps)
    
    % stimulus value
    stim_disp = disps(s);
    
    % GAUSSIANS
    % response to this stimulus value
    this_resp = g_resps(s,:);
    
    % euclidian distance between response to this disparity and all other
    % disparities
    resps_diff = sqrt(sum((g_resps - repmat(this_resp,size(g_resps,1),1)).^2,2));
    
    figure(fig_g_euc); hold on;
    subplot(5,5,s); hold on;
    plot(disps,resps_diff,'ko');
    plot([stim_disp stim_disp],[0 max(resps_diff(~isinf(resps_diff)))],'r-');
    xlabel('disparity (deg)'); ylabel('euclidian distance');
    
    
    % SINE WAVES
    % response to this stimulus value
    this_resp = s_resps(s,:);
    
    % euclidian distance between response to this disparity and all other
    % disparities
    resps_diff = sqrt(sum((s_resps - repmat(this_resp,size(s_resps,1),1)).^2,2));
    
    figure(fig_sine_euc); hold on;
    subplot(5,5,s); hold on;
    plot(disps,resps_diff,'ko');
    plot([stim_disp stim_disp],[0 max(resps_diff(~isinf(resps_diff)))],'r-');
    xlabel('disparity (deg)'); ylabel('euclidian distance');
    
    
    % V1
    
    % response to this stimulus value
    this_resp = v1_resps(s,:);
    
    % euclidian distance between response to this disparity and all other
    % disparities
    resps_diff = sqrt(sum((v1_resps - repmat(this_resp,size(v1_resps,1),1)).^2,2));
    
    figure(fig_v1_euc); hold on;
    subplot(5,5,s); hold on;
    plot(disps,resps_diff,'ko');
    plot([stim_disp stim_disp],[0 max(resps_diff(~isinf(resps_diff)))],'r-');
    xlabel('disparity (deg)'); ylabel('euclidian distance');
    
    
    % V2
    
    % response to this stimulus value
    this_resp = v2_resps(s,:);
    
    % euclidian distance between response to this disparity and all other
    % disparities
    resps_diff = sqrt(sum((v2_resps - repmat(this_resp,size(v2_resps,1),1)).^2,2));
    
    figure(fig_v2_euc); hold on;
    subplot(5,5,s); hold on;
    plot(disps,resps_diff,'ko');
    plot([stim_disp stim_disp],[0 max(resps_diff(~isinf(resps_diff)))],'r-');
    xlabel('disparity (deg)'); ylabel('euclidian distance');
    
    
    % MT
    
    % response to this stimulus value
    this_resp = mt_resps(s,:);
    
    % euclidian distance between response to this disparity and all other
    % disparities
    resps_diff = sqrt(sum((mt_resps - repmat(this_resp,size(mt_resps,1),1)).^2,2));
    
    figure(fig_mt_euc); hold on;
    subplot(5,5,s); hold on;
    plot(disps,resps_diff,'ko');
    plot([stim_disp stim_disp],[0 max(resps_diff(~isinf(resps_diff)))],'r-');
    xlabel('disparity (deg)'); ylabel('euclidian distance');
end

if(0)
    
    % check D0 calculation using an example Gaussian
    x = linspace(-5,5,1000);
    
    % define a Gaussian
    mu = 0;
    sigma = 1;
    gauss = (1/(sigma*sqrt(2*pi)))*exp(-0.5*(x-mu/sigma).^2);
    
    % add a bit of noise
    gauss = gauss + 0.1*abs(randn(1,numel(gauss)));
    
    % compute entropy
    input_pdf   = gauss/sum(gauss); % convert to pdf
    bin_val     = (input_pdf .* log(input_pdf./(x(2) - x(1)))); % Compute value at each bin
    bin_val(input_pdf == 0) = 0; % Make sure values = 0 for events with 0 probability
    entropy_num = -sum(bin_val) % Compute H(X) entropy
    
    entropy_exact = 0.5*log(2*pi*sigma^2) + 0.5
    
    % compute FI
    FI = sum( ( (diff(log(input_pdf)) ./(x(2) - x(1)) ).^2) .*input_pdf(1:end-1) ) % .*(x(2) - x(1)) )
    
    FI_D0 = 0.5*log( (2*pi*exp(1)) / FI );
    
    figure; hold on; title(['D0 = ' num2str(FI_D0,2)]);
    plot(gauss);
end


