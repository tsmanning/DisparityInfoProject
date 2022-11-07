clear all; close all;
addpath(genpath('./helper_functions'));
plotallfits = 0;

% load fitting results
V1 = load('resultsV1_final.mat');
V2 = load('resultsV2_final.mat');
MT = load('resultsMT_final.mat');



%% report how many cells in each population
display(['Num V1 units: ' num2str(numel(V1.E))]);
display(['Num V2 units: ' num2str(numel(V2.E))]);
display(['Num MT units: ' num2str(numel(MT.E))]);
display(['Num TOTAL units: ' num2str(numel(V1.E)+numel(V2.E)+numel(MT.E))]);

%% LOAD ADDITIONAL MT/V1 DATA AND APPEND TO DATA STRUCTURE

% load MT metadata
T_MT = readtable('metadataMT.xlsx');

% for each file name in experiments, find index from xls and get other data
for f = 1:length(MT.experiments)

    % get filename
    fname = MT.experiments{f}.fn;
    fname = fliplr(strtok(fliplr(fname),'_'));
    fname = fname(1:end-4);

    % find row in Greg's files
    this_row = find(strcmp(T_MT.FIleName,fname));

    % grab data from that row and assign it to this neuron
    % positive x is right in VF, positive y is up in VF
    MT.x_pos(f)             = T_MT.RF_x(this_row); % note this does not account for vergence, see Jenny's notes
    MT.y_pos(f)             = T_MT.RF_y(this_row);
    MT.GD_pref_disparity(f) = T_MT.PrefDisp(this_row);
    MT.GD_RF_diameter(f)    = T_MT.RFDiam_hand(this_row);
    MT.GD_R0(f)             = T_MT.R_0(this_row);
    MT.GD_A(f)              = T_MT.A(this_row);
    MT.GD_d0(f)             = T_MT.d_0(this_row);
    MT.GD_sigma(f)          = T_MT.sigma(this_row);
    MT.GD_f(f)              = T_MT.f(this_row);
    MT.GD_phi(f)            = T_MT.phi(this_row);

end

% grab V1 RF location data, much simpler
for v = 1:length(V1.experiments)

    % positive x is to the animals left in their visual field so we flip to match MT data, positive y is up
    V1.x_pos(v) = -V1.experiments{v}.x_pos;
    V1.y_pos(v) = V1.experiments{v}.y_pos;

end

% grab V2 RF location data, much simpler
for v = 1:length(V2.experiments)

    % positive x is to the animals left in their visual field so we flip to match MT data, positive y is up
    V2.x_pos(v) = -V2.experiments{v}.x_pos;
    V2.y_pos(v) = V2.experiments{v}.y_pos;

end

% calculate raw error R2 for each fit in each area
% also calculate some other tuning characteristics
areas   = {'V1','V2','MT'};
xg1     = -2 : 0.01 : 2; % fixed range to evaluate tuning function

for a = 1:length(areas)

    switch areas{a}

        case 'V1';  area = V1;
        case 'V2';  area = V2;
        case 'MT';  area = MT;
        otherwise;  error('invalid area')

    end

    for n = 1:length(area.experiments)

        % grab mean number of repeats per stimulus
        area.mean_repeats(n) = area.experiments{n}.mean_repeats;

        % mean response by disparity
        disps = area.experiments{n}.dat(:,1); % tested disparities
        resps = area.experiments{n}.dat(:,2); % measured responses;

        % calculate fitted spike rates for each tested disparity
        TF = area.P(n,1) + area.P(n,2)*exp( -(disps-area.P(n,3)).^2 / (2*area.P(n,4)^2) ) .* cos(2*pi*area.P(n,5)*(disps-area.P(n,3))+area.P(n,6));

        SSR = sum((resps - TF).^2);
        TSS = sum((resps - mean(resps)).^2);

        % R2
        area.r2(n)   = 1 - (SSR/TSS);

        % calculate fitted spike rates for full lattice of disparities
        g1 = area.P(n,1) + area.P(n,2)*exp( -(xg1-area.P(n,3)).^2 / (2*area.P(n,4)^2) ) .* cos(2*pi*area.P(n,5)*(xg1-area.P(n,3))+area.P(n,6));
        area.allresp(n,:) = g1;

        % max SPS
        area.maxsps(n) = max(g1);

        % preferred disparity
        pref_disp = xg1(g1 == max(g1));
        pref_disp = pref_disp(abs(pref_disp) == min(abs(pref_disp))); % when multiple, take the one with the smallest disparity
        area.pref_disp(n) = pref_disp;

        % take first derivative and find the max absolute slope
        deriv = diff(g1);
        area.max_slope(n) = max(abs(deriv));

        % use first derivative to determine whether the unit has even or
        % odd symmetry
        deriv_sign = sign(deriv); % is the derivative positive or negative?
        deriv_sign_change = diff(deriv_sign); % did the slope sign change direction?

        % if it changed from -1 to 1 an even number of times, it's
        % oddsymmetric
        area.oddsym(n) = mod(sum(abs(deriv_sign_change)==2),2) == 0;

    end

    switch areas{a}

        case 'V1';  V1 = area;
        case 'V2';  V2 = area;
        case 'MT';  MT = area;
        otherwise;  error('invalid area')

    end

end

% save these structures with extra data
save('resultsALL_final.mat','V1','V2','MT');

%% FITTING QUALITY

% histogram the mean repeats per stimulus disparity
figure; hold on;
edges = linspace(0,max([V1.mean_repeats V2.mean_repeats MT.mean_repeats]),20);
subplot(2,3,1); hold on; title('Mean repeats per stim for V1 fits');
h = histogram(V1.mean_repeats,edges);
h.FaceColor = ColorIt('b');
axis square; box on;
subplot(2,3,2); hold on; title('Mean repeats per stim for V2 fits');
h = histogram(V2.mean_repeats,edges);
h.FaceColor = ColorIt('g');
axis square; box on;
subplot(2,3,3); hold on; title('Mean repeats per stim for MT fits');
h = histogram(MT.mean_repeats,edges);
h.FaceColor = ColorIt('r');
axis square; box on;

% histogram the fitting errors -- note that these include the penality and
% are not normalized to mean rate
edges = linspace(0,max([V1.E' V2.E' MT.E']),20);
subplot(2,3,4); hold on; title('RMS error for V1 fits');
h = histogram(V1.E,edges);
h.FaceColor = ColorIt('b');
axis square; box on;
subplot(2,3,5); hold on; title('RMS error for V2 fits');
h = histogram(V2.E,edges);
h.FaceColor = ColorIt('g');
axis square; box on;
subplot(2,3,6); hold on; title('RMS error for MT fits');
h = histogram(MT.E,edges);
h.FaceColor = ColorIt('r');
axis square; box on;

saveas(gcf,'./plots/AssessFits/Fit_MeanRepeatsAndError.png');
saveas(gcf,'./plots/AssessFits/Fit_MeanRepeatsAndError.eps','epsc');


% histogram R2
figure; hold on;
edges = linspace(0,1,20);
subplot(1,3,1); hold on; title('V1 fits');
h = histogram(V1.r2,edges);
h.FaceColor = ColorIt('b');
h.FaceAlpha = 1;
ylim([0 250]);
xlabel('R2'); ylabel('frequency')
axis square; box on;
subplot(1,3,2); hold on; title('V2 fits');
h = histogram(V2.r2,edges);
h.FaceColor = ColorIt('g');
h.FaceAlpha = 1;
ylim([0 250]);
xlabel('R2'); ylabel('frequency')
axis square; box on;
subplot(1,3,3); hold on; title('MT fits');
h = histogram(MT.r2,edges);
h.FaceColor = ColorIt('r');
h.FaceAlpha = 1;
ylim([0 250]);
xlabel('R2'); ylabel('frequency')
axis square; box on;

% display median R2 for each area
display(['Median R2 V1 units: ' num2str(median(V1.r2))]);
display(['Median R2 V2 units: ' num2str(median(V2.r2))]);
display(['Median R2 MT units: ' num2str(median(MT.r2))]);

saveas(gcf,'./plots/AssessFits/Fit_R2.png');
saveas(gcf,'./plots/AssessFits/Fit_R2.eps','epsc');


%% COMPARE PARAMETERS ACROSS REGIONS
% these are the raw fitted parameters, which can be a bit confusing because
% they interact to determine the actual tuning curve shape (e.g., a
% wide/shallow envelope + high amplitude can produce the same max spike
% rate as a narrower envelope and a lower amplitude)
figure; hold on;
param_list = {'Offset(b)','Amplitude(a)','Envelope mean(d0)','Envelope std(sigma)','Frequency(f)','Phase(phi)'};
for p = 1:6

    subplot(2,3,p); hold on; title(param_list{p});
    distributionPlot([V1.P(:,p)],'xValues',3,'color',ColorIt('b'),'histOpt',1,'showMM',6,'xNames',{'V1'},'xyOri','flipped');
    distributionPlot([V2.P(:,p)],'xValues',2,'color',ColorIt('g'),'histOpt',1,'showMM',6,'xNames',{'V2'},'xyOri','flipped');
    distributionPlot([MT.P(:,p)],'xValues',1,'color',ColorIt('r'),'histOpt',1,'showMM',6,'xNames',{'MT'},'xyOri','flipped');
    if p == 1 || p == 2 || p == 4 || p == 5
        xlim([0 quantile([V1.P(:,p)' V2.P(:,p)' MT.P(:,p)'],.99)])
    end
    axis equal square; box on;

end

saveas(gcf,'./plots/AssessFits/Fit_Params.png');
saveas(gcf,'./plots/AssessFits/Fit_Params.eps','epsc');

%% COMMON SENSE PARAMETERS
% convert to a set of common sense parameters about the tuning curve
% baseline spike rate, max spike rate, preferred disparity, max abs derivative and circular phase
figure; hold on;

% baseline sps
subplot(2,3,1); hold on; title('Baseline spike rate (sps)');
%boxplot(V1.P(:,1),'Positions',-1,'Colors',ColorIt('b'),'Symbol','','Labels','V1')
%boxplot(V2.P(:,1),'Positions',0,'Colors',ColorIt('g'),'Symbol','','Labels','V2')
%boxplot(MT.P(:,1),'Positions',1,'Colors',ColorIt('r'),'Symbol','','Labels','MT')
distributionPlot([V1.P(:,1)],'xValues',3,'color',ColorIt('b'),'histOpt',1,'showMM',6,'xNames',{'V1'},'xyOri','flipped');
distributionPlot([V2.P(:,1)],'xValues',2,'color',ColorIt('g'),'histOpt',1,'showMM',6,'xNames',{'V2'},'xyOri','flipped');
distributionPlot([MT.P(:,1)],'xValues',1,'color',ColorIt('r'),'histOpt',1,'showMM',6,'xNames',{'MT'},'xyOri','flipped');
xlim([0 quantile([V1.P(:,1)' V2.P(:,1)' MT.P(:,1)'],.99)])
%xlim([-3 3])
axis square; box on;

% max sps
subplot(2,3,2); hold on; title('Maximum spike rate (sps)');
%boxplot(V1.maxsps,'Positions',-1,'Colors',ColorIt('b'),'Symbol','','Labels','V1')
%boxplot(V2.maxsps,'Positions',0,'Colors',ColorIt('g'),'Symbol','','Labels','V2')
%boxplot(MT.maxsps,'Positions',1,'Colors',ColorIt('r'),'Symbol','','Labels','MT')
distributionPlot(V1.maxsps','xValues',3,'color',ColorIt('b'),'histOpt',1,'showMM',6,'xNames',{'V1'},'xyOri','flipped');
distributionPlot(V2.maxsps','xValues',2,'color',ColorIt('g'),'histOpt',1,'showMM',6,'xNames',{'V2'},'xyOri','flipped');
distributionPlot(MT.maxsps','xValues',1,'color',ColorIt('r'),'histOpt',1,'showMM',6,'xNames',{'MT'},'xyOri','flipped');
xlim([0 quantile([V1.maxsps V2.maxsps MT.maxsps],.99)])
%xlim([-3 3])
axis square; box on;

% preferred disparity
subplot(2,3,3); hold on; title('Preferred disparity (deg)');
distributionPlot(V1.pref_disp','xValues',3,'color',ColorIt('b'),'histOpt',1,'showMM',6,'xNames',{'V1'},'xyOri','flipped');
distributionPlot(V2.pref_disp','xValues',2,'color',ColorIt('g'),'histOpt',1,'showMM',6,'xNames',{'V2'},'xyOri','flipped');
distributionPlot(MT.pref_disp','xValues',1,'color',ColorIt('r'),'histOpt',1,'showMM',6,'xNames',{'MT'},'xyOri','flipped');
xlim([quantile([V1.pref_disp V2.pref_disp MT.pref_disp],.01) quantile([V1.pref_disp V2.pref_disp MT.pref_disp],.99)])
axis square; box on;

% max absolute slope
subplot(2,3,4); hold on; title('Max slope (sps/deg)');
distributionPlot(V1.max_slope','xValues',3,'color',ColorIt('b'),'histOpt',1,'showMM',6,'xNames',{'V1'},'xyOri','flipped');
distributionPlot(V2.max_slope','xValues',2,'color',ColorIt('g'),'histOpt',1,'showMM',6,'xNames',{'V2'},'xyOri','flipped');
distributionPlot(MT.max_slope','xValues',1,'color',ColorIt('r'),'histOpt',1,'showMM',6,'xNames',{'MT'},'xyOri','flipped');
xlim([0 quantile([V1.max_slope V2.max_slope MT.max_slope],.95)])
axis square; box on;

% odd symmetry
subplot(2,3,5); hold on; title('Percentage odd symmetric');
barh(3,sum(V1.oddsym/numel(V1.oddsym)),'FaceColor',ColorIt('b'));
barh(2,sum(V2.oddsym/numel(V2.oddsym)),'FaceColor',ColorIt('g'));
barh(1,sum(MT.oddsym/numel(MT.oddsym)),'FaceColor',ColorIt('r'));
xlim([0 1])
axis square; box on;

saveas(gcf,'./plots/AssessFits/Fit_Characteristics.png');
saveas(gcf,'./plots/AssessFits/Fit_Characteristics.eps','epsc');

%% PLOT EXAMPLE FITS for figure

xg1 = [-2 : 0.01 : 2]; %fixed spatial range
for a = 1:length(areas)

    switch areas{a}

        case 'V1';  area = V1; n = 2; yrange = [30 65];
        case 'V2';  area = V2; n = 2; yrange = [30 80];
        case 'MT';  area = MT; n = 61; yrange = [20 80];
        otherwise;  error('invalid area')

    end

    g1 = area.P(n,1) + area.P(n,2)*exp( -(xg1-area.P(n,3)).^2 / (2*area.P(n,4)^2) ) .* cos(2*pi*area.P(n,5)*(xg1-area.P(n,3))+area.P(n,6));

    figure; hold on;
    scatter(area.experiments{n}.dat(:,1),area.experiments{n}.dat(:,2),'k','filled');
    plot(xg1,g1,'b-','linewidth',2);
    axis square; xticks([]); yticks([]); box on; ylim(yrange); xlim([-2 2]);
    saveas(gcf,['./plots/AssessFits/ExampleFits/ExampleFit_' areas{a} '_cell' num2str(n) '.png']);
    saveas(gcf,['./plots/AssessFits/ExampleFits/ExampleFit_' areas{a} '_cell' num2str(n) '.eps'],'epsc');

    figure; hold on;
    plot(xg1,area.FI(n,:),'k-','linewidth',2);
    axis square; xticks([]); yticks([]); box on;
    saveas(gcf,['./plots/AssessFits/ExampleFits/ExampleFit_' areas{a} '_cell' num2str(n) '_FI.png']);
    saveas(gcf,['./plots/AssessFits/ExampleFits/ExampleFit_' areas{a} '_cell' num2str(n) '_FI.eps'],'epsc');

end

%% COMPARE MT params TO GD FITS

% Compare the fits between our routine and the fits from DeAngelis & Uka
% Even though the fits match visually, there are some differences because
% - D&U fit to all data and we fit to mean
% - some parameter combos of Gabors result in very similar fits
figure; hold on;
subplot(2,3,1); hold on; title('offset');
scatter(MT.GD_R0,MT.P(:,1)','k','filled');
refline(1,0); xlabel('original fit'); ylabel('our fit');
axis equal square; box on;

subplot(2,3,2); hold on; title('amplitude');
scatter(MT.GD_A,MT.P(:,2)','k','filled');
refline(1,0); xlabel('original fit'); ylabel('our fit');
axis equal square; box on;

subplot(2,3,3); hold on; title('envelope mean');
scatter(MT.GD_d0,MT.P(:,3)','k','filled');
refline(1,0); xlabel('original fit'); ylabel('our fit');
axis equal square; box on;

subplot(2,3,4); hold on; title('envelope sigma');
scatter(MT.GD_sigma,MT.P(:,4)','k','filled');
refline(1,0); xlabel('original fit'); ylabel('our fit');
axis equal square; box on;

subplot(2,3,5); hold on; title('frequency');
scatter(MT.GD_f,MT.P(:,5)','k','filled');
refline(1,0); xlabel('original fit'); ylabel('our fit');
axis equal square; box on;

subplot(2,3,6); hold on; title('phase');
scatter(MT.GD_phi,MT.P(:,6)','k','filled');
refline(1,0); xlabel('original fit'); ylabel('our fit');
axis equal square; box on;

saveas(gcf,['./plots/AssessFits/Fit_MT_ComparisonToDandU.png']);

if(plotallfits)
    %% PLOT ALL FITS
    for a = 1:length(areas)

        switch areas{a}

            case 'V1';  area = V1;
            case 'V2';  area = V2;
            case 'MT';  area = MT;
            otherwise;  error('invalid area')

        end

        % counters
        pcnt = 1;
        fcnt = 1;
        for n = 1:length(area.experiments)

            

            if mod(n,66) == 1
                if fcnt > 1
                    saveas(gcf,['./plots/AssessFits/AllFits/Fit_' areas{a} '_' num2str(fcnt-1) '.png']);
                    close gcf;
                end
                figure; hold on;
                setupfig(18,10,10);
                pcnt = 1;
                fcnt = fcnt + 1;
            end

            %g1 = area.P(n,1) + area.P(n,2)*exp( -(xg1-area.P(n,3)).^2 / (2*area.P(n,4)^2) ) .* cos(2*pi*area.P(n,5)*(xg1-area.P(n,3))+area.P(n,6));

            subplot(6,11,pcnt); hold on; title([ num2str(n) ' R2 = ' num2str(area.r2(n),3)]);

            % data
            scatter(area.experiments{n}.dat(:,1),area.experiments{n}.dat(:,2),'k','filled');

            % our fit
            g1 = area.P(n,1) + area.P(n,2)*exp( -(xg1-area.P(n,3)).^2 / (2*area.P(n,4)^2) ) .* cos(2*pi*area.P(n,5)*(xg1-area.P(n,3))+area.P(n,6));

            plot(xg1,g1,'b-');
            box on;

            % if MT, also plot original fitted params in green
            if strcmp(areas{a},'MT')
                g1_GD = MT.GD_R0(n) + MT.GD_A(n)*exp( -(xg1-MT.GD_d0(n)).^2 / (2*MT.GD_sigma(n)^2) ) .* cos(2*pi*MT.GD_f(n)*(xg1-MT.GD_d0(n))+MT.GD_phi(n));
                plot(xg1,g1_GD,'g-');
            end

            % flag any bad fits
            if area.r2(n) < 0.2
                display([ areas{a} ' bad fit: ' num2str(n)]);
                plot(xg1,g1,'r-');
            end

            pcnt = pcnt + 1;

%             if abs(area.P(n,6) - pi/2) < .1
%                 keyboard
%             end


        end
        saveas(gcf,['./plots/AssessFits/AllFits/Fit_' areas{a} '_' num2str(fcnt-1) '.png']);

    end

end




