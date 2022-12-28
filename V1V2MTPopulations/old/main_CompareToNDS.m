% This script finds the best power law that relates the measured disparity
% distributions from the BORIS image dataset to the FI distributions from
% the neural datasets

clear all;
close all;

% Use original or resampled disparity histograms
% flag = '';
flag = '_resampledEccentricities';
%flag = '_resampledV1MT';

% Define path to saved distribution data
splPath = regexp(which('main_CompareToNDS'),filesep,'split');
imStatsDir = [filesep,fullfile(splPath{1:numel(splPath)-2}),filesep,'SceneStatsAnalysis/savedImageStats_BORISdataset/'];
fiDataDir  = [filesep,fullfile(splPath{1:numel(splPath)-1}),filesep,'analysisFiles',filesep];

% Toggle save figures on/off
saveOn = 0;


%% Collect data and define 
% Load in FI data
dat = load([fiDataDir,'results_GaussFit' flag '.mat']);

% Define resolution to which we're going to downsample FI (this resolution
% matches the scene statistics)
res = 52;
lb  = -2;
ub  = 2;edges_disp = linspace(lb,ub,res);
cntr_disp  = edges_disp(1:end-1) + diff(edges_disp(1:2))/2;

% Tasks for disparity statistics
tasks = {'Walk','Sando'};

% Define disparity sampling distributions (circular PD and KSDs)
sampling = {'Circ','V1','V2','MT'};

% % open fitting figure
% figRMSE = figure;

% Loop over each KSD
for s = 1:length(sampling)
    
    % Load in disparity histograms from BORIS image set
    sndo = load([imStatsDir,'dispHist' sampling{s} '_sando.mat']);
    wlk  = load([imStatsDir,'dispHist' sampling{s} '_walking.mat']);
    
    if s == 1
        Sando(s,:)  = sndo.dispHistCirc./sum(sndo.dispHistCirc);
        Walk(s,:)   = wlk.dispHistCirc./sum(wlk.dispHistCirc);
    elseif s == 2
        Sando(s,:)  = sndo.dispHistV1./sum(sndo.dispHistV1);    % sando scene stats
        Walk(s,:)   = wlk.dispHistV1./sum(wlk.dispHistV1);      % walking scene stats
        FI(s,:)     = dat.v1;                                   % V1 FI
        lfi(s,:,:)  = dat.lv1';                                  % bootstrapped V1 FI samples
        clr = 'b';
    elseif s == 3
        Sando(s,:)  = sndo.dispHistV2./sum(sndo.dispHistV2);
        Walk(s,:)   = wlk.dispHistV2./sum(wlk.dispHistV2);
        FI(s,:)     = dat.v2;
        lfi(s,:,:)  = dat.lv2';
        clr = 'g';
    elseif s == 4
        Sando(s,:)  = sndo.dispHistMT./sum(sndo.dispHistMT);
        Walk(s,:)   = wlk.dispHistMT./sum(wlk.dispHistMT);
        FI(s,:)     = dat.mt;
        lfi(s,:,:)  = dat.lmt';
        clr = 'r';
    end
    
    % If using the resampled V1, and MT data, all disparity stats should match V2
    if strcmp(flag,'_resampledV1MT') && s > 1
        sndo        = load([imStatsDir,'dispHistV2_sando.mat']);
        wlk         = load([imStatsDir,'dispHistV2_walking.mat']);
        Sando(s,:)  = sndo.dispHistV2./sum(sndo.dispHistV2);
        Walk(s,:)   = wlk.dispHistV2./sum(wlk.dispHistV2);
    end
    
    % Calcuate expected FI distributions under idealized infomax and discrimax
    % optimizations
    Sando_infomax(s,:)   = (Sando(s,:).^2)./sum(Sando(s,:).^2);
    Sando_discrimax(s,:) = (Sando(s,:).^0.5)./sum(Sando(s,:).^0.5);
    
    Walk_infomax(s,:)    = (Walk(s,:).^2)./sum(Walk(s,:).^2);
    Walk_discrimax(s,:)  = (Walk(s,:).^0.5)./sum(Walk(s,:).^0.5);

    % Only calculate best fitting power law using the neurally-defined KSDs
    if s > 1 

        % Define powers to test
        ps = linspace(0.1,3,200);

        % Loop over set of powers
        for p = 1:length(ps)

            Sando_pow = (Sando(s,:).^ps(p))./sum(Sando(s,:).^ps(p));
            Sando_err(p) = mean(abs(Sando_pow-FI(s,:))); %sqrt(mean((Sando_pow-FI(s,:)).^2));
            
            Walk_pow = (Walk(s,:).^ps(p))./sum(Walk(s,:).^ps(p));
            Walk_err(p) = mean(abs(Walk_pow-FI(s,:))); %sqrt(mean((Walk_pow-FI(s,:)).^2));
            
        end
        
        % Find power that most closely matches FI distribution and
        % disparity stats
        Sando_err_min(s) = ps((Sando_err == min(Sando_err)));
        Walk_err_min(s)  = ps((Walk_err == min(Walk_err)));
        
%         % Plot error function
%         figure(figRMSE); 
%         hold on; 
% 
%         plot(ps,Sando_err,'-','color',ColorIt(clr),'linewidth',2);
%         plot(ps,Walk_err,'--','color',ColorIt(clr),'linewidth',2);
        
        % Repeat with each of the bootstrapped disparity samples
    
        % number bootstraps
        nboots = size(lfi,3);
        
        % powers to test
        ps = linspace(0.1,4,200);
        
        % for each boot
        for n = 1:nboots
            
            % for each power
            for p = 1:length(ps)
                
                Sando_pow = (Sando(s,:).^ps(p))./sum(Sando(s,:).^ps(p));
                Sando_err(p) = mean(abs(Sando_pow-lfi(s,:,n))); %sqrt(mean((Sando_pow-FI(s,:)).^2));
                
                Walk_pow = (Walk(s,:).^ps(p))./sum(Walk(s,:).^ps(p));
                Walk_err(p) = mean(abs(Walk_pow-lfi(s,:,n))); %sqrt(mean((Walk_pow-FI(s,:)).^2));
                
            end
            
            % find minimum
            Sando_err_min_boot(s,n) = ps((Sando_err == min(Sando_err)));
            Walk_err_min_boot(s,n) = ps((Walk_err == min(Walk_err)));
            
        end

    end

    % Fit generalized Gaussian to scene probabilities, raw, squared, and halved
    params_Sando(s,:)           = fitGeneralizedLaplacian( cntr_disp, Sando(s,:) );
    params_Sando_infomax(s,:)   = fitGeneralizedLaplacian( cntr_disp, Sando_infomax(s,:) );
    params_Sando_discrimax(s,:) = fitGeneralizedLaplacian( cntr_disp, Sando_discrimax(s,:) );
    
    params_Walk(s,:)           = fitGeneralizedLaplacian( cntr_disp, Walk(s,:) );
    params_Walk_infomax(s,:)   = fitGeneralizedLaplacian( cntr_disp, Walk_infomax(s,:) );
    params_Walk_discrimax(s,:) = fitGeneralizedLaplacian( cntr_disp, Walk_discrimax(s,:) );
    
end


%% Plots

% % Add legend and save RMSE plot
% figure(figRMSE); 
% hold on;
% legend('Sando v1','Walk v1','Sando v2','Walk v2','Sando mt','Walk mt','location','northeastoutside');
% title('RMSE with FI')
 

%---------------------------------------------------------%
% Plot scene data histograms for the 3 area sampling regimes
% Food preparation
f2 = figure; 
f2.Position = [100 100 650 600];
hold on;

plot(cntr_disp, Sando(2,:),'color',ColorIt('b'),'linewidth',2);
plot(cntr_disp, Sando(3,:),'color',ColorIt('g'),'linewidth',2);
plot(cntr_disp, Sando(4,:),'color',ColorIt('r'),'linewidth',2);
%ylim([0 max(v1*1.1)]);
axis square; box on;
xlabel('horizontal disparity (deg)'); 
ylabel('Probability density');
legend('V1','V2','MT');

% Navigation
f3 = figure; 
f3.Position = [700 100 650 600];
hold on;

plot(cntr_disp, Walk(2,:),'color',ColorIt('b'),'linewidth',2);
plot(cntr_disp, Walk(3,:),'color',ColorIt('g'),'linewidth',2);
plot(cntr_disp, Walk(4,:),'color',ColorIt('r'),'linewidth',2);
%ylim([0 max(v1*1.1)]);
axis square; box on;
xlabel('horizontal disparity (deg)'); 
ylabel('Probability density');
legend('V1','V2','MT');


% %---------------------------------------------------------%
% % histogram of bootstrap errors
% f4 = figure; 
% f4.Position = [100 700 650 600];
% hold on;
% 
% for s = 2:4
%     subplot(2,3,s-1); 
%     hold on; 
%     
%     histogram(Sando_err_min_boot(s,:))
%     xlabel('best fit power');
%     title(['Sando bootstraps ' sampling(s)]);
%     
%     subplot(2,3,s-1+3); 
%     hold on; 
%     
%     histogram(Walk_err_min_boot(s,:))
%     xlabel('best fit power');
%     title(['Walk bootstraps ' sampling(s)]);
% end


%---------------------------------------------------------%
% Plot best fit power with error bars from disparity bootstrapping
% Food preparation
f5 = figure; 
f5.Position = [700 700 650 600];
hold on;

b = bar([1 2 3],Sando_err_min(2:4));
b.FaceColor = 'flat';
b.CData(1,:) = ColorIt('b');
b.CData(2,:) = ColorIt('g');
b.CData(3,:) = ColorIt('r');
er = errorbar([1 2 3],[Sando_err_min(2:4)],...
     [std(Sando_err_min_boot(2,:)) std(Sando_err_min_boot(3,:)) std(Sando_err_min_boot(4,:))],...
     [std(Sando_err_min_boot(2,:)) std(Sando_err_min_boot(3,:)) std(Sando_err_min_boot(4,:))]); 
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylabel('power');
set(gca,'xtick',[1 2 3],'xticklabel',{'V1','V2','MT'});
box on;

% Navigation
f6 = figure; 
f6.Position = [700 700 650 600];
hold on;

b = bar([1 2 3],Walk_err_min(2:4));
b.FaceColor = 'flat';
b.CData(1,:) = ColorIt('b');
b.CData(2,:) = ColorIt('g');
b.CData(3,:) = ColorIt('r');
er = errorbar([1 2 3],[Walk_err_min(2:4)],...
     [std(Walk_err_min_boot(2,:)) std(Walk_err_min_boot(3,:)) std(Walk_err_min_boot(4,:))],...
     [std(Walk_err_min_boot(2,:)) std(Walk_err_min_boot(3,:)) std(Walk_err_min_boot(4,:))]); 
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylabel('power');
set(gca,'xtick',[1 2 3],'xticklabel',{'V1','V2','MT'});
box on;


%---------------------------------------------------------%
% PLot a histogram of the best fit powers, fit a Gaussian to that
% histogram, then determine power CIs based on that Gaussian
for t = 1:2
    
    task = tasks{t};
    
    switch task
        case 'Walk'
            bins = linspace(0.25,1,15);
            err_min_boot = Walk_err_min_boot;
            
        case 'Sando'
            bins = linspace(0.5,4.5,15);
            err_min_boot = Sando_err_min_boot;
            
    end
    
    % Plot histograms of best fit powers
    f7{t} = figure; 
    hold on;

    h1 = histogram(err_min_boot(2,:),bins,'Normalization','pdf','FaceAlpha',1);
    h2 = histogram(err_min_boot(3,:),bins,'Normalization','pdf','FaceAlpha',1);
    h3 = histogram(err_min_boot(4,:),bins,'Normalization','pdf','FaceAlpha',1);
    
    h1.FaceColor = ColorIt('b');
    h2.FaceColor = ColorIt('g');
    h3.FaceColor = ColorIt('r');
    
    box on;
    
    % Fit Gaussians to power histograms and add to plot
    plt_points = linspace(min(bins),max(bins),200);
    
    [muHat_V1,sigmaHat_V1,muCI_V1,sigmaCI_V1] = normfit(err_min_boot(2,:));
    plot(plt_points,normpdf(plt_points,muHat_V1,sigmaHat_V1),'color',ColorIt('b'),'linewidth',3);
    
    [muHat_V2,sigmaHat_V2,muCI_V2,sigmaCI_V2] = normfit(err_min_boot(3,:));
    plot(plt_points,normpdf(plt_points,muHat_V2,sigmaHat_V2),'color',ColorIt('g'),'linewidth',3);
    
    [muHat_MT,sigmaHat_MT,muCI_MT,sigmaCI_MT] = normfit(err_min_boot(4,:));
    plot(plt_points,normpdf(plt_points,muHat_MT,sigmaHat_MT),'color',ColorIt('r'),'linewidth',3);
    
    xlim([min(bins) max(bins)]);
    
    % Report Mus and Cohen's D
%     display([ task ' muHat V1 : ' num2str(muHat_V1)]);
%     display([ task ' muHat V2 : ' num2str(muHat_V2)]);
%     display([ task ' muHat MT : ' num2str(muHat_MT)]);
%     
    D_V1_V2 = (muHat_V1 - muHat_V2) / sqrt( ( (nboots-1)*(sigmaHat_V1^2) + (nboots-1)*(sigmaHat_V2^2) ) / (2*nboots - 2) );
    D_V2_MT = (muHat_V2 - muHat_MT) / sqrt( ( (nboots-1)*(sigmaHat_V2^2) + (nboots-1)*(sigmaHat_MT^2) ) / (2*nboots - 2) );
    D_V1_MT = (muHat_V1 - muHat_MT) / sqrt( ( (nboots-1)*(sigmaHat_V1^2) + (nboots-1)*(sigmaHat_MT^2) ) / (2*nboots - 2) );
%     
%     display([ task ' D V1-V2 : ' num2str(D_V1_V2)]);
%     display([ task ' D V2-MT : ' num2str(D_V2_MT)]);
%     display([ task ' D V1-MT : ' num2str(D_V1_MT)]);
    
    %---------------------------------------------------------%
    % Plot mean power and CIs as a bar plot
    f8{t} = figure; 
    hold on;

    b = bar([1 2 3],[muHat_V1 muHat_V2 muHat_MT]);
    b.FaceColor = 'flat';
    b.CData(1,:) = ColorIt('b');
    b.CData(2,:) = ColorIt('g');
    b.CData(3,:) = ColorIt('r');
    er = errorbar([1 2 3],[muHat_V1 muHat_V2 muHat_MT],...
        abs([muHat_V1-muCI_V1(1) muHat_V2-muCI_V2(1) muHat_MT-muCI_MT(1)]),...
        abs([muHat_V1-muCI_V1(2) muHat_V2-muCI_V2(2) muHat_MT-muCI_MT(2)]));
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    ylabel('power');
    set(gca,'xtick',[1 2 3],'xticklabel',{'V1','V2','MT'});
    box on;  
    
end

%---------------------------------------------------------%
% Overlay FI with best fit power
% Food preparation
f9 = figure; 
hold on;

h(1) = plot(cntr_disp, FI(2,:),'color',ColorIt('b'),'linewidth',2);
plot(cntr_disp, FI(3,:),'color',ColorIt('g'),'linewidth',2);
plot(cntr_disp, FI(4,:),'color',ColorIt('r'),'linewidth',2);

h(2) = plot(cntr_disp, (Sando(2,:).^Sando_err_min(2))./sum(Sando(2,:).^Sando_err_min(2)),':','color',ColorIt('b'),'linewidth',2);
plot(cntr_disp, (Sando(3,:).^Sando_err_min(3))./sum(Sando(3,:).^Sando_err_min(3)),':','color',ColorIt('g'),'linewidth',2);
plot(cntr_disp, (Sando(4,:).^Sando_err_min(4))./sum(Sando(4,:).^Sando_err_min(4)),':','color',ColorIt('r'),'linewidth',2);
%ylim([0 max(v1*1.1)]);
axis square; box on;
xlabel('horizontal disparity (deg)'); ylabel('Probability density');
legend(h,'FI','best fit p');


% Navigation
f10 = figure; 
hold on;

h(1) = plot(cntr_disp, FI(2,:),'color',ColorIt('b'),'linewidth',2);
plot(cntr_disp, FI(3,:),'color',ColorIt('g'),'linewidth',2);
plot(cntr_disp, FI(4,:),'color',ColorIt('r'),'linewidth',2);

h(2) = plot(cntr_disp, (Walk(2,:).^Walk_err_min(2))./sum(Walk(2,:).^Walk_err_min(2)),':','color',ColorIt('b'),'linewidth',2);
plot(cntr_disp, (Walk(3,:).^Walk_err_min(3))./sum(Walk(3,:).^Walk_err_min(3)),':','color',ColorIt('g'),'linewidth',2);
plot(cntr_disp, (Walk(4,:).^Walk_err_min(4))./sum(Walk(4,:).^Walk_err_min(4)),':','color',ColorIt('r'),'linewidth',2);
%ylim([0 max(v1*1.1)]);
axis square; box on;
xlabel('horizontal disparity (deg)'); ylabel('Probability density');
legend(h,'FI','best fit p');


%% Save figures
if saveOn

%     % F1: RMSE between disparity distribution^p and empirical FI dists
%     saveas(figRMSE,['./plots/NDS_comparison/FINAL_RMSE_power_fits',flag,'.eps'],'epsc');

    % F2: Plots of true disparity distributions for each region (food prep)
    saveas(f2,['./plots/NDS_comparison/Disparity_p_eachRegion_Sando',flag,'.eps'],'epsc');

    % F3: Plots of true disparity distributions for each region (navigation)
    saveas(f3,['./plots/NDS_comparison/Disparity_p_eachRegion_Walk',flag,'.eps'],'epsc');

    % F5: Bar plot of best fit power law and bootstrapped errorbars (food prep)
    saveas(f5,['./plots/NDS_comparison/FINAL_power_Sando_errorbar',flag,'.eps'],'epsc'); 

    % F6: Bar plot of best fit power law and bootstrapped errorbars (navigation)
    saveas(f6,['./plots/NDS_comparison/FINAL_power_Walk_errorbar',flag,'.eps'],'epsc'); 

    % F7: Histograms and best fit Gaussians to histogram
    saveas(f7{1},['./plots/NDS_comparison/FINAL_power_Walk_histogram',flag,'.eps'],'epsc');
    saveas(f7{2},['./plots/NDS_comparison/FINAL_power_Sando_histogram',flag,'.eps'],'epsc');

    % F8: Bar plot of best fit power law and CIs based on Gaussian fit to histogram
    saveas(f8{1},['./plots/NDS_comparison/FINAL_power_Walk_distributions',flag,'.eps'],'epsc');
    saveas(f8{2},['./plots/NDS_comparison/FINAL_power_Sando_distributions',flag,'.eps'],'epsc');

    % F9: Best fit power law xformed disparity distribution vs. FI (food prep)
    saveas(f9,['./plots/NDS_comparison/FINAL_power_Sando_fits',flag,'.eps'],'epsc');

    % F10: Best fit power law xformed disparity distribution vs. FI (navigation)
    saveas(f10,['./plots/NDS_comparison/FINAL_power_Walk_fits',flag,'.eps'],'epsc');

end
