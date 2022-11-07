clear all;
close all;

flag = '';
%flag = '_resampledV1MT';

% load in FI data
dat = load(['results_GaussFit' flag '.mat']);

% resolution used in Tyler's scene stats
res = 52;
lb  = -2;
ub  = 2;edges_disp = linspace(lb,ub,res);
cntr_disp  = edges_disp(1:end-1) + diff(edges_disp(1:2))/2;

% load in Tyler's scene stats analysis
sampling = {'Circ','V1','V2','MT'};

% tasks for disparity statistics
tasks = {'Walk','Sando'};

% open fitting figure
figRMSE = figure;

for s = 1:length(sampling)
    
    % load data
    sndo = load(['/Users/emily/Dropbox/FisherDisparityProject/FisherDisparityProjectShared/SceneStatsAnalysis/DataFromTyler/dispHist' sampling{s} '_sando.mat']);
    wlk = load(['/Users/emily/Dropbox/FisherDisparityProject/FisherDisparityProjectShared/SceneStatsAnalysis/DataFromTyler/dispHist' sampling{s} '_walking.mat']);
    %sndo = load(['C:/Users/alexa/Dropbox/FisherDisparityProjectShared/SceneStatsAnalysis/DataFromTyler/dispHist' sampling{s} '_sando.mat']);
    %wlk = load(['C:/Users/alexa/Dropbox/FisherDisparityProjectShared/SceneStatsAnalysis/DataFromTyler/dispHist' sampling{s} '_walking.mat']);
    
    
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
    
    % if using the resampled V1, and MT data, all disparity stats should match V2
    if strcmp(flag,'_resampledV1MT') && s > 1
        sndo = load(['/Users/emily/Dropbox/FisherDisparityProject/FisherDisparityProjectShared/SceneStatsAnalysis/DataFromTyler/dispHistV2_sando.mat']);
        wlk = load(['/Users/emily/Dropbox/FisherDisparityProject/FisherDisparityProjectShared/SceneStatsAnalysis/DataFromTyler/dispHistV2_walking.mat']);
        Sando(s,:)  = sndo.dispHistV2./sum(sndo.dispHistV2);
        Walk(s,:)   = wlk.dispHistV2./sum(wlk.dispHistV2);
    end
        
    
    % calcuate infomax and discrimax
    Sando_infomax(s,:) = (Sando(s,:).^2)./sum(Sando(s,:).^2);
    Sando_discrimax(s,:) = (Sando(s,:).^0.5)./sum(Sando(s,:).^0.5);
    
    Walk_infomax(s,:) = (Walk(s,:).^2)./sum(Walk(s,:).^2);
    Walk_discrimax(s,:) = (Walk(s,:).^0.5)./sum(Walk(s,:).^0.5);
    
    % find the power law that, applied to each, best matches the population FI
    
    % powers to test
    ps = linspace(0.1,3,200);
    
    if s > 1
        
        % for each power
        for p = 1:length(ps)

            Sando_pow = (Sando(s,:).^ps(p))./sum(Sando(s,:).^ps(p));
            Sando_err(p) = mean(abs(Sando_pow-FI(s,:))); %sqrt(mean((Sando_pow-FI(s,:)).^2));
            
            Walk_pow = (Walk(s,:).^ps(p))./sum(Walk(s,:).^ps(p));
            Walk_err(p) = mean(abs(Walk_pow-FI(s,:))); %sqrt(mean((Walk_pow-FI(s,:)).^2));
            
        end
        
        % find minimum
        Sando_err_min(s) = ps((Sando_err == min(Sando_err)));
        Walk_err_min(s) = ps((Walk_err == min(Walk_err)));
        
        % plot error function
        figure(figRMSE); hold on; title('RMSE with FI')
        plot(ps,Sando_err,'-','color',ColorIt(clr),'linewidth',2);
        plot(ps,Walk_err,'--','color',ColorIt(clr),'linewidth',2);
        
        % repeat with bootstrapped sampling
    
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
    

    % fit generalized Gaussian to scene probabilities, raw, squared, and halved
    params_Sando(s,:)           = fitGeneralizedLaplacian( cntr_disp, Sando(s,:) );
    params_Sando_infomax(s,:)   = fitGeneralizedLaplacian( cntr_disp, Sando_infomax(s,:) );
    params_Sando_discrimax(s,:) = fitGeneralizedLaplacian( cntr_disp, Sando_discrimax(s,:) );
    
    params_Walk(s,:)           = fitGeneralizedLaplacian( cntr_disp, Walk(s,:) );
    params_Walk_infomax(s,:)   = fitGeneralizedLaplacian( cntr_disp, Walk_infomax(s,:) );
    params_Walk_discrimax(s,:) = fitGeneralizedLaplacian( cntr_disp, Walk_discrimax(s,:) );
    
    if(0)
    % plot probabilities
    figure; hold on;
    plot(cntr_disp, Sando(s,:),'color',ColorIt('k'),'linewidth',2);
    axis square; box on;
    xlabel('horizontal disparity (deg)'); ylabel('probability');
    saveas(gcf,['./plots/NDS_comparison/' flag 'Disparity_Sando_' sampling{s} '.eps'],'epsc');
    
    figure; hold on;
    plot(cntr_disp, Walk(s,:),'color',ColorIt('k'),'linewidth',2);
    axis square; box on;
    xlabel('horizontal disparity (deg)'); ylabel('probability');
    saveas(gcf,['./plots/NDS_comparison/' flag 'Disparity_Walk_' sampling{s} '.eps'],'epsc');
    
    % plot probabilities squared and to 1/2
    figure; hold on;
    plot(cntr_disp, Sando(s,:),'-','color',ColorIt('k'),'linewidth',2);
    plot(cntr_disp, Sando_infomax(s,:),'--','color',ColorIt('k'),'linewidth',2);
    plot(cntr_disp, Sando_discrimax(s,:),':','color',ColorIt('k'),'linewidth',2);
    axis square; box on;
    legend('p=1','infomax','discrimax');
    xlabel('horizontal disparity (deg)'); ylabel('predicted FI');
    saveas(gcf,['./plots/NDS_comparison/' flag 'Disparity_optimal_Sando_' sampling{s} '.eps'],'epsc');
    
    figure; hold on;
    plot(cntr_disp, Walk(s,:),'-','color',ColorIt('k'),'linewidth',2);
    plot(cntr_disp, Walk_infomax(s,:),'--','color',ColorIt('k'),'linewidth',2);
    plot(cntr_disp, Walk_discrimax(s,:),':','color',ColorIt('k'),'linewidth',2);
    axis square; box on;
    legend('p=1','infomax','discrimax');
    xlabel('horizontal disparity (deg)'); ylabel('predicted FI');
    saveas(gcf,['./plots/NDS_comparison/' flag 'Disparity_optimal_Walk_' sampling{s} '.eps'],'epsc');
    end
    
end

% add legend and save RMSE plot
figure(figRMSE); hold on;
legend('Sando v1','Walk v1','Sando v2','Walk v2','Sando mt','Walk mt','location','northeastoutside');
saveas(gcf,['./plots/NDS_comparison/' flag 'FINAL_RMSE_power_fits.eps'],'epsc'); 


% plot scene data histograms for the 3 area sampling regimes
figure; hold on;
plot(cntr_disp, Sando(2,:),'color',ColorIt('b'),'linewidth',2);
plot(cntr_disp, Sando(3,:),'color',ColorIt('g'),'linewidth',2);
plot(cntr_disp, Sando(4,:),'color',ColorIt('r'),'linewidth',2);
%ylim([0 max(v1*1.1)]);
axis square; box on;
xlabel('horizontal disparity (deg)'); ylabel('Probability density');
saveas(gcf,['./plots/NDS_comparison/' flag 'Disparity_p_eachRegion_Sando.eps'],'epsc');
legend('V1','V2','MT');

figure; hold on;
plot(cntr_disp, Walk(2,:),'color',ColorIt('b'),'linewidth',2);
plot(cntr_disp, Walk(3,:),'color',ColorIt('g'),'linewidth',2);
plot(cntr_disp, Walk(4,:),'color',ColorIt('r'),'linewidth',2);
%ylim([0 max(v1*1.1)]);
axis square; box on;
xlabel('horizontal disparity (deg)'); ylabel('Probability density');
saveas(gcf,['./plots/NDS_comparison/' flag 'Disparity_p_eachRegion_Walk.eps'],'epsc');
legend('V1','V2','MT');


% histogram of bootstrap errors
figure; hold on;
for s = 2:4
    subplot(2,3,s-1); hold on; title(['Sando bootstraps ' sampling(s)]);
    histogram(Sando_err_min_boot(s,:))
    xlabel('best fit power');
    
    subplot(2,3,s-1+3); hold on; title(['Walk bootstraps ' sampling(s)]);
    histogram(Walk_err_min_boot(s,:))
    xlabel('best fit power');
end
saveas(gcf,['./plots/NDS_comparison/' flag 'Bootstrap histograms_' sampling{s} '.eps'],'epsc');


% plot best fit power with bootstrapping error bars
figure; hold on;
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
saveas(gcf,['./plots/NDS_comparison/' flag 'FINAL_power_Sando_errorbar.eps'],'epsc'); 

figure; hold on;
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
saveas(gcf,['./plots/NDS_comparison/' flag 'FINAL_power_Walk_errorbar.eps'],'epsc'); 


% clean plot of best fit power histograms, gaussian fits, and CIs
for t = 1:2
    
    task = tasks{t};
    
    switch task
    
        case 'Walk'
            
            bins = linspace(0.25,1,15);
            err_min_boot = Walk_err_min_boot;
            
        case 'Sando'
            
            bins = linspace(0.5,4.5,15);
            err_min_boot = Sando_err_min_boot;
            
        otherwise
            
    end
    

    figure; hold on;
    h1 = histogram(err_min_boot(2,:),bins,'Normalization','pdf','FaceAlpha',1);
    h2 = histogram(err_min_boot(3,:),bins,'Normalization','pdf','FaceAlpha',1);
    h3 = histogram(err_min_boot(4,:),bins,'Normalization','pdf','FaceAlpha',1);
    
    h1.FaceColor = ColorIt('b');
    h2.FaceColor = ColorIt('g');
    h3.FaceColor = ColorIt('r');
    
    box on;
    
    % fit normal distribution and add to plot
    plt_points = linspace(min(bins),max(bins),200);
    
    [muHat_V1,sigmaHat_V1,muCI_V1,sigmaCI_V1] = normfit(err_min_boot(2,:));
    plot(plt_points,normpdf(plt_points,muHat_V1,sigmaHat_V1),'color',ColorIt('b'),'linewidth',3);
    
    [muHat_V2,sigmaHat_V2,muCI_V2,sigmaCI_V2] = normfit(err_min_boot(3,:));
    plot(plt_points,normpdf(plt_points,muHat_V2,sigmaHat_V2),'color',ColorIt('g'),'linewidth',3);
    
    [muHat_MT,sigmaHat_MT,muCI_MT,sigmaCI_MT] = normfit(err_min_boot(4,:));
    plot(plt_points,normpdf(plt_points,muHat_MT,sigmaHat_MT),'color',ColorIt('r'),'linewidth',3);
    
    xlim([min(bins) max(bins)]);
    
    % report Mus and Cohen's D
    display([ task ' muHat V1 : ' num2str(muHat_V1)]);
    display([ task ' muHat V2 : ' num2str(muHat_V2)]);
    display([ task ' muHat MT : ' num2str(muHat_MT)]);
    
    D_V1_V2 = (muHat_V1 - muHat_V2) / sqrt( ( (nboots-1)*(sigmaHat_V1^2) + (nboots-1)*(sigmaHat_V2^2) ) / (2*nboots - 2) );
    D_V2_MT = (muHat_V2 - muHat_MT) / sqrt( ( (nboots-1)*(sigmaHat_V2^2) + (nboots-1)*(sigmaHat_MT^2) ) / (2*nboots - 2) );
    D_V1_MT = (muHat_V1 - muHat_MT) / sqrt( ( (nboots-1)*(sigmaHat_V1^2) + (nboots-1)*(sigmaHat_MT^2) ) / (2*nboots - 2) );
    
    display([ task ' D V1-V2 : ' num2str(D_V1_V2)]);
    display([ task ' D V2-MT : ' num2str(D_V2_MT)]);
    display([ task ' D V1-MT : ' num2str(D_V1_MT)]);
    
    saveas(gcf,['./plots/NDS_comparison/' flag 'FINAL_power_' task '_histogram.eps'],'epsc');
    
    figure; hold on;
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
    saveas(gcf,['./plots/NDS_comparison/' flag 'FINAL_power_' task '_distributions.eps'],'epsc');
    
    
end


% overlay FI with best fit power
figure; hold on;

h(1) = plot(cntr_disp, FI(2,:),'color',ColorIt('b'),'linewidth',2);
plot(cntr_disp, FI(3,:),'color',ColorIt('g'),'linewidth',2);
plot(cntr_disp, FI(4,:),'color',ColorIt('r'),'linewidth',2);

h(2) = plot(cntr_disp, (Sando(2,:).^Sando_err_min(2))./sum(Sando(2,:).^Sando_err_min(2)),':','color',ColorIt('b'),'linewidth',2);
plot(cntr_disp, (Sando(3,:).^Sando_err_min(3))./sum(Sando(3,:).^Sando_err_min(3)),':','color',ColorIt('g'),'linewidth',2);
plot(cntr_disp, (Sando(4,:).^Sando_err_min(4))./sum(Sando(4,:).^Sando_err_min(4)),':','color',ColorIt('r'),'linewidth',2);
%ylim([0 max(v1*1.1)]);
axis square; box on;
xlabel('horizontal disparity (deg)'); ylabel('Probability density');
saveas(gcf,['./plots/NDS_comparison/' flag 'FINAL_power_Sando_fits.eps'],'epsc');
legend(h,'FI','best fit p');


figure; hold on;

h(1) = plot(cntr_disp, FI(2,:),'color',ColorIt('b'),'linewidth',2);
plot(cntr_disp, FI(3,:),'color',ColorIt('g'),'linewidth',2);
plot(cntr_disp, FI(4,:),'color',ColorIt('r'),'linewidth',2);

h(2) = plot(cntr_disp, (Walk(2,:).^Walk_err_min(2))./sum(Walk(2,:).^Walk_err_min(2)),':','color',ColorIt('b'),'linewidth',2);
plot(cntr_disp, (Walk(3,:).^Walk_err_min(3))./sum(Walk(3,:).^Walk_err_min(3)),':','color',ColorIt('g'),'linewidth',2);
plot(cntr_disp, (Walk(4,:).^Walk_err_min(4))./sum(Walk(4,:).^Walk_err_min(4)),':','color',ColorIt('r'),'linewidth',2);
%ylim([0 max(v1*1.1)]);
axis square; box on;
xlabel('horizontal disparity (deg)'); ylabel('Probability density');
saveas(gcf,['./plots/NDS_comparison/' flag 'FINAL_power_Walk_fits.eps'],'epsc');
legend(h,'FI','best fit p');


