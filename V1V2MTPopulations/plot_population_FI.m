function plot_population_FI(V1,V2,MT,flag,iterations)

%% FISHER INFORMATION

% plot all cells as heat map
figure; hold on;
subplot(1,3,1); hold on; title('V1');
imagesc(V1.FI); colorbar; axis tight;
subplot(1,3,2); hold on; title('V2');
imagesc(V2.FI); colorbar; axis tight;
subplot(1,3,3); hold on; title('MT');
imagesc(MT.FI); colorbar; axis tight;
saveas(gcf,['./plots/AnalyzeFI/FI_AllNeurons' flag '.png']); 

% FI as log
figure; hold on;
subplot(1,3,1); hold on; title('V1');
imagesc(log(V1.FI + eps)); colorbar; axis tight;
subplot(1,3,2); hold on; title('V2');
imagesc(log(V2.FI + eps)); colorbar; axis tight;
subplot(1,3,3); hold on; title('MT');
imagesc(log(MT.FI + eps)); colorbar; axis tight;
saveas(gcf,['./plots/AnalyzeFI/FI_AllNeurons_log' flag '.png']); 

% sqrt FI
figure; hold on;
subplot(1,3,1); hold on; title('V1');
imagesc(sqrt(V1.FI)); colorbar; axis tight;
subplot(1,3,2); hold on; title('V2');
imagesc(sqrt(V2.FI)); colorbar; axis tight;
subplot(1,3,3); hold on; title('MT');
imagesc(sqrt(MT.FI)); colorbar; axis tight;
saveas(gcf,['./plots/AnalyzeFI/FI_AllNeurons_sqrt' flag '.png']); 

% tuning curves
figure; hold on;
subplot(1,3,1); hold on; title('V1');
imagesc(V1.allresp); colorbar; axis tight;
subplot(1,3,2); hold on; title('V2');
imagesc(V2.allresp); colorbar; axis tight;
subplot(1,3,3); hold on; title('MT');
imagesc(MT.allresp); colorbar; axis tight;
saveas(gcf,['./plots/AnalyzeFI/FI_AllNeurons_tuningcurves' flag '.png']); 

% visualize population FI - RAW

% disparities corresponding to each column in FI matrices
x = [-2 : 0.01 : 2];

% downsample FI before plotting/fitting to match scene stats
if(0)
    
    % resolution used in Tyler's scene stats
    res = 52;
    lb  = -2;
    ub  = 2;edges_disp = linspace(lb,ub,res);
    cntr_disp  = edges_disp(1:end-1) + diff(edges_disp(1:2))/2;
    
    % find indices of x that intersect with x
    [y,ia,ib] = intersect(round(x,2),round(cntr_disp,2));
    
    % grab just those indices
    x = x(ia);
    V1.FI = V1.FI(:,ia);
    V2.FI = V2.FI(:,ia);
    MT.FI = MT.FI(:,ia);
    
end


% sum it up to get population FI
v1          = sum(V1.FI); v1 = v1 / sum(v1);
v2          = sum(V2.FI); v2 = v2 / sum(v2);
mt          = sum(MT.FI); mt = mt / sum(mt);

% PLOTS
figure; hold on; title('Raw FI');
plot(x, v1,'color',ColorIt('b'),'linewidth',2);
plot(x, v2,'color',ColorIt('g'),'linewidth',2);
plot(x, mt,'color',ColorIt('r'),'linewidth',2);
legend('V1','V2','MT');
axis square; box on;
xlabel('horizontal disparity (deg)'); ylabel('normalized FI');
saveas(gcf,['./plots/AnalyzeFI/FI_AllNorm' flag '.png']); 
saveas(gcf,['./plots/AnalyzeFI/FI_AllNorm' flag '.eps'],'epsc');

% plot data and laplacian fits
figure; hold on; title('Raw FI with Laplacian fits');
hp(1) = plot(x, v1,'color',ColorIt('b'),'linewidth',1);
hp(2) = plot(x, v2,'color',ColorIt('g'),'linewidth',1);
hp(3) = plot(x, mt,'color',ColorIt('r'),'linewidth',1);

legend('V1','V2','MT');

glV1            = fitGeneralizedLaplacian( x, v1 );
glV2            = fitGeneralizedLaplacian( x, v2 );
glMT            = fitGeneralizedLaplacian( x, mt );

gl = exp( -(abs(x)/glV1(1)).^glV1(2) ); gl = gl / sum(gl);
plot( x, gl, '--','color',ColorIt('b'),'linewidth',2 );
gl = exp( -(abs(x)/glV2(1)).^glV2(2) ); gl = gl / sum(gl);
plot( x, gl, '--','color',ColorIt('g'),'linewidth',2 );
gl = exp( -(abs(x)/glMT(1)).^glMT(2) ); gl = gl / sum(gl);
plot( x, gl, '--','color',ColorIt('r'),'linewidth',2 );

saveas(gcf,['./plots/AnalyzeFI/FI_rawwithfits' flag '.png']); 
saveas(gcf,['./plots/AnalyzeFI/FI_rawwithfits' flag '.eps'],'epsc'); 


% plot with error visualized with samples and fit laplacian to samples

iters   = iterations;
sams    = 400; % so the number of samples can be matched across regions

for n = 1:iters
    
    if ~mod(n,50)
        display(['running bootstraps iteration ' num2str(n)]);
    end
    
    % sample V1 with replacement
    littleV1FI          = V1.FI(randsample(size(V1.FI,1),sams,true),:);
    lv1(n,:)            = sum(littleV1FI); 
    lv1(n,:)            = lv1(n,:) / sum(lv1(n,:));
    params_V1All(n,:)   = fitGeneralizedLaplacian( x, lv1(n,:) );

    % sample V2
    littleV2FI          = V2.FI(randsample(size(V2.FI,1),sams,true),:);
    lv2(n,:)            = sum(littleV2FI); 
    lv2(n,:)            = lv2(n,:) / sum(lv2(n,:));
    params_V2All(n,:)   = fitGeneralizedLaplacian( x, lv2(n,:) );

    % sample MT
    littleMTFI          = MT.FI(randsample(size(MT.FI,1),sams,true),:);
    lmt(n,:)            = sum(littleMTFI); 
    lmt(n,:)            = lmt(n,:) / sum(lmt(n,:));
    params_MTAll(n,:)   = fitGeneralizedLaplacian( x, lmt(n,:) );
    
end

% compute variance
params_V1All_var         = (params_V1All(:,1).^2 .*gamma(3./params_V1All(:,2)) ) ./gamma(1./params_V1All(:,2));
params_V2All_var         = (params_V2All(:,1).^2 .*gamma(3./params_V2All(:,2)) ) ./gamma(1./params_V2All(:,2));
params_MTAll_var         = (params_MTAll(:,1).^2 .*gamma(3./params_MTAll(:,2)) ) ./gamma(1./params_MTAll(:,2));

% compute excess kurtosis
params_V1All_exkurt         = ( (gamma(5./params_V1All(:,2)).*gamma(1./params_V1All(:,2))) ./ (gamma(3./params_V1All(:,2)).^2) ) - 3;
params_V2All_exkurt         = ( (gamma(5./params_V2All(:,2)).*gamma(1./params_V2All(:,2))) ./ (gamma(3./params_V2All(:,2)).^2) ) - 3;
params_MTAll_exkurt         = ( (gamma(5./params_MTAll(:,2)).*gamma(1./params_MTAll(:,2))) ./ (gamma(3./params_MTAll(:,2)).^2) ) - 3;

% compute entropy
params_V1All_entropy         = 1./params_V1All(:,2) - log( params_V1All(:,2)./ (2*params_V1All(:,1).*gamma(1./params_V1All(:,2))) );
params_V2All_entropy         = 1./params_V2All(:,2) - log( params_V2All(:,2)./ (2*params_V2All(:,1).*gamma(1./params_V2All(:,2))) );
params_MTAll_entropy         = 1./params_MTAll(:,2) - log( params_MTAll(:,2)./ (2*params_MTAll(:,1).*gamma(1./params_MTAll(:,2))) );


% save it all for future use
save(['results_GaussFit' flag '.mat'],'x','v1','v2','mt',...
    'lv1','params_V1All','params_V1All_var','params_V1All_exkurt','params_V1All_entropy',...
    'lv2','params_V2All','params_V2All_var','params_V2All_exkurt','params_V2All_entropy',...
    'lmt','params_MTAll','params_MTAll_var','params_MTAll_exkurt','params_MTAll_entropy');

% plot all the lines
figure; hold on; 
plot(x, lv1, 'color', [ColorIt('b') 0.1],'linewidth',0.25);
plot(x, lv2, 'color', [ColorIt('g') 0.1],'linewidth',0.25);
plot(x, lmt, 'color', [ColorIt('r') 0.1],'linewidth',0.25);

plot(x, v1,'color',ColorIt('b'),'linewidth',2);
plot(x, v2,'color',ColorIt('g'),'linewidth',2);
plot(x, mt,'color',ColorIt('r'),'linewidth',2);
ylim([0 max(v1*1.1)]);
axis square; box on;
xlabel('horizontal disparity (deg)'); ylabel('Fisher Information (normalized)');
saveas(gcf,['./plots/AnalyzeFI/FI_withSampling' flag '.png']); 
saveas(gcf,['./plots/AnalyzeFI/FI_withSampling' flag '.eps'],'epsc'); 

% plot with shaded error bar

% bar plots of variance, excess kurtosis, and entropy
figure; hold on;
b = bar([1 2 3],[mean(params_V1All_var) mean(params_V2All_var) mean(params_MTAll_var)]);
b.FaceColor = 'flat';
b.CData(1,:) = ColorIt('b');
b.CData(2,:) = ColorIt('g');
b.CData(3,:) = ColorIt('r');
er = errorbar([1 2 3],[mean(params_V1All_var) mean(params_V2All_var) mean(params_MTAll_var)],...
    [std(params_V1All_var) std(params_V2All_var) std(params_MTAll_var)],[std(params_V1All_var) std(params_V2All_var) std(params_MTAll_var)]); 
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylabel('variance');
set(gca,'xtick',[1 2 3],'xticklabel',{'V1','V2','MT'});
box on;
saveas(gcf,['./plots/AnalyzeFI/FI_variance' flag '.png']); 
saveas(gcf,['./plots/AnalyzeFI/FI_variance' flag '.eps'],'epsc'); 

figure; hold on;
b = bar([1 2 3],[mean(params_V1All_exkurt) mean(params_V2All_exkurt) mean(params_MTAll_exkurt)]);
b.FaceColor = 'flat';
b.CData(1,:) = ColorIt('b');
b.CData(2,:) = ColorIt('g');
b.CData(3,:) = ColorIt('r');
er = errorbar([1 2 3],[mean(params_V1All_exkurt) mean(params_V2All_exkurt) mean(params_MTAll_exkurt)],...
    [std(params_V1All_exkurt) std(params_V2All_exkurt) std(params_MTAll_exkurt)],[std(params_V1All_exkurt) std(params_V2All_exkurt) std(params_MTAll_exkurt)]); 
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylabel('excess kurtosis');
set(gca,'xtick',[1 2 3],'xticklabel',{'V1','V2','MT'});
box on;
saveas(gcf,['./plots/AnalyzeFI/FI_exkurt' flag '.png']); 
saveas(gcf,['./plots/AnalyzeFI/FI_exkurt' flag '.eps'],'epsc'); 

figure; hold on;
b = bar([1 2 3],[mean(params_V1All_entropy) mean(params_V2All_entropy) mean(params_MTAll_entropy)]);
b.FaceColor = 'flat';
b.CData(1,:) = ColorIt('b');
b.CData(2,:) = ColorIt('g');
b.CData(3,:) = ColorIt('r');
er = errorbar([1 2 3],[mean(params_V1All_entropy) mean(params_V2All_entropy) mean(params_MTAll_entropy)],...
    [std(params_V1All_entropy) std(params_V2All_entropy) std(params_MTAll_entropy)],[std(params_V1All_entropy) std(params_V2All_entropy) std(params_MTAll_entropy)]); 
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylabel('entropy');
set(gca,'xtick',[1 2 3],'xticklabel',{'V1','V2','MT'});
box on;
saveas(gcf,['./plots/AnalyzeFI/FI_entropy' flag '.png']); 
saveas(gcf,['./plots/AnalyzeFI/FI_entropy' flag '.eps'],'epsc'); 
