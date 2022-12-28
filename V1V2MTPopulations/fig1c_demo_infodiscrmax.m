clear all; close all;

res = 52;
lb  = -2;
ub  = 2;edges_disp = linspace(lb,ub,res);
cntr_disp  = edges_disp(1:end-1) + diff(edges_disp(1:2))/2;

s = 0.3;
p = 1;

p1 = exp( -(abs(cntr_disp)/s).^p );
p1 = p1 / sum(p1);

p2 = p1.^2;
p2 = p2 / sum(p2);

ppt5 = p1.^0.5;
ppt5 = ppt5 / sum(ppt5);

% plot probabilities squared and to 1/2
figure; hold on;
plot(cntr_disp, p1,'-','color',ColorIt('k'),'linewidth',2);
plot(cntr_disp, p2,'--','color',ColorIt('k'),'linewidth',2);
plot(cntr_disp, ppt5,':','color',ColorIt('k'),'linewidth',2);
axis square; box on;
ylim([0 .27])
legend('p=1','infomax','discrimax');
xlabel('horizontal disparity (deg)'); ylabel('predicted FI');
saveas(gcf,['./plots/NDS_comparison/Disparity_demo_panel.eps'],'epsc');