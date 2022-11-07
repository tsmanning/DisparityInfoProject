% Demo of FI calculation

clear all
close all

numDisps = 1000;
disps = linspace(-2,2,numDisps);

p       = [20 60 0 0.4 0.7 0];
exGabor = p(1) + p(2)*exp( -(disps-p(3)).^2 / (2*p(4)^2) ) .* cos(2*pi*p(5)*(disps-p(3))+p(6));

dG1 = diff(exGabor(:,1:end-1),[],2);
dG2 = diff(exGabor(:,2:end),[],2);

dGabor = ((dG1 + dG2)/(2*diff(disps(1:2)))).^2;

FIarith = (dGabor.^2)./exGabor(2:end-1);

poisSupp = [0:round(max(exGabor))*2]';

noiseDist = nan(numel(poisSupp),numel(disps));

for ii = 1:numDisps
    
    noiseDist(:,ii) = poisspdf(poisSupp,exGabor(ii));

end


d1 = diff(log(noiseDist(:,1:end-1)),[],2);
d2 = diff(log(noiseDist(:,2:end)),[],2);

noiseDeriv = ((d1 + d2)/(2*diff(disps(1:2)))).^2;

FI = sum(noiseDeriv.*noiseDist(:,2:end-1));

%plot 
f1 = figure;
f1.Position = [100 100 750 600];

imagesc(disps,poisSupp,noiseDist); axis xy

colormap gray
cb1 = colorbar;
cb1.Label.String = 'p(m|\theta)';

set(gca,'PlotBoxAspectRatio',[1 1 1],'fontsize',20);
xlabel('\theta (disparity)');
ylabel('m');

f2 = figure;
f2.Position = [600 100 750 600];

imagesc(disps(2:end-1),poisSupp,noiseDeriv); axis xy

colormap gray
cb2 = colorbar;
cb2.Label.String = 'd/d\theta ln p(m|\theta)^2';

set(gca,'PlotBoxAspectRatio',[1 1 1],'fontsize',20);
xlabel('\theta (disparity)');
ylabel('m');

f3 = figure;
f3.Position = [1200 100 650 600];
hold on

plot(disps(2:end-1),FI/max(FI),'k','linewidth',2);
plot(disps(2:end-1),FIarith/max(FIarith),'r','linewidth',2);

set(gca,'PlotBoxAspectRatio',[1 1 1],'fontsize',20);
xlabel('\theta (disparity)');
ylabel('FI');