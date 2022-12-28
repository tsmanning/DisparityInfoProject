% Demo of screen to retinal conversion

clear all
close all

disps = linspace(-2,2,100);

cortAreas = {'V1','MT'};
titles = {'V1/V2','MT'};

rfCents = [0 -10 -20];

colors = [1 0 0; 0.5 0 0.5; 0 0 1];

f = figure;
f.Position = [100 100 1100 600];

for ii = 1:numel(cortAreas)
    
    subplot(1,2,ii);
    hold on

    plot([-2 2],[-2 2],'--k');
    plot([0 0],[-2 2],'--k');
    plot([-2 2],[0 0],'--k');

    for jj = 1:numel(rfCents)

        
        retDisps{ii,jj} = screen2retDisp(disps,rfCents(jj),cortAreas{ii});

        p(jj) = plot(disps,retDisps{ii,jj},'color',colors(jj,:),'linewidth',2);

    end

    set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'ylim',[-2 2],'fontsize',20);
    xlabel('Screen disparity');
    ylabel('Retinal disparity');
    title(titles{ii});
    legend(p,{'0','+/-10','+/-20'},'location','southeast');

end

p       = [5 60 0 0.4 0.7 0];
exGabor = p(1) + p(2)*exp( -(disps-p(3)).^2 / (2*p(4)^2) ) .* cos(2*pi*p(5)*(disps-p(3))+p(6));

if 1 

    f2 = figure;
    f2.Position  = [100 500 1100 600];

    for ii = 1:numel(cortAreas)

        subplot(1,2,ii);
        hold on

        for jj = 1:numel(rfCents)

            p2(jj) = plot(retDisps{ii,jj},exGabor,'color',colors(jj,:),'linewidth',2);

        end

        set(gca,'plotboxaspectratio',[1 1 1],'xlim',[-2 2],'ytick',[],'fontsize',20);
        xlabel('Retinal disparity');
        ylabel('Cell response (sp/s)');
        title(titles{ii});
        legend(p2,{'0','+/-10','+/-20'},'location','northwest');

    end

end