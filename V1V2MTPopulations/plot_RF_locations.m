function plot_RF_locations(V1,V2,MT,flag,topDir)

%% Plot receptive field locations from the current sample of cells from each area

saveOn = 1;

f1 = figure; 
f1.Position = [100 100 650 600];
hold on; 
title('Receptive field locations');
scatter(V1.x_pos,V1.y_pos,[],ColorIt('b'),'filled');
scatter(V2.x_pos,V2.y_pos,[],ColorIt('g'),'filled');
scatter(MT.x_pos,MT.y_pos,[],ColorIt('r'),'filled');
plot([-12 12],[0 0],'k:');
plot([0 0],[-12 12],'k:');
axis([-12 12 -12 12]); axis equal tight;
legend('V1','V2','MT');
xlabel('horizontal eccentricity (deg)'); 
ylabel('vertical eccentricity (deg)'); 
box on;
set(gca,'fontsize',20);

% Save images
if saveOn

    saveas(f1,[topDir,'plots/RFLocations/RFlocations' flag '.svg']);

end

end
