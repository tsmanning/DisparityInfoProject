clear all; close all;

% load combo data file
load('./dataV1V2Combo/AllDTData.mat');
V1V2 = DTdata;

% counters for V1/V2 units
V1cnt = 0;
V2cnt = 0;

% plot RF locations and V1/V2 identity
figure; hold on;

for x = 1:length(V1V2)
    
    % get RF coordinates
    rf(x,:) = V1V2(x).rf;
    
    if V1V2(x).visualarea == 1
        h(1) = scatter(rf(x,1),rf(x,2),20,ColorIt('r'),'filled');
        V1cnt = V1cnt + 1;
    elseif V1V2(x).visualarea == 2
        h(2) = scatter(rf(x,1),rf(x,2),20,ColorIt('b'),'filled');
        V2cnt = V2cnt + 1;
    end
    
end

axis equal square; box on;
refline(0,0);
xline(0);
xlabel('azimuth?'); ylabel('elevation?');
legend(h,'V1','V2');
saveas(gcf,['./plots/V1V2_VFs.png']); 

display(['V1:' num2str(V1cnt) '   V2:' num2str(V2cnt)]);
