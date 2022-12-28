% Check deviation from disparity on retina vs disparity defined on screen 

% ultimately the error between assuming the horopter and a frontoparallel
% screen are coplanar; Greg used a +/-20deg screen at 57cm. 

% Half of vergeance angle when fixated on a central target (rotation state
% of fovea)
%
% thetaL = atan2d(xL,z);
% thetaR = atan2d(xR,z);

% Angle of incidence for ray from point (at xP) on frontoparallel screen
%
% phiL = atan2d(xP-xL,z);
% phiR = atan2d(-xP+xR,z);

% Putting it all together, get disparity of all points on screen along 0deg
% vertical line
% thisF = @(z,xL,xR,xP) - atan2d(xP-xL,z) - atan2d(xL,z) - atan2d(xR-xP,z) + atan2d(xR,z);

clear all; close all;

% Eye positions
posL = @(z,xL,xP) atan2d(xP-xL,z) - atan2d(xL,z);
posR = @(z,xR,xP) atan2d(xP-xR,z) - atan2d(xR,z);

% Delta screen pos for given disparity (disp/2)
dotSep = @(z,disps,xR) 2*(z*tand(disps/2 + atan2d(xR,z)) - xR);

% Position from eccentricity
scrPos = @(z,angs) z*tand(angs);

%%%%%%%%%%%%%%%%%%%%

screenD = 0.89;
% screenD = 0.57;
angs    = linspace(-20,20,500);
% disps   = linspace(-2,2,500);
% xP      = screenD*tand(angs);
z       = screenD;

% Rough average is macaques is about 30mm
ipd = 0.03;
xL  = -ipd/2;
xR  = ipd/2;

% dispErr = thisF(z,xL,xR,xP);

% close all

% a = posR(z,xR,xP) - posL(z,xL,xP) + 2*atan2d(xR,z) - 2*atan2d(xL,z);
% 
% f = figure;
% f.Position = [100 100 650 600];
% hold on
% 
% plot(angs,dispErr,'k','linewidth',2);
% % plot(angs,a,'r','linewidth',2);
% xlabel('Stimulus eccentricity');
% ylabel('Disparity Error');
% set(gca,'fontsize',20);

% f2 = figure;
% f2.Position = [700 100 650 600];
% hold on;
% 
% % plot(-angs,-posL(z,xL,xP)-angs,'r','linewidth',2);
% plot(angs,posL(z,xL,xP)-angs+2*atan2d(xL,z),'r','linewidth',2);
% plot(angs,posR(z,xR,xP)-angs+2*atan2d(xR,z),'b','linewidth',2);
% xlabel('Stimulus eccentricity');
% ylabel('Retinal position error');
% set(gca,'fontsize',20);

% f3 = figure;
% f3.Position = [100 700 650 600];
% hold on;
% 
% plot(disps,dotSep(z,disps,xR),'k','linewidth',2);
% xlabel('Disparity');
% ylabel('\Delta position');
% set(gca,'fontsize',20);

%%%%%%%%%%%%%%%%%%%
qDisp = [-2 -1 0 1 2];

for ii = 1:numel(qDisp)
    % Get dot positions on screen
    xPos = scrPos(z,angs) + dotSep(z,qDisp(ii),xR)/2*[-1;1];

    % Get eye positions
    eyeL = posL(z,xL,xPos(1,:));
    eyeR = posR(z,xR,xPos(2,:));

    % Get retinal disparities
    retDisps(ii,:) = eyeR - eyeL + 2*atan2d(xR,z) - 2*atan2d(xL,z);
end

colors = [flipud(linspace(0,1,numel(qDisp))') zeros(numel(qDisp),1) linspace(0,1,numel(qDisp))'];

%% Plot
f4 = figure;
f4.Position = [700 700 650 600];
hold on;

% plot(angs,a,'--k','linewidth',2);
for ii = 1:numel(qDisp)
    plot(angs,retDisps(ii,:)-qDisp(ii),'color',colors(ii,:),'linewidth',2);
    text(-3,0.5-0.035*(ii-1),[num2str(qDisp(ii)),'\circ disp'],'fontsize',20,'color',colors(ii,:));
end
xlabel('Eccentricity');
ylabel('Disparity error (Screen - Retinal)');
set(gca,'fontsize',20,'ylim',[-0.1 0.6]);
plot([-20 20],[0 0],'--k');

saveas(gcf,'disparityError_Tyler89cm.png')
% saveas(gcf,'disparityError_Tyler.eps','epsc')




