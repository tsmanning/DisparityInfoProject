clear all; close all;

% sign convenrtion: leftward in VF = negative, crossed disparity = negative

% screen distance in Greg's study (m)
screen_dist = 0.57;

% Rough average in macaque (in m)
ipd = 0.03;
xL  = -ipd/2;
xR  = ipd/2;

% vergence angle associated with center of screen
x_cntr = 0;
z_cntr = screen_dist;

theta_L = atand((x_cntr - xL)/z_cntr);
theta_R = atand((x_cntr - xR)/z_cntr);

center_verg = theta_R - theta_L;

% calculate disparities for a range of eccentricities (deg)
ecc = linspace(-20,20,101);

for ii = 1:length(ecc)

    this_ecc = ecc(ii);

    % vergence angle associated with a point at this eccentricity
    x_side = screen_dist*tand(this_ecc);

    theta_L = atand((x_side - xL)/z_cntr);
    theta_R = atand((x_side - xR)/z_cntr);

    side_verg = theta_R - theta_L;

    % disparity
    screen_disparity(ii) = side_verg - center_verg;

    % repeat for a set of non-zero disparities (deg)
    % negative is crossed, positive is uncrossed
    disp = [-2 -1 0 1 2];
    %disp = -2;

    for jj = 1:length(disp)

        this_disp = disp(jj);

        % on screen shift
        screen_shift = abs(screen_dist*tand(this_disp/2));

        % make it crossed or uncrossed
        if sign(this_disp) == -1
            % crossed
            delta_xL = screen_shift;
            delta_xR = -screen_shift;
        else
            % uncrossed
            delta_xL = -screen_shift;
            delta_xR = +screen_shift;
        end

        % vergence angle at center of screen
        theta_L = atand((x_cntr + delta_xL - xL)/z_cntr);
        theta_R = atand((x_cntr + delta_xR - xR)/z_cntr);

        center_verg_disp = theta_R - theta_L;

        % vergence angle at side of screen
        theta_L = atand((x_side + delta_xL - xL)/z_cntr);
        theta_R = atand((x_side + delta_xR - xR)/z_cntr);

        side_verg = theta_R - theta_L;

        % disparity
        other_disparity(jj,ii) = side_verg - center_verg_disp;

    end

end

figure; hold on;
plot(ecc,screen_disparity)
plot(ecc,other_disparity)

%% Plot
f4 = figure;
f4.Position = [700 700 650 600];
hold on;

colors = [flipud(linspace(0,1,numel(disp))') zeros(numel(disp),1) linspace(0,1,numel(disp))'];

% plot(angs,a,'--k','linewidth',2);
for ii = 1:numel(disp)
    plot(ecc,other_disparity(ii,:),'color',colors(ii,:),'linewidth',2);
    text(-3,0.5-0.035*(ii-1),[num2str(disp(ii)),'\circ disp'],'fontsize',20,'color',colors(ii,:));
end
xlabel('Eccentricity');
ylabel('Disparity error (Screen - Retinal)');
set(gca,'fontsize',20);

saveas(gcf,'disparityError_Emily.png')
saveas(gcf,'disparityError_Emily.eps','epsc')