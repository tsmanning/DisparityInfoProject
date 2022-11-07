function plot_tuning_characteristics(V1,V2,MT)


%% PREFERRED DISPARITIES AND BANDWIDTH

% plot preferred disparity histogram and bandwidth as a function of
% disparity (as well as RF size for MT)
histo_lattice = linspace(-2,2,21);

% envelope mean
figure; hold on;

subplot(3,2,1); hold on; title('Envelope mean');
hh(1) = histogram(V1.P(:,3),histo_lattice,'normalization','probability');
xlabel('envelope mean (deg)'); ylabel('freq');

subplot(3,2,1); hold on;
hh(2) = histogram(V2.P(:,3),histo_lattice,'normalization','probability');

subplot(3,2,1); hold on;
hh(3) = histogram(MT.P(:,3),histo_lattice,'normalization','probability');

legend('V1','V2','MT');

subplot(3,2,2); hold on; title('Envelope mean');
hh(1) = histogram(V1.P(:,3),histo_lattice,'normalization','probability');
xlabel('envelope mean (deg)'); ylabel('freq');

subplot(3,2,2); hold on; 
hh(2) = histogram(V2.P(:,3),histo_lattice,'normalization','probability');

subplot(3,2,2); hold on;
hh(3) = histogram(MT.P(MT.overlapping_RF_inds,3),histo_lattice,'normalization','probability');

legend('V1','V2','MT-overlap');

% peak response
subplot(3,2,3); hold on; title('Peak response');

hh(1) = histogram(V1.pref_disp,histo_lattice,'normalization','probability');
xlabel('Peak response (deg)'); ylabel('freq');

subplot(3,2,3); hold on;
hh(2) = histogram(V2.pref_disp,histo_lattice,'normalization','probability');

subplot(3,2,3); hold on;
hh(3) = histogram(MT.pref_disp,histo_lattice,'normalization','probability');

legend('V1','V2','MT');

subplot(3,2,4); hold on;  title('Peak response');

hh(1) = histogram(V1.pref_disp,histo_lattice,'normalization','probability');
xlabel('Peak response (deg)'); ylabel('freq');

subplot(3,2,4); hold on;
hh(2) = histogram(V2.pref_disp,histo_lattice,'normalization','probability');

subplot(3,2,4); hold on;
hh(3) = histogram(MT.pref_disp(MT.overlapping_RF_inds),histo_lattice,'normalization','probability');

legend('V1','V2','MT-overlap');

% sigma
subplot(3,2,5); hold on; title('Envelope sigma');

hh(1) = histogram(V1.P(:,4),histo_lattice,'normalization','probability');
xlabel('Envelope sigma (deg)'); ylabel('freq');

subplot(3,2,5); hold on;
hh(2) = histogram(V2.P(:,4),histo_lattice,'normalization','probability');

subplot(3,2,5); hold on;
hh(3) = histogram(MT.P(:,4),histo_lattice,'normalization','probability');

legend('V1','V2','MT');

subplot(3,2,6); hold on; title('Peak response');

hh(1) = histogram(V1.P(:,4),histo_lattice,'normalization','probability');
xlabel('Peak response (deg)'); ylabel('freq');

subplot(3,2,6); hold on;
hh(2) = histogram(V2.P(:,4),histo_lattice,'normalization','probability');

subplot(3,2,6); hold on;
hh(3) = histogram(MT.P(MT.overlapping_RF_inds,4),histo_lattice,'normalization','probability');

legend('V1','V2','MT-overlap');

saveas(gcf,['./plots/TuningCharacteristics/TuningCharacteristics.png']); 