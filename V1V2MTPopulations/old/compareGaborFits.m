% two parameters that produce very similar fits, even though amplitude and
% freq are quite different

clear all; close all;
% offset, amplitude, gaussian mean, gaussian sigma, freq, phase
param_set_1 = [168.1104     1.9219e+04  -0.0911     0.5318      0.0028      1.5710];
param_set_2 = [168.2300     168.8200    -0.0760     0.7380      0.3160      1.6360];

% support
x = [-5 : 0.01 : 5]; 

% evaluate different parts
gauss_1 = param_set_1(1) + param_set_1(2)*exp( -(x-param_set_1(3)).^2 / (2*param_set_1(4)^2) );
gauss_2 = param_set_2(1) + param_set_2(2)*exp( -(x-param_set_2(3)).^2 / (2*param_set_2(4)^2) );

wave_1 = cos(2*pi*param_set_1(5)*(x-param_set_1(3))+param_set_1(6));
wave_2 = cos(2*pi*param_set_2(5)*(x-param_set_2(3))+param_set_2(6));

% evaluate full Gabor
full_1 = param_set_1(1) + param_set_1(2)*exp( -(x-param_set_1(3)).^2 / (2*param_set_1(4)^2) ) .* cos(2*pi*param_set_1(5)*(x-param_set_1(3))+param_set_1(6));
full_2 = param_set_2(1) + param_set_2(2)*exp( -(x-param_set_2(3)).^2 / (2*param_set_2(4)^2) ) .* cos(2*pi*param_set_2(5)*(x-param_set_2(3))+param_set_2(6));

% plot
figure; hold on;
subplot(1,3,1); hold on; title('gaussian envelope');
plot(x,gauss_1,'k-');
plot(x,gauss_2,'r-');
legend('set 1','set 2');

subplot(1,3,2); hold on; title('sine wave');
plot(x,wave_1,'k-');
plot(x,wave_2,'r-');

subplot(1,3,3); hold on; title('Gabor');
plot(x,full_1,'k-');
plot(x,full_2,'r-');