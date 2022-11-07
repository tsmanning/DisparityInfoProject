function[FI, P, S, X, E] = fit1DGabor( tag, dat, FI, P, S, X, E, n, N )

% number of iterations for fitting
niters = 200;

% convert disparity and mean responses to row vectors
x_orig  = dat(:,1)';
resp_orig = dat(:,2)'; %

% interpolate disparities and responses to constrain fits
xi = zeros(1,2*length(x_orig)-1);
xi(1:2:end) = x_orig;
for i = 2 : 2 : length(xi)-1
    xi(i) = (xi(i-1)+xi(i+1))/2.0; % interpolate between each measured disparity
end

ri = interp1(x_orig,resp_orig,xi,'linear');

ri = abs(ri); % all responses should be positive


% fit 1D Gabor to 1D data [Gabor: r0 + a*exp(-(d-d0)^2 / (2*s^2)) * cos(2*pi*f*(d-d0)+phi)]
% parameters: [pedestal amplitude Gaussian_mean Gaussian_sigma frequency phase]

% initialize optimization
options = optimset('Display','off');
minp = NaN;
minerr = inf;

% % start by constraining the frequency of the cos to zero (fitting with a Gaussian, set phase = 0 also)
% 
% lb = [0   0      -1.75 0    0 0]; % lower bounds
% ub = [500 500    1.75  5    0 0]; % upper bounds
% 
% for k = 1 : niters
% 
%     % random seed
%     seed = [uniform_sample_in_range(lb(1),ub(1)) ...
%         uniform_sample_in_range(lb(2),ub(2)) ...
%         uniform_sample_in_range(lb(3),ub(3)) ...
%         uniform_sample_in_range(lb(4),ub(4)) ...
%         uniform_sample_in_range(lb(5),ub(5)) ...
%         uniform_sample_in_range(lb(6),ub(6))];
% 
%     p = fmincon( 'errfun1DGabor', seed, [], [], [], [], lb, ub, [], options, [xi' ri'] );
%     err = errfun1DGabor( p, [xi' ri'] );
% 
%     if( k == 1 || err < minerr )
% 
%         fprintf( '%d, Gaussian iteration %d: found better fit: %.2f %.2f %.2f %.2f %.2f %.2f (%f) --> %.2f %.2f %.2f %.2f %.2f %.2f (%f)\n', ...
%             n, k, minp, minerr, p, err );
% 
%         minp = p;
%         minerr = err;
% 
%     end
% 
% end
% 
% % calculate R2 of best Gaussian parameters
% 
% % calculate fitted spike rates for each tested disparity
% TF = minp(1) + minp(2)*exp( -(x_orig-minp(3)).^2 / (2*minp(4)^2) ) .* cos(2*pi*minp(5)*(x_orig-minp(3))+minp(6));
% 
% SSR = sum((resp_orig - TF).^2);
% TSS = sum((resp_orig - mean(resp_orig)).^2);
% r2   = 1 - (SSR/TSS);
% 
% fprintf('Gaussian R2 = %.2f \n', r2)

% if this isn't a great fit, try a Gabor
% parameters: [pedestal amplitude Gaussian_mean Gaussian_sigma frequency phase]

%lb = [0   0      -1.75 0 .1*pi -2*pi]; % lower bounds
%ub = [500 500    1.75  5 3   2*pi]; % upper bounds

lb = [0   0      -1.75 0 0 -2*pi]; % lower bounds
ub = [500 500    1.75  5 4.5   2*pi]; % upper bounds

% lb = [0   0      -1.75 0 0 -2*pi]; % lower bounds
% ub = [500 500    1.75  5 1   2*pi]; % upper bounds

%if r2 < 0.95

    for k = 1 : niters

        % random seed
        seed = [uniform_sample_in_range(lb(1),ub(1)) ...
            uniform_sample_in_range(lb(2),ub(2)) ...
            uniform_sample_in_range(lb(3),ub(3)) ...
            uniform_sample_in_range(lb(4),ub(4)) ...
            uniform_sample_in_range(lb(5),ub(5)) ...
            uniform_sample_in_range(lb(6),ub(6))];

        p = fmincon( 'errfun1DGabor', seed, [], [], [], [], lb, ub, [], options, [xi' ri'] );
        err = errfun1DGabor( p, [xi' ri'] );

        % if error is better than best Gaussian, go with it
        if( err < minerr )

            fprintf( '%d, Gabor iteration %d: found better fit: %.2f %.2f %.2f %.2f %.2f %.2f (%f) --> %.2f %.2f %.2f %.2f %.2f %.2f (%f)\n', ...
                n, k, minp, minerr, p, err );

            minp = p;
            minerr = err;

        end

    end

%end

% calculate R2 of best Gabor parameters

% calculate fitted spike rates for each tested disparity
TF = minp(1) + minp(2)*exp( -(x_orig-minp(3)).^2 / (2*minp(4)^2) ) .* cos(2*pi*minp(5)*(x_orig-minp(3))+minp(6));

SSR = sum((resp_orig - TF).^2);
TSS = sum((resp_orig - mean(resp_orig)).^2);
r2   = 1 - (SSR/TSS);

fprintf('Final R2 = %.2f \n', r2)


% store final parameters
p = minp;

% visualize data and model
rng = 2;
xg1 = [-rng : 0.01 : rng]; %fixed spatial range
g1 = p(1) + p(2)*exp( -(xg1-p(3)).^2 / (2*p(4)^2) ) .* cos(2*pi*p(5)*(xg1-p(3))+p(6));


figure(1);
subplot(2,3,1);
cla; plot(xg1,g1,'b'); hold on;  plot(xi,ri,'ro');  hold off;
xlabel('horizontal disparity');
ylabel('model response');
set( gca, 'Ylim', [0 1.2*max(ri)] );
grid on;
h = title( sprintf('neuron = %d | err = %.3f', n, minerr) );

% model derivative
syms f(y)
f(y) =  p(1) + p(2)*exp( -(y-p(3)).^2 / (2*p(4)^2) ) .* cos(2*pi*p(5)*(y-p(3))+p(6));
Df = diff(f,y); % symbolic derivative
dg1 = eval(Df(xg1)); % numeric derivative

% compute and visualize Fisher information
fi = dg1.^2 ./ g1; % Fisher
subplot(2,3,3);
plot(xg1,fi,'b');
xlabel( 'horizontal disparity' ); ylabel( 'Fisher information' );
grid on;
set( gca, 'Xlim', [-rng rng] );

% visualize model deriviative
subplot(2,3,2);
plot(xg1,dg1,'b'); %symbolic deriviatve
hold on;
plot(xg1(1:end-1),100*diff(g1),'c:');
hold off;
xlabel('horizontal disparity');
ylabel('derivative of model response');
set( gca, 'Xlim', [-rng rng] );
grid on;

% Store across all neurons
if(n == 1)
    FI = zeros( N, length(xg1) ); % Fisher: first neuron, initialize
    P  = zeros( N, 6 ); % model parameters
    S(length(ri)).s = []; % spike count
    X(length(xi)).x = []; % horizontal disparity
    E  = zeros( N, 1 );
end

FI(n,:) = fi;
P(n,:)  = p;
S(n).s  = ri;
X(n).x  = xi;
E(n)    = minerr;

% visualize
subplot(2,3,5); imagesc( FI(1:n,:) );
h = colorbar; h.Label.String = 'Fisher information';
xlabel( 'horizontal disparity' ); ylabel( 'neuron' );

subplot(2,3,6); plot( xg1, sum(FI(1:n,:)), 'b' );
xlabel( 'horizontal disparity' ); ylabel( 'Fisher information' );
set( gca, 'Ylim', [0 1.1*max(sum(FI(1:n,:)))] );
set( gca, 'Xlim', [-rng rng] );

drawnow;

% print figure
cmd = sprintf( 'print -dpng plots%s/fig%d.png', tag, n );
eval( cmd );
