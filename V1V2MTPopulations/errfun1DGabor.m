% GABOR MODEL
% r0 + a*exp(-(d-d0)^2 / (2*s^2)) * cos(2*pi*f*(d-d0)+phi)
% a <= 1.5*(max_neural_response - min_neural_response)
% f > 0
% half-way rectify

function[err] = errfun1DGabor( p, dat )
    
    % evaluate 1D Gabor
    x = dat(:,1);
    g = p(1) + p(2)*exp( -(x-p(3)).^2 / (2*p(4)^2) ) .* cos(2*pi*p(5)*(x-p(3))+p(6));
    
    % error
    err = mean( (dat(:,2) - g).^2 ); % following DeAngelis and Uka
    
    % evaluate over a finely sampled lattice
    x = [-2:0.05:2]; 
    g = p(1) + p(2)*exp( -(x-p(3)).^2 / (2*p(4)^2) ) .* cos(2*pi*p(5)*(x-p(3))+p(6)); % re-evaluate over [-2,2] horizontal disparity
    
    % gently penalize fits for which the max sps is less than the sum of
    % the baseline and the amplitude, because this means that the cos
    % component has been set to a very low freq even though the neuron has
    % a higher response rate
    %err = err + 0.001*( ((p(1)+p(2)) / max(g)) - 1 );

    % penalize error if response goes below a baseline
    baseline = 0.05;
    if( any( g < baseline ) )
        err = 10000000*err;
    end

%     % for V1/V2 refits, don't allow max spike rate to exceed max observed
%     % rate by more than a factor of 2
%     max_fac = 2; %lower min for V1/V2 fits with low spike rates
%     if( any( g > max_fac*(max(dat(:,2))) ) )
%         err = 10000000*err;
%     end
    
    % penalize error if this is a degenerate gabor with very low SF + high amplitude
    if p(5) <= .25
        err = 10000000*err;
    end