% Run Fano factor check to ID how far each cortical area departs from the Poisson assumptions under 
% our FI calculation and make sure there aren't substantial difference from
% each other

% flag indicating whether to correct the stimulus disparities to account
% for planar screen (horopter deviation and foreshortening)
correct_screen_disparity = 1;

% flag whether to include all data or only data included in paper
subsample = 1;

FanoFacCheck(correct_screen_disparity,subsample);