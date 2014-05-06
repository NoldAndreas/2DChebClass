function rmPaths()
% rmPaths()
%   Removes the paths added by addPaths
%
% INPUTS:
%   <none>
%
% OUTPUTS:
%   <none>

rmpath('Inputs');
rmpath('Potentials');
rmpath('PreProcessor');
rmpath('Stochastic');
rmpath(['Stochastic' filesep 'InitialGuess'])
rmpath(['Stochastic' filesep 'HI']) 
rmpath('Plotting');

end