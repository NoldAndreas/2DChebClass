function [plotPosMask,plotTimes]=getPlotPos(optsStoc)
%[plotPosMask,plotTimes]=getPlotPos(optsStoc)
%   calculates approximations to the requested plotTimes which are integer
%   multiples of dt and also returns a mask for their positions in the full
%   time vector (0:dt:tMax) used by stochastic calculations
%
% INPUTS:
%   optsStoc    -- struct of size [1 1] containing
%                  plotTimes     (vector of desired plot times)
%                  tMax          (maximum time)
%                  tSteps        (number of time steps in stoc calculation)
%
% OUTPUTS:
%   plotMask    -- logical vector giving the locations of the approximate
%                   plot times in the full time vector
%   plotTimes   -- vector of approximate plot times


% get times
plotTimes=optsStoc.plotTimes;
tMax=optsStoc.tMax;
tSteps=optsStoc.tSteps;
dt=tMax/tSteps;
times=(0:tSteps)*dt;

% find the place in times closest to each of plotTimes
plotPosMask(floor(plotTimes/dt)+1)=1;

% change to logical so we can use it as a mask
plotPosMask=logical(plotPosMask);

% determine the plotting times, which should be close to the requested
% plotting times as long as dt is small
plotTimes=times(plotPosMask);

