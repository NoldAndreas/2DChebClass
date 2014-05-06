function x=samplepdf(optsPhys,optsStoc)
% x=samplepdf(optsPhys,optsStoc)
%   performs slice sampling of equilibrium position distribution for a 
%   given potential using slicesample.m
%
% INPUTS: 
%  optsPhys  -- a structure of size [1 1] containing at least:
%         geom          (geometry; 'spherical'/'planar')
%         dim           (dimension; 1 or 3)
%         kBT           (temperature)
%         nParticles    (number of particles)
%         sigma         (particle diameter)
%         m             (particle mass)
%         V1DV1         (string identifying potential and gradient
%                         function)
%         t             (time at which to evaluate potential)
%         [pot params]  (parameters for potential)
%
% optsStocVR -- a structure of size [1 1] containing at least
%         nSamples      (number of samples for intial/final distribution)
%         thin          (option for slicesampling, integer, default 1)
%         burnin        (            "           , integer, default 0)
%     
% OUTPUTS:
%  x           -- a matrix of positions of size 
%                 (nSamples, dim x nParticles)


% get temperature which scales the probabilities
kBT=optsPhys.kBT;
% and time at which to evaluate the potential
t=optsPhys.t;

% and options for sampling
nSamples=optsStoc.nSamples;
thin=optsStoc.thin;
burnin=optsStoc.burnin;
   
% choose the initial data on a grid, this doesn't matter except that we
% need initial data where the particles don't overlap
getInitialGuess=str2func(optsStoc.initialGuess);
initialGuess=getInitialGuess(optsPhys);

%temp=getV(initialGuess,t,optsPhys);

% slice sample takes the matrix the other way around
initialGuess=initialGuess';

tic
% myslicesample is just slicesample with a progress bar
x = myslicesample(initialGuess,nSamples,'logpdf',@logpdf,'thin',thin,'burnin',burnin);  
toc


%--------------------------------------------------------------------------
% logpdf function generated from V
%--------------------------------------------------------------------------
function logp = logpdf(x)
    % note that this needs to be here in order to inherit the parameters

    % get potential (where we don't pass a time argument) and calculate the
    % log of the probability, which should be a single number
    V=getV(x,t,optsPhys);
    logp=-V/kBT;

end % logpdf

end % end samplepdf