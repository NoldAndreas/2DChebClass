function optsStoc=getOptsStoc(optsStruct)
% optsStoc=getOptsStoc(optsStruct,fileStruct)
%  Collates stochastic options into the output structure
%
% INPUTS:
%   optsStruct  -- struct of size [1 1] containing
%                   tMax:          end time (scalar)
%                   tSteps:        number if time steps (nStoc,1)
%                   plotTimes:     integer number of save times
%                   nRuns:         integer number of runs
%                   initialGuess:  initial guess file (string)
%                   loadSamples:   t/f to load samples
%                   fixedInitial:  t/f for fixed initial condition
%                   sampleFinal:   t/f to sample final distribution
%                   nSamples:      integer number of samples
%                   sampleShift:   number of samples to initially skip over
%                   sampleSkip:    number of samples to skip each time
%                   thin:          see slicesample
%                   burnin:        see slicesample
%                   poolsize:      number of processors for parallel computing
%                   doStoc:        t/f do stoc calculation {nStoc,1}
%                   loadStoc:      t/f load stoc calculation {nStoc,1}
%                   saveStoc:      t/f save stoc calculation {nStoc,1}
%                   stocHI:        t/f whether stoc calculation has HI {nStoc,1}
%                   stocHIType:    string for HI type {'HI1', ..., 'HInStoc')
%                   stocType:      'r'/'rv' without/with inertia {'x1',...,'xnStoc')
%                   stocName:      name of stoc calculation {'name1',...,'namenStoc')
%                   noise:         t/f whether to include noise
%                   MBp:           t/f to sample momentum from Maxwell-Boltzmann.
%       
%   fileStruct  -- struct of size [1 1] containing
%                   saveFileI:     string save file for initial eq
%                   saveFileF:     string save file for final eq
%                   saveFilePIF:   string save file for momentum eq
%                   saveFileStoc:  string save file for stoc data
%                   saveFile:      empty placeholder
%
% OUTPUTS:
%   optsStoc  -- struct of size [1 1] containing all the above data

if( optsStruct.anyStoc ) 
    optsStoc=struct('tMax',optsStruct.tMax, ...
                    'tSteps',optsStruct.tSteps, ...
                    'plotTimes',optsStruct.plotTimes, ...
                    'nRuns',optsStruct.nRuns, ...
                    'initialGuess',optsStruct.initialGuess, ...
                    'loadSamples',optsStruct.loadSamples, ...
                    'fixedInitial',optsStruct.fixedInitial, ...
                    'sampleFinal',optsStruct.sampleFinal,...
                    'nSamples',optsStruct.nSamples, ...
                    'sampleShift',optsStruct.sampleShift, ...
                    'sampleSkip',optsStruct.sampleSkip, ...
                    'thin',optsStruct.thin, ...
                    'burnin',optsStruct.burnin, ...
                    'poolsize',optsStruct.poolsize, ...
                    'doStoc',optsStruct.doStoc,...
                    'loadStoc',optsStruct.loadStoc, ...
                    'saveStoc',optsStruct.saveStoc, ...
                    'HI',optsStruct.stocHI, ...    
                    'HIType',optsStruct.stocHIType, ...    
                    'useDivergence',optsStruct.stocUseDivergence, ...    
                    'type',optsStruct.stocType, ...
                    'stocName',optsStruct.stocName, ...
                    'noise',optsStruct.noise, ...
                    'useNewHS',optsStruct.useNewHS, ...
                    'MBp',optsStruct.MBp);
%                     'saveFileI',fileStruct.saveFileI, ...
%                     'saveFileF', fileStruct.saveFileF, ...
%                     'saveFilePIF',fileStruct.saveFilePIF, ...
%                     'saveFile',fileStruct.saveFile, ...
%                     'saveFileStoc',fileStruct.saveFileStoc,...
%                     'saveFileStocMeans',fileStruct.saveFileStocMeans);

    % cut down to the calculations we want to do
    optsStoc=optsStoc(cell2mat(optsStruct.doStoc));
else
    optsStoc=[];
end

end