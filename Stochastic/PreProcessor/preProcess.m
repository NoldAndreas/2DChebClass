function optsStruct=preProcess(optsStruct)
% optsStruct=preProcess(optsStruct)
%   Determines which calculations we're doing,
%   adds default parameters for any missing,
%   reformats some parameters,
%   determines what plotting we're doing
%
% INPUTS:
%   optsStruct  -- struct of size [1 1] containing
%
%                  PHYSICAL PARAMETERS:
%
%                   nParticlesS: number of particles (nSpecies,1)
%                   kBT:         temperature (scalar)
%                   D0S:         diffusion coefficient (nSpecies,1) 
%                   gammaS:      friction coefficient (nSpecies,1) 
%                   mS:          particle masses (nSpecies,1) 
%                   geom:        geometry (string)
%                   stocDim:     stochastic (i.e. full) dimension (integer)
%                   DDFTDim:     DDFT (symmetry adapted) dimension (integer)
%                   tMax:        end time (scalar)
%
%                  POTENTIALS:
%
%                   V1DV1:           name of 1-body potential ({string} or
%                                     {string a,string b,string c})
%                   potParamsNames:  names of V1 parameters {{'p1',...,'pn'}}
%                                     or {{'p1',...,'pn'},{...},{...}}
%                   <V1 Parameters>: called p1S, ..., pnS of form
%                                     {(nSpecies,1)} or 
%                                     {(nSpecies,1),(nSpecies,1),(nSpecies,1)}
%
%                   V2DV2:           name of 2-body potential (string)
%                   potParams2Names: names of V2 parameters {{'p1',...,'pn'}}
%                   <V2 Parameters>: called p1S, ..., pnS of form (nSpecies,nSpecies)
%
%                   HIParamsNames:   names of HI parameters {{'p1',...,'pn'}}
%                   <HI Parameters>: called p1S, ..., pnS of form (nSpecies,nSpecies)
%
%                  STOCHASTIC OPTIONS (optional):
%
%                   stocName:   names of stoc calculations 
%                               ({string 1, ..., string nStoc}}
%                   nRuns:      number of runs (scalar)
%                   nSamples:   number of initial and final samples (scalar)   
%                   stocType:   type of stoc calculation -- 'r' overdamped,
%                               'rv' with inertia ({type 1, ... type nStoc})
%                   stocHI:     whether stoc calculations include HI
%                               ({t/f 1, ..., t/f nStoc})
%                   stocHIType: name of HI interaction ({HI 1, ..., HI nStoc})
%                   tSteps:     number of time steps for each calculation
%                               {N 1, ..., N nStoc}
%
%                  DDFT OPTIONS (optional):
%
%                   DDFTName:       names of DDFT calculations 
%                                    ({string 1, ..., string nDDFT}}
%                   DDFTDir:        directory to add to path (string or 
%                                    {string 1, ..., string nDDFT}
%                   DDFTCode:       names of DDFT code (string or 
%                                    {string 1, ..., string nDDFT}
%                   DDFTParamsNames: names of DDFT parameters
%                                    ({{'p1', ..., 'pn'}} or
%                                     {{'p1', ..., 'pn'}, ..., {}) of
%                                     length nDDFT)
%                   DDFTParamsSaveNames: subset of DDFTParamsNames to use 
%                                        in constructing file names, same
%                                        format
%                   <DDFT Parameters>: called p1, ..., pn
%
% OUTPUTS:
%   optsStruct  -- struct of size [1 1] containing all necessary parameters
%                   for the calculations

%--------------------------------------------------------------------------
% Determine which calculations we're doing
%--------------------------------------------------------------------------

if(isfield(optsStruct,'doStoc'))
    optsStruct.anyStoc=any( cell2mat(optsStruct.doStoc) );
    optsStruct.nStoc=length(optsStruct.stocName);
else
    optsStruct.anyStoc=false;
end

if(isfield(optsStruct,'doDDFT'))
    optsStruct.anyDDFT=any( cell2mat(optsStruct.doDDFT) );
    optsStruct.nDDFT=length(optsStruct.DDFTName);
else
    optsStruct.anyDDFT=false;
end

if (~optsStruct.anyStoc && ~optsStruct.anyDDFT)
    return;
end

%--------------------------------------------------------------------------
% Assign defaults to any missing values
%--------------------------------------------------------------------------

optsStruct=setDefaults(optsStruct);

%--------------------------------------------------------------------------
% Calculate constants
%--------------------------------------------------------------------------

optsStruct.nParticles=sum(optsStruct.nParticlesS);
optsStruct.nSpecies=length(optsStruct.nParticlesS);

%--------------------------------------------------------------------------
% Duplicate required physical options for stochastic calculations
%--------------------------------------------------------------------------

optsStruct=repPhysVars(optsStruct);

%--------------------------------------------------------------------------
% Set up V1 parameters
%--------------------------------------------------------------------------

optsStruct=setV1Params(optsStruct);

%--------------------------------------------------------------------------
% Set up V2 parameters
%--------------------------------------------------------------------------

optsStruct=setV2Params(optsStruct);

%--------------------------------------------------------------------------
% Set up HI parameters
%--------------------------------------------------------------------------

optsStruct=setHIParams(optsStruct);

%--------------------------------------------------------------------------
% Set up DDFT parameters
%--------------------------------------------------------------------------

optsStruct=setDDFTParams(optsStruct);

%--------------------------------------------------------------------------
% Determine which plots we're doing
%--------------------------------------------------------------------------

optsStruct.anyPlots=any([optsStruct.doPdfs,optsStruct.doMovieGif,optsStruct.doMovieSwf,optsStruct.doMovieAvi, ...
              optsStruct.doInitialFinal,optsStruct.doMeans,optsStruct.doCustom,optsStruct.doEquilibria]);
optsStruct.anyPlotsP=any([optsStruct.doPdfsP,optsStruct.doMovieGifP,optsStruct.doMovieSwfP, ...
               optsStruct.doInitialFinalP,optsStruct.doCustomP]);

end