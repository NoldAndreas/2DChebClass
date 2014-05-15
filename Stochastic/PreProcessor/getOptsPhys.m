function [optsPhys,optsPhysDDFT]=getOptsPhys(optsStruct)
% [optsPhysGIF,optsPhysDDFTGIF]=getOptsPhys(optsStruct)
%   Get structures containing physical parameters and options for both
%   stochastic and DDFT calculations
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
%                   tMax:        end time (scalar)
%
%                  POTENTIALS:
%
%                   V1DV1:           name of 1-body potential ({string} or
%                                     {string a,string b,string c})
%                   potParamsNames:  names of V1 parameters of form
%                                     {{'p1',...,'pn'},{...},{...}}
%                   V2DV2:           name of 2-body potential (string)
%                   potParams2Names: names of V2 parameters {{'p1',...,'pn'}}
%                   potNames       : save names of potentials of form
%                                     {nameG, nameI, nameF}
%
%                  STOCHASTIC OPTIONS (optional):
%
%                   anyStoc:    doing any stoc calculations (t/f)
%                     IF SO
%                       D0:          diffusion coefficient (nParticles*stocDim,1) 
%                       gamma:       friction coefficient (nParticles*stocDim,1)
%                       m:           particle masses (nParticles*stocDim,1) 
%                       stocDim:     stochastic (i.e. full) dimension (integer)                 
%                       potParams:   structure containing p1S, ..., pnS of form
%                                     {(nSpecies,1),(nSpecies,1),(nSpecies,1)}
%                                     and p1, ..., pn of form
%                                     {pG,pI,pF} with pX (nParticles,1)
%                       potParams2:  structure containing p1S, ..., pnS of form 
%                                     (nSpecies,nSpecies)
%                                     and p1, ..., pn of form (nParticles,nParticles)
%                       HIParams:    structure containing p1S, ..., pnS of form 
%                                     (nSpecies,1)
%                                     and p1, ..., pn of form (nParticles*stocDim,1)                                        
%
%                  DDFT OPTIONS (optional):
%
%                   anyDDFT:    doing any DDFT calculations (t/f)
%                     IF SO
%                       DDFTDim:     DDFT (symmetry adapted) dimension (integer)
%                       potParamsDDFT:   structure containing p1S, ..., pnS of form
%                                         {(nSpecies,1),(nSpecies,1),(nSpecies,1)}
%                                         and p1, ..., pn = p1S, ..., pnS
%                       potParams2DDFT:  structure containing p1S, ..., pnS of form 
%                                          (nSpecies,nSpecies)
%                                          and p1, ..., pn = p1S, ..., pnS
%                       HIParamsDDFT:    structure containing p1S, ..., pnS of form 
%                                          (nSpecies,1)
%                                          and p1, ..., pn = p1S, ..., pnS
%
% OUTPUTS:
%   Structures containing relevant input data:
%
%   optsPhysGIF      -- struct of size [1 3] containing all physica parameters
%                        for stochastic calculations (or DDFT if ~anyStoc)
%   optsPhysGIFDDFT  -- struct of size [nDDFT 3] containing all physical parameters
%                        for DDFT calculations


%--------------------------------------------------------------------------
% Set global options common to both stochastic and DDFT
%--------------------------------------------------------------------------

optsPhysGlobal=struct('nParticles',optsStruct.nParticles, ...
                'nParticlesS',optsStruct.nParticlesS, ... 
                'kBT',optsStruct.kBT, ...
                'mS',optsStruct.mS, ...
                'tMax',optsStruct.tMax, ...
                'D0S',optsStruct.D0S, ...
                'gammaS',optsStruct.gammaS, ...
                'geom',optsStruct.geom, ...
                'potNames',{optsStruct.potNames});

%--------------------------------------------------------------------------
% Set stochastic physical options
%--------------------------------------------------------------------------
            
if(optsStruct.anyStoc)

    optsPhys=struct('m',optsStruct.m, ...
                'D0',optsStruct.D0, ...
                'gamma',optsStruct.gamma, ...
                'dim',optsStruct.stocDim, ...
                'V1DV1',optsStruct.V1DV1, ...
                'V2DV2',optsStruct.V2DV2, ...
                'type','stoc');
 
    % merge in global parameters and parameters for potentials
    optsPhys=mergeStruct(optsPhysGlobal, ...
                            optsPhys, ...
                            optsStruct.potParams, ...
                            optsStruct.potParams2, ...
                            optsStruct.HIParams); 
end

%--------------------------------------------------------------------------
% Set DDFT plotting options
%--------------------------------------------------------------------------

if(optsStruct.anyDDFT)
    
    nDDFT=optsStruct.nDDFT;
    
    % note we use the species data for standard physical parameters                             
    optsPhysDDFTGlobal=struct('dim',optsStruct.DDFTDim, ...
                                    'type','DDFT');  

    optsPhysDDFTGlobal = mergeStruct(optsPhysGlobal, ...  
                                    optsPhysDDFTGlobal );
                                   
    for iDDFT=1:nDDFT
        optsPhysDDFT(iDDFT) = optsPhysDDFTGlobal;  %#ok

    end
                       
    for iDDFT=1:nDDFT
        optsPhysDDFT(iDDFT).V1 = optsStruct.potParamsDDFT{iDDFT}; %#ok
        optsPhysDDFT(iDDFT).V2 = optsStruct.potParams2DDFT; %#ok
        optsPhysDDFT(iDDFT).HI = optsStruct.HIParamsDDFT; %#ok
    end
    
    if(~optsStruct.anyStoc)
        optsPhys=optsPhysDDFT(1);
    end
    
else
    
    optsPhysDDFT=[];
    
end

end