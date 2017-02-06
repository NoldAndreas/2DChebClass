function optsStruct=repPhysVars(optsStruct)
% optsStruct=repPhysVars(optsStruct)
%   Replicate compulsory physical variables for stochastic calculations so
%   that they're the correct dimension (nParticles*stocDim,1)
%
% INPUTS:
%   optsStruct  -- struct of size [1 1] containing
%                   D0S:         diffusion coefficient (nSpecies,1) 
%                   gammaS:      friction coefficient (nSpecies,1) 
%                   mS:          particle masses (nSpecies,1) 
%                   anyStoc:     are we doing any stoc calculations (t/f)
%                   nParticles:  total number of particles (integer)
%                   nParticlesS: number of particles (nSpecies,1)
%                   nSpecies:    number of species (integer)
%                   stocDim:     stochastic (i.e. full) dimension (integer)
%
% OUTPUTS:
%   optsStruct  -- struct of size [1 1] with all fields of input plus
%                   D0:         diffusion coefficient (nParticles*stocDim,1) 
%                   gamma:      friction coefficient (nParticles*stocDim,1)
%                   m:          particle masses (nParticles*stocDim,1) 



if(optsStruct.anyStoc)
    
    vars={'m','gamma','D0'};
    varsS={'mS','gammaS','D0S'};

    nParticles=optsStruct.nParticles;   %#ok
    nParticlesS=optsStruct.nParticlesS;
    nSpecies=optsStruct.nSpecies;
    stocDim=optsStruct.stocDim;
    
    endPosDim=cumsum(nParticlesS)*stocDim;
    startPosDim= [1;endPosDim(1:end-1)+1];
    
    for iVar=1:length(vars)
        eval(['optsStruct.' vars{iVar} '=zeros(nParticles*stocDim,1);'])       
        for iSpecies=1:nSpecies
            mask=startPosDim(iSpecies) : endPosDim(iSpecies); %#ok
            eval(['optsStruct.' vars{iVar} '( mask ) = optsStruct.' varsS{iVar} '(iSpecies);'])
        end
        
    end
    
end

end