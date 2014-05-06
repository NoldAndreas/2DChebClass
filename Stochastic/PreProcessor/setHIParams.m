function optsStruct=setHIParams(optsStruct)
% optsStruct=setHIParams(optsStruct)
%   Replicate HI parameter variables for stochastic calculations so
%   that they're the correct dimension (nParticles*stocDim,1) 
%   and DDFT calculations so they agree with the species-wise definitions
%
% INPUTS:
%   optsStruct  -- struct of size [1 1] containing
%                   HIParamsNames:   names of HI parameters {{'p1',...,'pn'}}
%                   <HI Parameters>: called p1S, ..., pnS of form (nSpecies,1)
%                   anyStoc:         are we doing any stoc calculations (t/f)
%                    IF SO
%                      nParticles:   total number of particles (integer)
%                      nParticlesS:  number of particles (nSpecies,1)
%                      nSpecies:     number of species (integer)
%                      stocDim:      dimension of stoc calculation (integer)
%                   anyDDFT:         are we doing any DDFT calculations (t/f)
%
% OUTPUTS:
%   optsStruct  -- struct of size [1 1] with all fields of input plus
%                     IF anyStoc
%                       <HI parameters>:  called p1, ..., pn, of form 
%                                         (nParticles*stocDim,1)
%                       HIParams:         structure containing p1S, ..., pnS
%                                          and p1, ..., pn
%                     IF anyDDFT
%                       HIParams:         structure containing p1S, ..., pnS
%                                          and p1, ..., pn, = p1s, ..., pnS


%--------------------------------------------------------------------------
% Duplicate HI parameters for stochastic calculation
%--------------------------------------------------------------------------

if(optsStruct.anyStoc)
    
    nParticlesS=optsStruct.nParticlesS;
    nParticles=sum(nParticlesS); %#ok
    nSpecies=length(nParticlesS);
    stocDim=optsStruct.stocDim;
    
    endPosDim=cumsum(nParticlesS)*stocDim;
    startPosDim= [1;endPosDim(1:end-1)+1];

    vars=optsStruct.HIParamsNames;
   
    nVars=length(vars);
    
    for iVar=1:nVars
        % set up empty variable
        eval(['optsStruct.' vars{iVar} '=zeros(nParticles*stocDim,1);']) 
    
        % assign variable of size (nParticles*stocDim,1) with correct value
        % for each species
        for iSpecies=1:nSpecies
            mask=startPosDim(iSpecies):endPosDim(iSpecies); %#ok
            eval(['optsStruct.' vars{iVar} '( mask ) = optsStruct.' vars{iVar} 'S(iSpecies);'])
        end
    end
  
end

%--------------------------------------------------------------------------
% Set up stochastic HI parameters in HIParams
%--------------------------------------------------------------------------

% remove extra braces
vars=optsStruct.HIParamsNames;

nVars=length(vars);

if(optsStruct.anyStoc)
    addVarText=[];

    % append each variable and its value to the command
    % var is (nParticles*stocDim,1) as above, varS as in input
    for iVar=1:nVars
        if(iVar==1)
            comma = [];
        else
            comma = ',';
        end
        
        % add var with value var
        addVarText=cat(2,addVarText, [comma '''' vars{iVar} ''',optsStruct.'  vars{iVar}]);
        % add varS with value varS
        addVarText=cat(2,addVarText, [',''' vars{iVar} 'S'',optsStruct.' vars{iVar} 'S']);
    end
    
    HIParamsCmd=['optsStruct.HIParams = struct(' addVarText ');'];

    eval(HIParamsCmd);

end
    
%--------------------------------------------------------------------------
% Set up DDFT HI parameters in HIParamsDDFT
%--------------------------------------------------------------------------

if(optsStruct.anyDDFT)
    addVarText=[];

    % HIParamsDDFT has var and varS with the same values as varS in input
    for iVar=1:nVars
        if(iVar==1)
            comma = [];
        else
            comma = ',';
        end
        
        % add var with value varS
        addVarText=cat(2,addVarText, [comma '''' vars{iVar} ''',optsStruct.'  vars{iVar} 'S']);
        % add varS with value varS
        addVarText=cat(2,addVarText, [',''' vars{iVar} 'S'',optsStruct.'  vars{iVar} 'S']);
    end
    
    HIParamsCmd=['optsStruct.HIParamsDDFT = struct(' addVarText ');'];

    eval(HIParamsCmd);
    
end

end