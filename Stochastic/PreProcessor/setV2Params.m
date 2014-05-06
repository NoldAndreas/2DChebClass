function optsStruct=setV2Params(optsStruct)
% optsStruct=setV2Params(optsStruct)
%   Replicate V2 parameter variables for stochastic calculations so
%   that they're the correct dimension (nParticles,nParticles) 
%   and DDFT calculations so they agree with the species-wise definitions
%
% INPUTS:
%   optsStruct  -- struct of size [1 1] containing
%                   potParams2Names: names of V2 parameters {{'p1',...,'pn'}}
%                   <V2 Parameters>: called p1S, ..., pnS of form (nSpecies,nSpecies)
%                   anyStoc:         are we doing any stoc calculations (t/f)
%                    IF SO
%                      nParticles:   total number of particles (integer)
%                      nParticlesS:  number of particles (nSpecies,1)
%                      nSpecies:     number of species (integer)
%                   anyDDFT:         are we doing any DDFT calculations (t/f)
%
% OUTPUTS:
%   optsStruct  -- struct of size [1 1] with all fields of input plus
%                     IF anyStoc
%                       <V2 parameters>:  called p1, ..., pn, of form 
%                                          (nParticles,nParticles)
%                       potParams2:       structure containing p1S, ..., pnS
%                                          and p1, ..., pn
%                     IF anyDDFT
%                       potParams2DDFT:   structure containing p1S, ..., pnS
%                                          and p1, ..., pn, = p1s, ..., pnS

%--------------------------------------------------------------------------
% Duplicate V2 parameters for stochastic calculation
%--------------------------------------------------------------------------

if(optsStruct.anyStoc)

    nParticles=optsStruct.nParticles;   %#ok
    nParticlesS=optsStruct.nParticlesS;
    nSpecies=optsStruct.nSpecies;
    
    endPos=cumsum(nParticlesS);
    startPos= [1;endPos(1:end-1)+1];

    vars=optsStruct.potParams2Names;
    
    nVars=length(vars);
    
    for iVar=1:nVars
        % set up empty variable
        eval(['optsStruct.' vars{iVar} '=zeros(nParticles,nParticles);'])    
        
        % assign variable of size (nParticles,nParticles) with the correct 
        % value for each pair of species
        for iSpecies=1:nSpecies
            rowMask = startPos(iSpecies):endPos(iSpecies);      %#ok
            for jSpecies=1:nSpecies
                colMask = startPos(jSpecies):endPos(jSpecies);  %#ok      
                eval(['optsStruct.' vars{iVar} '(rowMask,colMask) = optsStruct.' vars{iVar} 'S(iSpecies,jSpecies);'])
                
            end 
        end
        
    end
    
end
    
%--------------------------------------------------------------------------
% Set up stochastic V2 parameters in potParams2
%--------------------------------------------------------------------------

% remove extra braces
vars=optsStruct.potParams2Names;

nVars=length(vars);

if(optsStruct.anyStoc)
    addVarText=[];

    % append each variable and its value to the command
    % var is (nParticles,nParticles) as above, varS as in input
    
    if(iVar==1)
        comma = [];
    else
        comma = ',';
    end
    
    for iVar=1:nVars
        % add var with value var
        addVarText=cat(2,addVarText, [comma '''' vars{iVar} ''',optsStruct.'  vars{iVar}]);
        % add varS with value varS
        addVarText=cat(2,addVarText, [',''' vars{iVar} 'S'',optsStruct.' vars{iVar} 'S']);
    end

    potParamsCmd=['optsStruct.potParams2 = struct(' addVarText ');'];

    eval(potParamsCmd);
  
end
    
%--------------------------------------------------------------------------
% Set up DDFT V2 parameters in potParams2DDFT
%--------------------------------------------------------------------------

if(optsStruct.anyDDFT)
    
    addVarText=['''V2DV2'', ''' optsStruct.V2DV2 '''' ];

    % potParams2DDFT has var and varS with the same values as varS in input
    for iVar=1:nVars
        % add var with value varS
        addVarText=cat(2,addVarText, [',''' vars{iVar} ''',optsStruct.'  vars{iVar} 'S']);
        % add varS with value varS
        addVarText=cat(2,addVarText, [',''' vars{iVar} 'S'',optsStruct.' vars{iVar} 'S']);
    end
   
    potParamsCmd=['optsStruct.potParams2DDFT = struct(' addVarText ');'];

    eval(potParamsCmd);
    
end

end

