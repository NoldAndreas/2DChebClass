function optsStruct=setICParams(optsStruct)
% optsStruct=setUParams(optsStruct)
%   Replicate U parameter variables for stochastic calculations so
%   that they're the correct dimension (nParticles,1) and DDFT calculations
%   so they agree with the species-wise definitions
%
% INPUTS:
%   optsStruct  -- struct of size [1 1] containing
%                   flowParamsNames:  names of U parameters {{'p1',...,'pn'}}
%                                     or {{'p1',...,'pn'},{...},{...}}
%                   <U Parameters>: called p1S, ..., pnS of form
%                                     {(nSpecies,1)} or 
%                                     {(nSpecies,1),(nSpecies,1),(nSpecies,1)}
%                   anyStoc:         are we doing any stoc calculations (t/f)
%                    IF SO
%                      nParticles:   total number of particles (integer)
%                      nParticlesS:  number of particles (nSpecies,1)
%                      nSpecies:     number of species (integer)
%                   anyDDFT:         are we doing any DDFT calculations (t/f)
%
% OUTPUTS:
%   optsStruct  -- struct of size [1 1] with all fields of input plus
%                     <U parameters>:    called p1S, ..., pnS of form
%                                          {(nSpecies,1),(nSpecies,1),(nSpecies,1)}
%                     IF anyStoc
%                       <U parameters>:  called p1, ..., pn, of form 
%                                          {pG,pI,pF} with pX (nParticles,1)
%                       potParams:        structure containing p1S, ..., pnS
%                                          and p1, ..., pn
%                     IF anyDDFT
%                       potParamsDDFT:    structure containing p1S, ..., pnS
%                                          and p1, ..., pn, = p1s, ..., pnS

%--------------------------------------------------------------------------
% Find list of IC parameters
%--------------------------------------------------------------------------

vars=optsStruct.ICParamsNames;

%--------------------------------------------------------------------------
% Set up DDFT IC parameters in potParamsDDFT
%--------------------------------------------------------------------------

if(optsStruct.anyDDFT)
     
    nDDFT=optsStruct.nDDFT;
        
    for iDDFT=1:nDDFT
    
        N=optsStruct.PhysArea{iDDFT}.N;
        
        N = prod(N); %#ok
               
        nVars=length(vars);
        varVals=cell(nVars,1);
        
        addVarText=['''ICrho'', ''' optsStruct.ICrho '''' ];
        
        % potParamsDDFT has var (varS replicated as column of N points for each species)
        % and varS (as input values)
        for iVar=1:nVars
            varVals{iVar}=optsStruct.([vars{iVar} 'S']);

            % if there's only one value, use a scalar, else replicate over
            % N points for each value (species)
                        
            if(length(varVals{iVar})==1)
                addVarText=cat(2,addVarText, [',''' vars{iVar} ''',varVals{' num2str(iVar) '}']);  
            else
                addVarText=cat(2,addVarText, [',''' vars{iVar} ''',repmat( varVals{' num2str(iVar) '}.'',N,1)']);        
            end
            
            addVarText=cat(2,addVarText, [',''' vars{iVar} 'S'',optsStruct.'  vars{iVar} 'S']);
        end

        ICParamsCmd=['optsStruct.ICParamsDDFT{iDDFT} = struct(' addVarText ');'];
        
        eval(ICParamsCmd)
    end

end
   

%--------------------------------------------------------------------------
% Set up plotting U parameters in potParamsPlot (only for V_G)
%--------------------------------------------------------------------------

if(optsStruct.dim==1)
    NPlot=optsStruct.NPlot^2; %#ok
elseif(optsStruct.dim==2)
    NPlot=optsStruct.NPlot1*optsStruct.NPlot2;  %#ok
end

nVars=length(vars);
varVals=cell(nVars,1);
addVarText=[];

% potParamsPlot has var (varS replicated as column of NPlot points for each species)
% and varS (as input values)
for iVar=1:nVars
    varVals{iVar}=optsStruct.([vars{iVar} 'S']).';           


    % if there's only one value, use a scalar, else replicate over
    % NPlot points for each value (species)
    
    if(iVar==1)
        comma = [];
    else
        comma = ',';
    end
    
    if(length(varVals{iVar})==1)
        addVarText=cat(2,addVarText, [comma '''' vars{iVar} ''',varVals{' num2str(iVar) '}']);  
    else
        addVarText=cat(2,addVarText, [comma '''' vars{iVar} ''',repmat( varVals{' num2str(iVar) '},NPlot,1)']);        
    end

    addVarText=cat(2,addVarText, [',''' vars{iVar} 'S'',optsStruct.'  vars{iVar} 'S']);
end


flowParamsCmd=['optsStruct.flowParamsPlot = struct(' addVarText ');'];

eval(flowParamsCmd)


end
