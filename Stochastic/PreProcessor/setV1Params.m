function optsStruct=setV1Params(optsStruct)
% optsStruct=setV1Params(optsStruct)
%   Replicate V1 parameter variables for stochastic calculations so
%   that they're the correct dimension (nParticles,1) and DDFT calculations
%   so they agree with the species-wise definitions
%
% INPUTS:
%   optsStruct  -- struct of size [1 1] containing
%                   potParamsNames:  names of V1 parameters {{'p1',...,'pn'}}
%                                     or {{'p1',...,'pn'},{...},{...}}
%                   <V1 Parameters>: called p1S, ..., pnS of form
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
%                     <V1 parameters>:    called p1S, ..., pnS of form
%                                          {(nSpecies,1),(nSpecies,1),(nSpecies,1)}
%                     IF anyStoc
%                       <V1 parameters>:  called p1, ..., pn, of form 
%                                          {pG,pI,pF} with pX (nParticles,1)
%                       potParams:        structure containing p1S, ..., pnS
%                                          and p1, ..., pn
%                     IF anyDDFT
%                       potParamsDDFT:    structure containing p1S, ..., pnS
%                                          and p1, ..., pn, = p1s, ..., pnS

%--------------------------------------------------------------------------
% Find list of V1 parameters
%--------------------------------------------------------------------------

vars=optsStruct.potParamsNames;


if(optsStruct.anyStoc)

    %----------------------------------------------------------------------
    % Duplicate V1 parameters for stochastic calculation
    %----------------------------------------------------------------------
    
    nParticles=optsStruct.nParticles;   %#ok
    nParticlesS=optsStruct.nParticlesS;
    nSpecies=optsStruct.nSpecies;
    
    endPos=cumsum(nParticlesS);
    startPos= [1;endPos(1:end-1)+1];
              
    for iVar=1:length(vars)
               
        % if only a scalar then store scalar
        if( length( optsStruct.([vars{iVar} 'S']) ) ==1)

            eval(['optsStruct.' vars{iVar} '=optsStruct.' vars{iVar} 'S;'])

        else
            % otherwise replicate by species

            % construct empty vector in each cell
            eval(['optsStruct.' vars{iVar} '=zeros(nParticles,1);'])       

            % assign values as per species
            for iSpecies=1:nSpecies
                mask=startPos(iSpecies) : endPos(iSpecies); %#ok
                eval(['optsStruct.' vars{iVar} '( mask ) = optsStruct.' vars{iVar} 'S(iSpecies);'])
            end

        end
    end
    
    %----------------------------------------------------------------------
    % Set up stochastic V1 parameters in potParams
    %----------------------------------------------------------------------

    % append each variable and its value to the command
    nVars=length(vars);
    addVarText=[];
    
    % var is (nParticles,1) as above, varS as in input
    for iVar=1:nVars
        % add var with value var
        
        if(iVar==1)
            comma = [];
        else
            comma = ',';
        end
        
        addVarText=cat(2,addVarText, [comma '''' vars{iVar} ''',optsStruct.'  vars{iVar}]);
        % add varS with value varS
        addVarText=cat(2,addVarText, [',''' vars{iVar} 'S'',optsStruct.'  vars{iVar} 'S']);
    end


    potParamsCmd=['optsStruct.potParams = struct(' addVarText ');'];
    eval(potParamsCmd)
    
end

%--------------------------------------------------------------------------
% Set up DDFT V1 parameters in potParamsDDFT
%--------------------------------------------------------------------------

if(optsStruct.anyDDFT)
     
    nDDFT=optsStruct.nDDFT;
        
    for iDDFT=1:nDDFT
    
        N=optsStruct.PhysArea{iDDFT}.N;
        
        N = prod(N); %#ok
               
        nVars=length(vars);
        varVals=cell(nVars,1);
        
        addVarText=['''V1DV1'', ''' optsStruct.V1DV1 '''' ];
        
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

        potParamsCmd=['optsStruct.potParamsDDFT{iDDFT} = struct(' addVarText ');'];
        
        eval(potParamsCmd)
    end

end
   

%--------------------------------------------------------------------------
% Set up plotting V1 parameters in potParamsPlot (only for V_G)
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


potParamsCmd=['optsStruct.potParamsPlot = struct(' addVarText ');'];

eval(potParamsCmd)


end
