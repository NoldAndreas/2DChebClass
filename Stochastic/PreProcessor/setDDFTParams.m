function optsStruct=setDDFTParams(optsStruct)
% optsStruct=setDDFTParams(optsStruct)
%   Assign DDFT parameters from DDFTParamsNames to DFTParams
%
% INPUTS:
%   optsStruct  -- struct of size [1 1] containing
%                   anyDDFT:          are we doing any DDFT calculations (t/f)
%                   DDFTParamsNames:  names of DDFT parameters {{'p1',...,'pn'}}
%                                      or {{'p1',...,'pn'}, ..., {...}} (of
%                                      length nDDFT)
%
% OUTPUTS:
%   optsStruct  -- struct of size [1 1] with all fields of input plus
%                     IF anyDDFT
%                       DDFTParams:   structure containing p1, ..., pn

if(optsStruct.anyDDFT)

    % find unique variable names over all DDFT calculations
    
    vars=optsStruct.DDFTParamsNames;
    nVarLists=length(vars);
    allVars=[];
    
    for iVarList=1:nVarLists
        allVars=[allVars,vars{iVarList}]; %#ok
    end
    
    uniqueVars=unique(allVars);
    
    % append each variable and its value to the command
    nVars=length(uniqueVars);
    addVarText=[];
    
    for iVar=1:nVars
        if(iVar>1)
            addVarText=cat(2,addVarText, ',');    
        end
        % add var with value varS
        addVarText=cat(2,addVarText, ['''' uniqueVars{iVar} ''',optsStruct.' uniqueVars{iVar}]);
    end

    % construct DDFTParams
    DDFTParamsCmd=['optsStruct.DDFTParams = struct(' addVarText ');'];
  
    eval(DDFTParamsCmd);
       
end

end