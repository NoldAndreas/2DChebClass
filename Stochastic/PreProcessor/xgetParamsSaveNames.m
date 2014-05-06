function optsStruct=getParamsSaveNames(optsStruct)
% optsStruct=getParamsSaveNames(optsStruct)
%   Construct parameters for V1, V2 and HI for stochastic and DDFT
%   calculations from the given species-wise parameters
%
% INPUTS:
%   optsStruct  -- struct of size [1 1] containing
%                   potParamsNames:  names of V1 parameters {{'p1',...,'pn'}}
%                                     or {{'p1',...,'pn'},{...},{...}}
%                   potParams2Names: names of V2 parameters {{'p1',...,'pn'}}
%                   HIParamsNames:   names of HI parameters {{'p1',...,'pn'}}
%
% OUTPUTS:
%   optsStruct  -- struct of size [1 1] with all fields of input plus
%                   potParamsSaveNames:  names of save parameters for V1
%                                         {{'p1S',...,'pnS'},{...},{...}}
%                   potParams2SaveNames: names of save parameters for V1
%                                         {'p1S',...,'pnS'}
%                   HIParamsSaveNames:   names of save parameters for HI
%                                         {'p1S',...,'pnS'}

%--------------------------------------------------------------------------
% Construct save parameters for V1
%--------------------------------------------------------------------------

if(length(optsStruct.potParamsNames)==1)
    temp=optsStruct.potParamsNames{1};
    optsStruct.potParamsNames={temp, temp, temp};
end

for iGIF=1:3
    % remove extra braces
    vars=optsStruct.potParamsNames{iGIF};

    % add 'S' to the end of each parameter
    nVars=length(vars);
    varsS=cell(1,nVars);
    for iVars=1:nVars
        varsS{iVars}=[vars{iVars} 'S'];
    end
    
    optsStruct.potParamsSaveNames{iGIF}=varsS;
end

%--------------------------------------------------------------------------
% Construct save parameters for V2
%--------------------------------------------------------------------------

% remove extra braces
vars=optsStruct.potParams2Names{1};

% add 'S' to the end of each parameter
nVars=length(vars);
varsS=cell(1,nVars);
for iVars=1:nVars
    varsS{iVars}=[vars{iVars} 'S'];
end

optsStruct.potParams2SaveNames=varsS;

%--------------------------------------------------------------------------
% Construct save parameters for HI
%--------------------------------------------------------------------------

% remove extra braces
vars=optsStruct.HIParamsNames{1};

% add 'S' to the end of each parameter
nVars=length(vars);
varsS=cell(1,nVars);
for iVars=1:nVars
    varsS{iVars}=[vars{iVars} 'S'];
end

optsStruct.HIParamsSaveNames=varsS;

end