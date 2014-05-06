function optsStruct=checkPotentials(optsStruct)
% optsStruct=checkPotentials(optsStruct)
%   Check that V1DV1 potential and potNames are the correct size ({1,3})
%
% INPUTS:
%   optsStruct  -- struct of size [1 1] containing
%                   V1DV1:    name of 1-body potential ({string} or
%                              {string a,string b,string c})
%                   potNames: save name of potentials ({string} or
%                              {string a,string b,string c})
% OUTPUTS:
%   optsStruct  -- struct of size [1 1] with correctly formatted V1DV1 and
%                   potNames


if(length(optsStruct.V1DV1)==1)
    V1DV1=optsStruct.V1DV1{1};
    optsStruct.V1DV1={V1DV1,V1DV1,V1DV1};
end

if(length(optsStruct.potNames)==1)
    potNames=optsStruct.potNames{1};
    optsStruct.potNames={potNames,potNames,potNames};
end

end