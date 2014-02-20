function [VBack,DVBack] = getVBackDVBack1D(x,t,optsPhys)

% get name of function
getV1DV1=str2func(optsPhys.V1DV1);

% retrieve the part we're interested in
[VBack_S,~]=getV1DV1(x,t,optsPhys);
VBack      = VBack_S.V;

if(isfield(VBack_S,'DV') && nargout == 2)
    DVBack = VBack_S.DV;
end

end

