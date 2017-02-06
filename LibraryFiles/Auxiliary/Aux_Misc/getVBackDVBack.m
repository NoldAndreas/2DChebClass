function [VBack,DVBack]=getVBackDVBack(x1,x2,optsPhys)

% get name of function
getV1DV1=str2func(optsPhys.V1DV1);

% retrieve the parts we're interested in
[VBack_S,~] = getV1DV1(x1,x2,0,optsPhys);

VBack  = VBack_S.V;
DVBack = VBack_S.grad;