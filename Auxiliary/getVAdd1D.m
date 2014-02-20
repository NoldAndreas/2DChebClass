function VAdd=getVAdd1D(x,t,optsPhys)

% get name of function
getV1DV1=str2func(optsPhys.V1DV1);

% retrieve the part we're interested in
[h1,VAdd_S]=getV1DV1(x,t,optsPhys);
VAdd = VAdd_S.V;