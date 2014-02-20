function [VAdd,dVAdd] = getVAdd(x1,x2,t,optsPhys)

% get name of function
getV1DV1=str2func(optsPhys.V1DV1);

% retrieve the part we're interested in
[h1,VAdd_S]=getV1DV1(x1,x2,t,optsPhys);
VAdd      = VAdd_S.V;
if(isfield(VAdd_S,'dy1') && isfield(VAdd_S,'dy2') && nargout == 2)
    dVAdd.dy1 = VAdd_S.dy1;
    dVAdd.dy2 = VAdd_S.dy2;
end
