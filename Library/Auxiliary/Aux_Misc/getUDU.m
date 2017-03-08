function [U,DU]=getUDU(x1,x2,t,optsPhys)

% get name of function
getUDU=str2func(optsPhys.UDU);

U_S = getUDU(x1,x2,t,optsPhys);

U  = U_S.U;
DU = U_S.div;