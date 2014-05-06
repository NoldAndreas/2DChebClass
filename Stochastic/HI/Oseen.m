function D= Oseen(x,optsPhys)

D0=optsPhys.D0;
sigmaH=optsPhys.sigmaH;
dim=optsPhys.dim;

D=diag(D0)*(eye(length(x))+makeOseen(x,sigmaH,dim));


end
