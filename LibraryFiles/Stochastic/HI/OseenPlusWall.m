function D = OseenPlusWall(x,optsPhys)

D0=optsPhys.D0;
sigmaH=optsPhys.sigmaH;
dim=optsPhys.dim;

if(dim~=3)
    fprintf(1,'Wall only implemented for dimension 3!');
    pause;
end

D=diag(D0)*(eye(length(x))+makeOseen(x,sigmaH,dim) + wallT(x,sigmaH));


end
